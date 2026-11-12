"""Screen ordered waypoints and accumulate visits from the same parent indices."""
from collections import Counter
import csv

import nibabel as nib
from nibabel.affines import apply_affine, voxel_sizes
import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.spatial import cKDTree

GATE_ROLES = ("Tang_pons", "Tang_medulla", "Navigator_mRt", "Navigator_sMRt", "Boyne_GRN")
ORDINARY = {
    "Tang_pons": ("Tang_pons",),
    "Tang_pons_medulla": ("Tang_pons", "Tang_medulla"),
    "Navigator_mRt": ("Navigator_mRt",),
    "Navigator_sMRt": ("Navigator_mRt", "Navigator_sMRt"),
    "Navigator_GRN": ("Navigator_mRt", "Boyne_GRN"),
}
WAYPOINTS = {
    "CST_PLIC_CP": ("PLIC", "CP"),
    "CST_PLIC_CP_Tang_pons": ("PLIC", "CP", "pons"),
    "CST_PLIC_plus2_CP": ("PLIC_plus2", "CP"),
    "CST_PLIC_plus2_CP_Tang_pons": ("PLIC_plus2", "CP", "pons"),
}


def sample_line(line):
    """Keep original vertices and subdivide longer segments to at most 0.25 mm."""
    delta = np.diff(line, axis=0)
    lengths = np.linalg.norm(delta, axis=1)
    if lengths.max() <= .25:
        return line
    parts = np.maximum(1, np.ceil(lengths / .25).astype(int))
    ends = parts.cumsum()
    seg = np.repeat(np.arange(len(parts)), parts)
    offsets = np.arange(ends[-1]) - np.repeat(ends - parts, parts)
    return np.vstack((line[seg] + delta[seg] * (offsets / parts[seg])[:, None], line[-1]))


def ordered(positions, sequence):
    cursor = -1
    for role in sequence:
        hits = positions[role]
        index = hits.searchsorted(cursor, side="right")
        if index == len(hits):
            return False
        cursor = hits[index]
    return True


def read_tck_header(path):
    header = {}
    with path.open("rb") as stream:
        if stream.readline(100).rstrip() != b"mrtrix tracks":
            raise ValueError(f"Not an MRtrix TCK: {path}")
        while stream.tell() < 1024 * 1024:
            raw = stream.readline(16384)
            if raw.rstrip() == b"END":
                location = header.get("file", "").split()
                if len(location) != 2 or location[0] != "." or int(location[1]) < stream.tell():
                    raise ValueError("TCK must contain its own binary data")
                return header
            if not raw or not raw.endswith(b"\n"):
                break
            key, sep, value = raw.partition(b":")
            if sep:
                header[key.decode("ascii")] = value.decode("utf-8").strip()
    raise ValueError(f"Invalid TCK header: {path}")


def seed_records(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = (line for line in stream if line.strip() and not line.lstrip().startswith("#"))
        for row in csv.reader(rows):
            yield int(row[0]), np.array([float(row[i]) for i in (2, 3, 4)])


def slab_indices(shape, affine, width):
    i, j, k = np.ogrid[:shape[0], :shape[1], :shape[2]]
    z = affine[2, 0]*i + affine[2, 1]*j + affine[2, 2]*k + affine[2, 3]
    slab = np.floor(z / width + .5).astype(np.int64)
    return slab - slab.min()


def core_mask(visits, slabs):
    active = visits > 0
    bins = slabs[active]
    values = visits[active].astype(np.int64)
    peak = np.zeros(int(slabs.max()) + 1, dtype=np.int64)
    np.maximum.at(peak, bins, values)
    core = np.zeros_like(active)
    core[active] = 100 * values >= 15 * peak[bins]
    return core


class Screening:
    def __init__(self, ref, ribbon, seeds, support, gate_support, masks, laterality,
                 plic=None):
        self.shape = np.asarray(ref.shape)
        self.inverse = np.linalg.inv(ref.affine)
        spacing = voxel_sizes(ref.affine)
        if not ribbon.any() or not seeds.any() or not support.any():
            raise ValueError("Ribbon, seed mask and tracking support must be nonempty")
        self.seed_near = distance_transform_edt(~seeds, sampling=spacing) <= 2
        self.ribbon = cKDTree(apply_affine(ref.affine, np.argwhere(ribbon)))
        self.support = support
        self.laterality = laterality
        gate_names = [f"{role}_{side}" for role in GATE_ROLES for side in "LR"]
        self.bits = {name: 1 << i for i, name in enumerate(gate_names)}
        self.regular = np.zeros(ref.shape, dtype=np.uint16)
        self.gate_counts = {}
        for name, bit in self.bits.items():
            native = masks[name]
            used = native & gate_support
            self.regular[used] |= bit
            self.gate_counts[name] = {"native": int(native.sum()), "gate_support": int(used.sum()),
                                      "tracking_support": int((native & support).sum())}
        self.plic = None
        if plic is not None:
            if not all(plic[name].any() for name in ("PLIC", "CP", "pons")):
                raise ValueError("The PLIC profile requires nonempty PLIC, CP and original pons")
            expanded = (distance_transform_edt(~plic["PLIC"], sampling=spacing) <= 2) & support
            self.plic = np.zeros(ref.shape, dtype=np.uint8)
            for i, mask in enumerate((plic["PLIC"], expanded, plic["CP"], plic["pons"])):
                self.plic[mask] |= 1 << i
            self.gate_counts["PLIC_plus2"] = {"tracking_support": int(expanded.sum())}

    def classify(self, line, side, ribbon_mm, counts):
        start = np.rint(apply_affine(self.inverse, line[0])).astype(int)
        if ((start < 0) | (start >= self.shape)).any() or not self.seed_near[tuple(start)]:
            return (), None
        counts["seed_near"] += 1
        if self.ribbon.query(line[0])[0] > ribbon_mm:
            return (), None
        counts["cortical_start"] += 1
        if line[0, 2] - line[-1, 2] < 10:
            return (), None
        counts["descending"] += 1
        indices = np.floor(apply_affine(self.inverse, sample_line(line)) + .5).astype(int)
        inside = ((indices >= 0) & (indices < self.shape)).all(axis=1)
        ijk = tuple(indices[inside].T)
        selected = []
        flags = self.regular[ijk]
        if flags.any():
            anatomical_side = self.laterality(line[0])
            counts["seed_vs_template_side_mismatch"] += int(anatomical_side != side)
            positions = {}
            for role in GATE_ROLES:
                sides = anatomical_side if role.startswith("Tang") else "LR"
                bits = sum(self.bits[f"{role}_{s}"] for s in sides)
                positions[role] = np.flatnonzero(flags & bits)
            selected.extend(name for name, sequence in ORDINARY.items() if ordered(positions, sequence))
        if self.plic is not None:
            flags = self.plic[ijk]
            if flags.any():
                roles = ("PLIC", "PLIC_plus2", "CP", "pons")
                positions = {name: np.flatnonzero(flags & (1 << i)) for i, name in enumerate(roles)}
                selected.extend(name for name, sequence in WAYPOINTS.items() if ordered(positions, sequence))
        if not selected:
            return (), None
        occupied = np.unique(np.ravel_multi_index(ijk, tuple(self.shape))[self.support[ijk]])
        return selected, occupied

    def replay(self, parent, seed_log, side, ribbon_mm, saved):
        names = (*ORDINARY, *WAYPOINTS) if self.plic is not None else tuple(ORDINARY)
        shape = tuple(self.shape)
        visits = {name: np.zeros(int(self.shape.prod()), dtype=np.uint32) for name in names}
        indices = {name: [] for name in names}
        header = read_tck_header(parent)
        if int(header.get("count", -1)) != saved:
            raise ValueError("Parent TCK header count differs from receipt")
        tracks = nib.streamlines.load(parent, lazy_load=True).streamlines
        counts = Counter(dict.fromkeys(("seed_near", "cortical_start", "descending",
                                        "seed_vs_template_side_mismatch"), 0))
        total = 0
        max_error = 0.
        for index, (line, (seed_index, seed)) in enumerate(zip(tracks, seed_records(seed_log), strict=True)):
            total = index + 1
            if total > saved or len(line) < 2 or not np.isfinite(line).all():
                raise ValueError(f"Invalid parent streamline {index}")
            error = float(np.linalg.norm(line[0] - seed))
            if seed_index != index or not np.isfinite(error) or error > .05:
                raise ValueError(f"TCK and seed log disagree at streamline {index}")
            max_error = max(max_error, error)
            gates, occupied = self.classify(line, side, ribbon_mm, counts)
            for gate in gates:
                indices[gate].append(index)
                visits[gate][occupied] += 1
        if total != saved:
            raise ValueError("Parent TCK count differs from the receipt")
        if self.plic is not None:
            if (not set(indices["CST_PLIC_CP"]).issubset(indices["CST_PLIC_plus2_CP"])
                    and not np.any((self.plic & 1 > 0) & (self.plic & 2 == 0))):
                raise ValueError("PLIC expansion membership inconsistent with original support")
        return dict(counts, max_seed_error_mm=max_error), {
            key: np.asarray(value, dtype=np.int64) for key, value in indices.items()
        }, {key: value.reshape(shape) for key, value in visits.items()}
