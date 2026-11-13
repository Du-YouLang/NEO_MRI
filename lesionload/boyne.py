"""Reproduce the archived Boyne Method C subject-count maps from pinned assets."""

from collections import Counter
import gzip
import struct

import nibabel as nib
import numpy as np

from anatomy_review import file_record

ORIGINS = ("Whole", "M1", "PMd", "SMA", "preSMA", "ACC", "mPFC", "aPFC",
           "FEFarea", "Other", "Unassigned")
SIDES = {"L": "Left", "R": "Right"}
EXPECTED = {"L": 222374, "R": 208408}
AFFINE = np.array([[-1, 0, 0, 90], [0, 1, 0, -126], [0, 0, 1, -72], [0, 0, 0, 1]])
SHAPE = (182, 218, 182)
STRIDES = np.array([SHAPE[1] * SHAPE[2], SHAPE[2], 1])


def traversed_cells(first, last, starts, ends):
    """Yield intermediate cells; endpoints are assigned by the caller.

    Exact crossing order in voxel-center coordinates, with simultaneous
    crossings advanced together. All vertices must lie inside the grid.
    """
    gap = np.abs(ends - starts).sum(axis=1)
    rows = np.flatnonzero(gap > 1)
    if not len(rows):
        return
    start_point = first[rows].astype(np.float64)
    delta = last[rows].astype(np.float64) - start_point
    cell = starts[rows].copy()
    end = ends[rows]
    step = np.sign(end - cell)
    with np.errstate(divide="ignore", invalid="ignore"):
        boundary = (cell + .5 * step - start_point) / delta
        increment = np.abs(1.0 / delta)
    for _ in range(int(gap[rows].max()) - 1):
        pending = (cell != end)
        crossing = np.where(pending, boundary, np.inf)
        soonest = crossing.min(axis=1)
        move = pending & (np.abs(np.where(pending, crossing, 2.) - soonest[:, None]) <= 1e-12)
        cell += move * step
        boundary += np.where(move, increment, 0.)
        live = (cell != end).any(axis=1)
        if not live.any():
            return
        yield rows[live], cell[live]


def index_cells(points):
    return np.floor(points + .5).astype(np.int32)


def set_visits(visited, points, cells, ends, bits):
    offsets = np.r_[0, ends[:-1] + 1]
    counts = ends - offsets + 1
    line_bits = np.repeat(bits, counts)
    np.bitwise_or.at(visited, cells @ STRIDES, line_bits)
    active = np.ones(len(points), dtype=bool)
    active[ends] = False
    rows = np.flatnonzero(active)
    for which, intermediate in traversed_cells(points[rows], points[rows + 1],
                                                cells[rows], cells[rows + 1]):
        np.bitwise_or.at(visited, intermediate @ STRIDES, line_bits[rows[which]])


def roi_lookup(paths):
    lookup = np.zeros(SHAPE, dtype=np.uint32)
    sources = {side: [] for side in SIDES}
    for number, origin in enumerate(ORIGINS[1:-1]):
        for offset, side in enumerate(SIDES):
            name = f"{origin}_{side}.inDil2.nii.gz"
            path = paths[name]
            image = nib.load(path)
            data = np.asanyarray(image.dataobj)
            if data.shape == SHAPE + (1,):
                data = data[..., 0]
            if data.shape != SHAPE or not np.array_equal(image.affine, AFFINE) or not np.isin(data, (0, 1)).all():
                raise ValueError(f"Invalid cortical ROI geometry or values: {path}")
            lookup[data > 0] |= 1 << (2 * number + offset)
            sources[side].append(file_record(path))
    return lookup, sources

def voxelize(side, path, roi):
    atlas = np.zeros((len(ORIGINS), int(np.prod(SHAPE))), dtype=np.uint8)
    visited = np.zeros(atlas.shape[1], dtype=np.uint16)
    counts = Counter()
    patterns = Counter()
    endpoint_patterns = Counter()
    group_lines = Counter()
    members = Counter()
    chunks = []
    endpoint_limits = {end: [float("inf"), -float("inf")] for end in ("superior", "inferior")}
    label = None
    nlines = 0
    unassigned = 1 << (len(ORIGINS) - 2)
    side_offset = list(SIDES).index(side)
    origin_bits = [1 << (2 * i + side_offset) for i in range(len(ORIGINS) - 2)]
    assignment = np.full(1 << (2 * len(origin_bits)), unassigned, dtype=np.uint16)
    for i, flag in enumerate(origin_bits):
        assignment[flag] = 1 << i
    side_mask = sum(origin_bits)

    def flush_chunk():
        if not chunks:
            return
        lengths = np.fromiter((len(points) for points in chunks), np.int32)
        starts = np.r_[0, lengths.cumsum()[:-1]]
        ends = starts + lengths - 1
        points = np.concatenate(chunks)
        # DSI Studio writes voxel-index times spacing with its LPS storage order.
        # The declared affine is center-based; y is reversed into LAS NIfTI.
        points[:, 1] = SHAPE[1] - 1 - points[:, 1]
        if not np.isfinite(points).all():
            raise ValueError("Nonfinite streamline coordinates")
        upper_first = points[starts, 2] > points[ends, 2]
        upper = points[np.where(upper_first, starts, ends)]
        lower = points[np.where(upper_first, ends, starts)]
        if np.any(upper[:, 2] <= lower[:, 2]):
            raise ValueError(f"Indeterminate superior endpoint in {side}")
        for name, endpoints in (("superior", upper), ("inferior", lower)):
            z = endpoints[:, 2] + AFFINE[2, 3]
            endpoint_limits[name][0] = min(endpoint_limits[name][0], float(z.min()))
            endpoint_limits[name][1] = max(endpoint_limits[name][1], float(z.max()))
        cells = index_cells(points)
        if np.any(cells < 0) or np.any(cells >= SHAPE):
            raise ValueError("A CRST point falls outside its declared template lattice")
        vertex_labels = roi[tuple(cells.T)]
        endpoint_labels = vertex_labels[np.where(upper_first, starts, ends)] & side_mask
        endpoint_patterns.update(endpoint_labels.tolist())
        # The author's ROI/ROA filtering tests saved vertices. The earliest
        # label hit in their region-major L/R order wins after exclusions.
        membership = np.bitwise_or.reduceat(vertex_labels, starts)
        patterns.update(membership.tolist())
        first = membership & (-membership)
        assigned = assignment[first]
        assigned_list = assigned.tolist()
        counts.update(assigned_list)
        for bit in set(assigned_list):
            members[bit] |= 1 << label
        set_visits(visited, points, cells, ends, assigned)
        chunks.clear()

    def flush_subject():
        flush_chunk()
        atlas[0] += visited != 0
        for index in range(1, len(ORIGINS)):
            atlas[index] += (visited & (1 << (index - 1))) != 0
        visited.fill(0)

    with gzip.open(path, "rb") as source:
        header = source.read(1000)
        if len(header) != 1000:
            raise ValueError("Truncated TRK header")
        dims = struct.unpack_from("<3h", header, 6)
        mm = struct.unpack_from("<3f", header, 12)
        fmt = (struct.unpack_from("<h", header, 36)[0], struct.unpack_from("<h", header, 238)[0],
               struct.unpack_from("<i", header, 996)[0], header[948:951])
        actual_count = struct.unpack_from("<i", header, 988)[0]
        grid = np.asarray(struct.unpack_from("<16f", header, 440)).reshape(4, 4)
        flip = np.eye(4)
        flip[1, 1], flip[1, 3] = -1, SHAPE[1] - 1
        if (header[:6] != b"TRACK\x00" or dims != SHAPE or mm != (1., 1., 1.)
                or fmt != (0, 1, 1000, b"LPS") or actual_count != EXPECTED[side]
                or not np.array_equal(grid, AFFINE @ flip)):
            raise ValueError(f"Unexpected DSI track header: {path}")
        for _ in range(actual_count):
            prefix = source.read(4)
            if len(prefix) != 4:
                raise ValueError("Truncated streamline length")
            length = struct.unpack("<i", prefix)[0]
            if not 2 <= length <= 100000:
                raise ValueError("Invalid streamline point count")
            record = source.read(12 * length + 4)
            if len(record) != 12 * length + 4:
                raise ValueError("Truncated streamline record")
            group = struct.unpack_from("<f", record, 12 * length)[0]
            if not group.is_integer() or not 0 <= group < 50:
                raise ValueError("Invalid DSI subject cluster")
            group = int(group)
            if group != label:
                if label is not None:
                    flush_subject()
                if group in group_lines:
                    raise ValueError("Source has noncontiguous subject clusters")
                label = group
            group_lines[group] += 1
            nlines += 1
            chunks.append(np.frombuffer(record, dtype="<f4", count=3 * length).reshape(-1, 3))
            if len(chunks) == 640:
                flush_chunk()
            if nlines % 50000 == 0:
                print(f"{side}: {nlines:,}/{actual_count:,} streamlines", flush=True)
        flush_subject()
        if source.read(1) or set(group_lines) != set(range(50)) or sum(counts.values()) != actual_count:
            raise ValueError(f"Incomplete cluster or streamline accounting: {side}")
    if not np.array_equal(atlas[0] > 0, np.any(atlas[1:] > 0, axis=0)):
        raise ValueError("The whole-track support differs from the union of branch assignments")
    if np.any(atlas[1:] > atlas[0]):
        raise ValueError("A branch has higher subject support than its whole tract")
    origin_counts = {origin: dict(streamlines=counts[1 << i],
                                  contributing_subjects=members[1 << i].bit_count())
                     for i, origin in enumerate(ORIGINS[1:])}
    origin_counts["Whole"] = dict(streamlines=actual_count, contributing_subjects=50)
    summary = dict(**origin_counts, no_available_roi_hit=patterns[0],
                   multiple_roi_hits=sum(n for bit, n in patterns.items() if bit.bit_count() > 1),
                   excluded_by_earlier_opposite_roi=counts[unassigned] - patterns[0],
                   endpoint_only_no_match=endpoint_patterns[0],
                   endpoint_only_multiple_matches=sum(n for bit, n in endpoint_patterns.items() if bit.bit_count() > 1),
                   endpoint_only_unique_by_region={origin: endpoint_patterns[flag]
                                                   for origin, flag in zip(ORIGINS[1:-1], origin_bits, strict=True)},
                   endpoint_ras_z_ranges_mm=endpoint_limits,
                   cluster_counts=dict(sorted(group_lines.items())), trk_header_affine=grid.tolist(),
                   trk_to_las_voxel=flip.tolist())
    return atlas.reshape((len(ORIGINS), *SHAPE)), summary
