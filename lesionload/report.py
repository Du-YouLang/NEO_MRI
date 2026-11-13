"""Write an offline C–F viewer linked by RAS z."""

import base64
from html import escape
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image

def encode(data):
    return base64.b64encode(data).decode("ascii")


def ras_display(data):
    """Input is RAS x,y,z. Display is left-to-right, anterior on top."""
    return data.transpose(2, 1, 0)[:, ::-1, :]


def contours(volume):
    """Exact voxel-edge runs for every axial mask slice, including internal holes."""
    result = []
    for mask in ras_display(volume):
        if not mask.any():
            result.append("")
            continue
        horizontal = np.diff(np.pad(mask, ((1, 1), (0, 0))), axis=0)
        vertical = np.diff(np.pad(mask, ((0, 0), (1, 1))), axis=1).T
        segments = []
        for direction, boundaries in enumerate((horizontal, vertical)):
            for line in np.flatnonzero(boundaries.any(axis=1)):
                endpoints = np.flatnonzero(np.diff(np.pad(boundaries[line], 1)))
                for start, end in endpoints.reshape((-1, 2)):
                    segments.append((line, start, line, end) if direction else (start, line, end, line))
        result.append(encode(np.array(segments, dtype="<u2").tobytes()))
    return result


def background_tiles(background):
    positive = background[background > 0]
    ceiling = float(np.percentile(positive, 99)) if positive.size else 1.0
    raster = ras_display(np.rint(np.clip(background / ceiling, 0, 1) * 255).astype(np.uint8))
    nz, ny, nx = raster.shape
    columns = 14
    atlas = Image.new("LA", (columns * nx, ((nz + columns - 1) // columns) * ny))
    for k, layer in enumerate(raster):
        alpha = (layer > 0).astype(np.uint8) * 255
        atlas.paste(Image.fromarray(np.stack((layer, alpha), axis=-1)),
                    ((k % columns) * nx, (k // columns) * ny))
    stream = io.BytesIO()
    atlas.save(stream, format="PNG", optimize=True)
    return columns, encode(stream.getvalue())


def format_value(value, percent=False):
    if value is None:
        return "—"
    if isinstance(value, list):
        return "–".join(format_value(item, percent) for item in value)
    return f"{value:.2%}" if percent else f"{value:.4g}"


def render(path, method, affine, background, bank, lesion, covered, tables):
    if not np.isfinite(background).all():
        raise ValueError("Viewer background is not finite")
    sides = ("Left", "Right")
    origins = sorted({name.split("_", 1)[1] for name in bank if name.startswith(sides)})
    if not origins:
        origins = ["Union"]
    summary = {row["region"]: row for row in tables["whole"]}
    axial = {(row["region"], row["z_ras_mm"]): row for row in tables["axial"]}
    diffusion = {}
    for row in tables["fa_axial"]:
        if method == "Boyne2025" and row["group"] != "covered_all":
            continue
        key = row["region"], row["z_ras_mm"]
        numerator, count = diffusion.get(key, (0., 0))
        n = row["valid_fa_voxels"]
        diffusion[key] = numerator + (row["FA_mean"] or 0.) * n, count + n
    reverse = tuple(slice(None, None, int(np.sign(affine[i, i]))) for i in range(3))
    nx, ny, nz = lesion.shape
    lower = np.minimum(affine[:3, 3], affine[:3, 3] + np.diag(affine)[:3] * (np.array(lesion.shape) - 1))
    levels = lower[2] + np.arange(nz)
    regions = {}
    atlas = {}
    for side in sides:
        for origin in origins:
            name = f"{side}_{origin}"
            tract = bank.get(name)
            if origin == "Union":
                tract = bank["Union"]
                name = "Union"
            if tract is None:
                tract = np.zeros_like(lesion)
            atlas[f"{side}_{origin}"] = contours(tract[reverse])
            whole = summary.get(name, {})
            values = [whole.get(key, 0) for key in ("hits_voxels", "atlas_voxels", "covered_voxels")]
            layers = []
            for z in levels:
                row = axial.get((name, z))
                numerator, n = diffusion.get((name, z), (0., 0))
                layers.append(None if row is None else [
                    row["hits_voxels"], row["atlas_voxels"], row["covered_voxels"],
                    numerator / n if n else None, n,
                ])
            regions[f"{side}_{origin}"] = {"whole": values, "slices": layers}
    columns, image = background_tiles(background[reverse])
    payload = dict(shape=lesion.shape, zmin=float(lower[2]), x_edge=float(lower[0] - .5),
                   y_edge=float(lower[1] + ny - .5), columns=columns, background=image,
                   cases={"current": {"regions": regions, "lesion": contours(lesion[reverse]),
                                      "has_fa": bool(tables["fa_axial"])}},
                   atlas=atlas, origins=origins, method=method)
    headers = ["分束", "N", "C", "H", "U", "H/C", "全束比例区间", "轴位峰值区间"]
    if method == "HCPA":
        headers.append("Max-WLL 区间（体素等效）")
    rows = []
    for name, row in summary.items():
        cells = [escape(name)] + [str(row[key]) for key in (
            "atlas_voxels", "covered_voxels", "hits_voxels", "unknown_voxels")]
        cells += [format_value(row[key], percent=True) for key in (
            "fraction_covered", "fraction_bounds", "peak_fraction_bounds")]
        if method == "HCPA":
            cells.append(format_value(row["max_wll_bounds_voxel"]))
        rows.append("<tr>" + "".join(f"<td>{value}</td>" for value in cells) + "</tr>")
    table = ("<div class='scroll'><table><thead><tr>" + "".join(f"<th>{h}</th>" for h in headers)
             + "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>")
    template = Path(__file__).with_name("viewer.html").read_text(encoding="utf-8")
    content = template.replace("__METHOD__", escape(method)).replace("__SUMMARY__", table)
    content = content.replace("__DATA__", json.dumps(payload, separators=(",", ":"), allow_nan=False)
                              .replace("<", "\\u003c"))
    with path.open("x", encoding="utf-8") as stream:
        stream.write(content)
