"""Generate a small geometry fixture without subject data or native tools."""
from pathlib import Path

import nibabel as nib
import numpy as np

from .workflow import record, report, screen, write


def demo(output):
    root = Path(output).resolve()
    root.mkdir(parents=True)
    inputs = root / "inputs"
    inputs.mkdir()
    shape = (52, 44, 57)
    affine = np.eye(4)
    affine[:3, 3] = [-26, -22, -32]
    xyz = np.ogrid[-26:26, -22:22, -32:25]
    support = (xyz[0]/24)**2 + (xyz[1]/20)**2 + ((xyz[2]+2)/33)**2 <= 1

    def volume(name, array):
        path = inputs / f"{name}.nii.gz"
        nib.save(nib.Nifti1Image(array, affine), path)
        item = record(path)
        item["path"] = path.relative_to(root).as_posix()
        return item

    def sphere(center, radius):
        return sum((axis - value)**2 for axis, value in zip(xyz, center, strict=True)) <= radius**2

    b0 = support * (75 + 15 * np.cos(xyz[0] / 5) + 8 * np.cos(xyz[2] / 3))
    reference = volume("synthetic_b0", b0.astype(np.float32))
    support_file = volume("support", support.astype(np.uint8))
    ribbon = np.zeros(shape, dtype=bool)
    gates = {}
    runs = []
    rng = np.random.default_rng(17)
    for side, sign in (("L", -1), ("R", 1)):
        x = sign * 8
        centers = {"Tang_pons": (x, -3, -8), "Tang_medulla": (x, -3, -18),
                   "Navigator_mRt": (x, 7, -5), "Navigator_sMRt": (x, 7, -11),
                   "Boyne_GRN": (x, 7, -19)}
        for name, center in centers.items():
            gates[f"{name}_{side}"] = volume(f"{name}_{side}", sphere(center, 2).astype(np.uint8))
        plic = volume(f"PLIC_{side}", sphere((x, -3, 8), 1).astype(np.uint8))
        cp = volume(f"CP_{side}", sphere((x, -3, 1), 1.5).astype(np.uint8))
        for origin, y in (("M1", -3), ("PMv", 4), ("SMA", -8), ("preSMA", 0)):
            seeds = sphere((x, y, 15), 2)
            ribbon |= seeds
            seed_mask = volume(f"seed_{side}_{origin}", seeds.astype(np.uint8))
            for model in ("Tensor_Det", "ST_iFOD2"):
                run_id = f"{model}_{side}_{origin}"
                parent = inputs / f"{run_id}.tck"
                log = inputs / f"{run_id}.seeds.txt"
                tracks = []
                for index in range(24):
                    jx, jy = rng.uniform(-.24, .24, 2)
                    start = [x + jx, y + jy, 15.]
                    family = index % 5
                    if family in (0, 1):
                        shift = sign * (2.1 if family == 1 and model == "ST_iFOD2" else 0.)
                        path = [[x + shift, -3, 8], [x, -3, 1], [x, -3, -8],
                                [x, -3, -18], [x + jx, -3 + jy, -23]]
                    elif family == 2:
                        path = [[x, 7, -5], [x, 7, -11], [x, 7, -19], [x, 7, -23]]
                    elif family == 3:
                        path = [[x + sign * 8, 10, 3], [x + sign * 9, 10, -14]]
                    else:
                        path = [[x, y, 11]]
                    tracks.append(np.asarray([start, *path], dtype=np.float32))
                tractogram = nib.streamlines.Tractogram(tracks, affine_to_rasmm=np.eye(4))
                nib.streamlines.TckFile(tractogram).save(str(parent))
                with log.open("w", encoding="utf-8") as stream:
                    for i, line in enumerate(tracks):
                        stream.write(f"{i},0,{line[0,0]:.8f},{line[0,1]:.8f},{line[0,2]:.8f}\n")
                run = {"id": run_id, "origin": origin, "side": side, "model": model,
                       "seed_mask": seed_mask, "tracking_support": support_file, "attempts": 32,
                       "saved": len(tracks), "ribbon_distance_mm": 2.0}
                for key, path in (("parent", parent), ("seed_log", log)):
                    run[key] = {**record(path), "path": path.relative_to(root).as_posix()}
                if model == "ST_iFOD2":
                    run["plic_waypoints"] = {"PLIC": plic, "CP": cp, "pons": gates[f"Tang_pons_{side}"]}
                runs.append(run)
    config = root / "synthetic.local.json"
    write(config, {"schema": "neo-tractography-1", "synthetic": True, "reference": reference,
                   "gate_support": support_file, "ribbon": volume("ribbon", ribbon.astype(np.uint8)),
                   "gates": gates, "slab_mm": 3, "runs": runs})
    results = root / "screened"
    screen(config, results)
    html = root / "report.html"
    report(results, html)
    return {"config": str(config), "results": str(results), "report": str(html)}
