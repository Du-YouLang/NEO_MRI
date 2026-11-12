"""Render a selected run and gate from verified saved results."""
from itertools import islice
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import nibabel as nib
from nibabel.affines import voxel_sizes
import numpy as np

from .geometry import checked, load_reference
from .workflow import indices_for, read, result_array, result_run, selected_lines

VIEWS = {"覆盖计数": "visits.npz", "逐层 15% 核心": "core_15pct.npz"}
LAYERS = ("原始父流线", "筛选后流线", *VIEWS)
M2 = "M2（SMA / preSMA）"


def density_panels(fig, ref, volume, label):
    image = np.asarray(ref.dataobj)
    spacing = voxel_sizes(ref.affine)
    axes = fig.subplots(1, 3)
    maximum = max(1, int(volume.max()))
    for axis, ax in enumerate(axes):
        other = [i for i in range(3) if i != axis]
        weights = np.count_nonzero(volume, axis=tuple(other))
        index = int(weights.argmax()) if weights.any() else ref.shape[axis] // 2
        slice_image = np.take(image, index, axis=axis).T
        slice_visits = np.take(volume, index, axis=axis).T
        aspect = spacing[other[1]] / spacing[other[0]]
        ax.imshow(slice_image, origin="lower", cmap="gray", aspect=aspect)
        overlay = ax.imshow(np.ma.masked_where(slice_visits == 0, slice_visits), origin="lower",
                            cmap="inferno", alpha=.8, aspect=aspect, vmin=0,
                            vmax=maximum, interpolation="nearest")
        ax.set_title(f"{('Sagittal', 'Coronal', 'Axial')[axis]} {index}")
        ax.set_axis_off()
    fig.colorbar(overlay, ax=axes, shrink=.67, label=label)
    return fig


def report_preview(ref, visits):
    return density_panels(Figure(figsize=(9, 3), layout="constrained"), ref, visits, "Visits")


def render(folder, run_id, gate, view):
    folder = Path(folder)
    result, run, root = result_run(folder, run_id, gate)
    ref = load_reference(checked(result["reference"], folder))
    stats = run["gates"][gate]
    if view in VIEWS:
        volume = result_array(root, run, VIEWS[view], gate)
        fig = density_panels(plt.figure(figsize=(11, 4), layout="constrained"), ref, volume,
                             "Core mask" if view == LAYERS[3] else "Visits")
    else:
        if view not in LAYERS:
            raise ValueError(f"Unknown display mode: {view}")
        parent = checked(run["parent"], folder)
        lines = (nib.streamlines.load(parent, lazy_load=True).streamlines if view == LAYERS[0]
                 else selected_lines(parent, indices_for(root, run, gate)))
        fig = plt.figure(figsize=(8, 5), layout="constrained")
        ax = fig.add_subplot(111, projection="3d")
        for line in islice(lines, 120):
            ax.plot(*line.T, alpha=.52, color="#166a80", linewidth=1)
        ax.set(xlabel="R (mm)", ylabel="A (mm)", zlabel="S (mm)")
        ax.set_box_aspect((1, 1, 1))
        ax.set_title("First 120 available streamlines at most")
    fig.suptitle(f"{run_id} | {gate}\n{stats['selected']} selected / {run['saved']} saved / {run['attempts']} attempts")
    return fig


def interactive(folder):
    """Cascading controls. Rendering reads artifacts only."""
    import ipywidgets as widgets
    from IPython.display import display

    folder = Path(folder)
    runs = read(folder / "results.json")["runs"]
    model = widgets.Dropdown(description="模型", options=sorted({r["model"] for r in runs}))
    origin = widgets.Dropdown(description="起点")
    side = widgets.Dropdown(description="侧别")
    run = widgets.Dropdown(description="运行")
    gate = widgets.Dropdown(description="筛选门")
    view = widgets.Dropdown(description="图层", options=LAYERS, value=LAYERS[1])
    show = widgets.Button(description="查看", button_style="primary")
    panel = widgets.Output()
    state = {"refreshing": False}

    def choices(widget, options):
        previous = widget.value
        widget.options = options
        widget.value = previous if previous in options else options[0]

    def refresh(_change=None):
        if state["refreshing"]:
            return
        state["refreshing"] = True
        try:
            matching = [r for r in runs if r["model"] == model.value]
            origins = sorted({r["origin"] for r in matching})
            m2_sides = sorted({
                r["side"] for r in matching if r["origin"] == "SMA"
            } & {r["side"] for r in matching if r["origin"] == "preSMA"})
            if m2_sides:
                origins.append(M2)
            choices(origin, origins)
            if origin.value == M2:
                matching = [r for r in matching if r["origin"] in {"SMA", "preSMA"} and r["side"] in m2_sides]
            else:
                matching = [r for r in matching if r["origin"] == origin.value]
            choices(side, sorted({r["side"] for r in matching}))
            matching = [r for r in matching if r["side"] == side.value]
            choices(run, [r["id"] for r in matching])
            selected = next(r for r in matching if r["id"] == run.value)
            choices(gate, list(selected["gates"]))
        finally:
            state["refreshing"] = False

    def draw(_button=None):
        with panel:
            panel.clear_output(wait=True)
            fig = render(folder, run.value, gate.value, view.value)
            display(fig)
            plt.close(fig)

    for control in (model, origin, side, run):
        control.observe(refresh, names="value")
    show.on_click(draw)
    refresh()
    display(widgets.VBox([widgets.HBox([model, origin, side]),
                          widgets.HBox([run, gate, view, show]), panel]))
    draw()
    return {"model": model, "origin": origin, "side": side, "run": run, "gate": gate,
            "view": view, "show": show}
