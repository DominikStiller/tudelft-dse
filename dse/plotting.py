import glob
import os
from pathlib import Path
from typing import Union

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib.font_manager import fontManager

# Load fonts
for f in glob.glob(os.path.abspath(os.path.dirname(__file__)) + "/../fonts/*.ttf"):
    fontManager.addfont(f)

assert (
    "Bookman Old Style" in fontManager.get_font_names()
), "Could not find plot font. Contact Dominik."

sb.set(
    context="paper",
    style="ticks",
    font="serif",
    rc={
        "lines.linewidth": 1.2,
        "axes.titleweight": "bold",
        "font.serif": "Bookman Old Style",
        "mathtext.fontset": "custom",
        "mathtext.it": "Cambria Math:italic",
        "mathtext.cal": "Cambria Math",
        "mathtext.rm": "Cambria Math",
        "figure.dpi": 300,
    },
)


def save_plot(results_folder: Union[Path, str], name: str, fig=None, type="pdf"):
    if isinstance(results_folder, str):
        results_folder = Path(results_folder)

    plots_folder = "plots" / results_folder
    plots_folder.mkdir(parents=True, exist_ok=True)

    if fig is None:
        fig = plt.gcf()
    fig.savefig(
        os.path.join(plots_folder, f"{name}.{type}"),
        dpi=450,
        bbox_inches="tight",
        pad_inches=0.01,
    )


def format_plot(
    xlocator=None,
    ylocator=None,
    tight_layout=True,
    zeroline=False,
    xgrid=True,
    ygrid=True,
):
    fig = plt.gcf()
    for ax in fig.axes:
        if zeroline:
            ax.axhline(0, linewidth=1.5, c="black")

        xlocator_ax = xlocator
        if not xlocator_ax:
            if ax.get_xscale() == "log":
                xlocator_ax = matplotlib.ticker.LogLocator(base=10, subs="auto", numticks=100)
            else:
                xlocator_ax = matplotlib.ticker.AutoMinorLocator()

        ylocator_ax = ylocator
        if not ylocator_ax:
            if ax.get_yscale() == "log":
                ylocator_ax = matplotlib.ticker.LogLocator(base=10, subs="auto", numticks=100)
            else:
                ylocator_ax = matplotlib.ticker.AutoMinorLocator()

        ax.get_xaxis().set_minor_locator(xlocator_ax)
        ax.get_yaxis().set_minor_locator(ylocator_ax)

        if xgrid:
            ax.grid(visible=True, which="major", linewidth=1.0)
        if ygrid:
            ax.grid(visible=True, which="minor", linewidth=0.5, linestyle="-.")

    if tight_layout:
        fig.tight_layout(pad=0.1, h_pad=0.4, w_pad=0.4)
