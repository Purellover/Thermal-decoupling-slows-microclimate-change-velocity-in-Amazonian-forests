#!/usr/bin/env python3
"""
Figure 3. Climate change velocity maps and hotspot rose diagrams.

Workflow
1. Read microclimate and macroclimate velocity rasters.
2. Read hotspot polygons.
3. Insert precomputed rose diagrams for each hotspot.
4. Export one figure for each emissions scenario.

This script assumes that the hotspot rose diagrams have already been generated as PDF files.
Update the paths in CONFIG before running.
"""

from pathlib import Path
import warnings

import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "outputs"

CONFIG = {
    "fig3_dir": DATA_DIR / "fig3",
    "shape_file": DATA_DIR / "shapefiles" / "amazon_bio.shp",
    "hotspot_file": DATA_DIR / "shapefiles" / "hotpoint.shp",
    "output_dir": OUTPUT_DIR / "fig3",
    "scenarios": ["ssp245", "ssp585"],
    "scenario_labels": {"ssp245": "SSP2-4.5", "ssp585": "SSP5-8.5"},
    "velocity_min": 0,
    "velocity_max": 10,
}

CMAP_VELOCITY = LinearSegmentedColormap.from_list(
    "velocity", ["#4575b4", "#74add1", "#abd9e9", "#e0f3f8",
                 "#ffffbf", "#fee090", "#fdae61", "#f46d43",
                 "#d73027", "#a50026"], N=256
)
HOTSPOT_COLORS = ["#2166ac", "#762a83", "#1b7837", "#b2182b"]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 30,
    "axes.linewidth": 1.0,
})


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def read_raster(path):
    with rasterio.open(path) as src:
        array = src.read(1).astype(np.float32)
        nodata = src.nodata
        extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
    if nodata is not None:
        array[array == nodata] = np.nan
    array[array <= -9990] = np.nan
    return array, extent


def pdf_to_array(pdf_path):
    """Convert the first page of a PDF rose diagram to a NumPy image array."""
    try:
        import fitz
        doc = fitz.open(str(pdf_path))
        page = doc[0]
        pix = page.get_pixmap(matrix=fitz.Matrix(2.0, 2.0))
        image = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
        doc.close()
        return image[:, :, :3] if image.shape[2] == 4 else image
    except Exception:
        return None


def draw_velocity_map(ax, array, extent, amazon_gdf, hotspot_gdf, label, title):
    im = ax.imshow(array, cmap=CMAP_VELOCITY,
                   vmin=CONFIG["velocity_min"], vmax=CONFIG["velocity_max"],
                   extent=extent, origin="upper", aspect="auto", interpolation="nearest")
    amazon_gdf.boundary.plot(ax=ax, color="#111111", linewidth=1.0)

    x_range = extent[1] - extent[0]
    y_range = extent[3] - extent[2]
    fig = ax.get_figure()
    bbox = ax.get_position()
    ax_width = fig.get_size_inches()[0] * bbox.width
    ax_height = fig.get_size_inches()[1] * bbox.height
    radius_screen_inches = 0.60
    radius_x = radius_screen_inches / (ax_width / x_range)
    radius_y = radius_screen_inches / (ax_height / y_range)

    for i, (_, row) in enumerate(hotspot_gdf.iterrows()):
        lon = row.geometry.centroid.x
        lat = row.geometry.centroid.y
        color = HOTSPOT_COLORS[i]

        ax.add_patch(Ellipse((lon, lat), 2 * radius_x, 2 * radius_y,
                             fill=False, edgecolor=color, linewidth=4.5,
                             linestyle="--", transform=ax.transData, zorder=5))
        ax.add_patch(Ellipse((lon, lat), 0.6 * radius_x, 0.6 * radius_y,
                             facecolor=color, edgecolor="white", linewidth=1.0,
                             transform=ax.transData, zorder=7))
        ax.text(lon, lat, str(i + 1), ha="center", va="center",
                fontsize=30, fontweight="bold", color="white", zorder=8)

    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.02, 0.97, label, transform=ax.transAxes,
            fontsize=42, fontweight="bold", va="top")
    ax.set_title(title, fontsize=33, fontweight="normal", pad=6)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return im


def draw_rose_row(axes, images, prefix):
    for i, (ax, image) in enumerate(zip(axes, images)):
        color = HOTSPOT_COLORS[i]
        if image is not None:
            ax.imshow(image, aspect="auto")
        else:
            ax.text(0.5, 0.5, f"Hotspot {i + 1}", ha="center", va="center",
                    transform=ax.transAxes, fontsize=18)
        for spine in ax.spines.values():
            spine.set_edgecolor(color)
            spine.set_linewidth(2.5)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{prefix}{i + 1}", fontsize=27, color=color, fontweight="normal", pad=5)


def make_figure(scenario):
    fig3_dir = CONFIG["fig3_dir"]
    label = CONFIG["scenario_labels"][scenario]

    micro_path = fig3_dir / scenario / "10km" / f"climate_velocity_{scenario}_10km.tif"
    macro_path = fig3_dir / "macro" / scenario / "10km" / f"climate_velocity_{scenario}_10km.tif"

    micro, micro_extent = read_raster(micro_path)
    macro, macro_extent = read_raster(macro_path)

    amazon_gdf = gpd.read_file(CONFIG["shape_file"]).to_crs("EPSG:4326")
    hotspot_gdf = gpd.read_file(CONFIG["hotspot_file"]).to_crs("EPSG:4326")

    rose_dir = fig3_dir / "rose_hotspot" / scenario
    rose_micro, rose_macro = [], []
    for i in range(1, 5):
        rose_micro.append(pdf_to_array(rose_dir / f"rose_micro_hotspot{i}_Hotspot_{i}_{scenario}.pdf"))
        rose_macro.append(pdf_to_array(rose_dir / f"rose_macro_hotspot{i}_Hotspot_{i}_{scenario}.pdf"))

    fig = plt.figure(figsize=(22, 22))
    gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.04,
                           top=0.97, bottom=0.05, left=0.04, right=0.97,
                           height_ratios=[1.0, 1.8, 1.0])

    gs_top = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[0], wspace=0.10)
    gs_map = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1], wspace=0.04)
    gs_bottom = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[2], wspace=0.10)

    axes_micro_rose = [fig.add_subplot(gs_top[0, i]) for i in range(4)]
    ax_micro = fig.add_subplot(gs_map[0, 0])
    ax_macro = fig.add_subplot(gs_map[0, 1])
    axes_macro_rose = [fig.add_subplot(gs_bottom[0, i]) for i in range(4)]

    draw_rose_row(axes_micro_rose, rose_micro, "a")
    draw_rose_row(axes_macro_rose, rose_macro, "b")

    draw_velocity_map(ax_micro, micro, micro_extent, amazon_gdf, hotspot_gdf,
                      "a", f"Microclimate — {label}")
    draw_velocity_map(ax_macro, macro, macro_extent, amazon_gdf, hotspot_gdf,
                      "b", f"Macroclimate — {label}")

    cbar_ax = fig.add_axes([0.18, 0.015, 0.64, 0.018])
    sm = plt.cm.ScalarMappable(
        cmap=CMAP_VELOCITY,
        norm=mcolors.Normalize(vmin=CONFIG["velocity_min"], vmax=CONFIG["velocity_max"]),
    )
    cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cb.set_label("Climate Change Velocity (km/yr)", fontsize=33, labelpad=6)
    cb.ax.tick_params(labelsize=26)
    cb.set_ticks(np.arange(CONFIG["velocity_min"], CONFIG["velocity_max"] + 1, 1))

    CONFIG["output_dir"].mkdir(parents=True, exist_ok=True)
    out_pdf = CONFIG["output_dir"] / f"Fig3_{scenario.upper()}.pdf"
    out_png = CONFIG["output_dir"] / f"Fig3_{scenario.upper()}.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_pdf}")


def main():
    for scenario in CONFIG["scenarios"]:
        make_figure(scenario)


if __name__ == "__main__":
    main()
