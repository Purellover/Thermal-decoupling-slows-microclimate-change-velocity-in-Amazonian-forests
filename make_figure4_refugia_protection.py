#!/usr/bin/env python3
"""
Figure 4. Microclimate refugia, climate residence time and protected-area coverage.

Workflow
1. Read protected-area polygons, refugia polygons and microclimate velocity rasters.
2. Compute climate residence time for each protected area as equivalent diameter divided by median velocity.
3. Calculate the overlap between refugia and protected areas.
4. Plot refugia maps, residence-time distributions and refugia area/patch counts.
5. Export a supplementary coverage-gap figure.

Update the paths in CONFIG before running.
"""

from pathlib import Path
import warnings

import geopandas as gpd
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask as rasterio_mask
from shapely.geometry import mapping
from shapely.ops import unary_union

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "outputs"

CONFIG = {
    "fig3_dir": DATA_DIR / "fig3",
    "fig4_dir": DATA_DIR / "fig4",
    "pa_dir": DATA_DIR / "fig4" / "protected_areas",
    "output_dir": OUTPUT_DIR / "fig4",
    "amazon_shape": DATA_DIR / "shapefiles" / "amazon_bio.shp",
    "projected_crs": "EPSG:5880",
    "scenarios": ["ssp245", "ssp585"],
    "scenario_labels": {"ssp245": "SSP2-4.5", "ssp585": "SSP5-8.5"},
}

COLORS = {
    "ssp245": "#4575b4",
    "ssp585": "#d73027",
    "refugia": "#15803d",
    "pa_fill": "#fef3c7",
    "pa_edge": "#ca8a04",
    "covered": "#74add1",
    "uncovered": "#f46d43",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 26,
    "axes.titlesize": 30,
    "axes.labelsize": 27,
    "xtick.labelsize": 27,
    "ytick.labelsize": 27,
    "legend.fontsize": 24,
    "axes.linewidth": 0.8,
})


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def format_int(value):
    return f"{int(value):,}".replace(",", "\u202f")


def read_protected_areas():
    files = [
        CONFIG["pa_dir"] / "ANP_with_diameter1.shp",
        CONFIG["pa_dir"] / "ANP_with_diameter2.shp",
        CONFIG["pa_dir"] / "ANP_with_diameter3.shp",
    ]
    layers = [gpd.read_file(path) for path in files]
    pa = gpd.GeoDataFrame(pd.concat(layers, ignore_index=True),
                          geometry="geometry", crs=layers[0].crs)
    pa_projected = pa.to_crs(CONFIG["projected_crs"])
    pa_projected["geometry"] = pa_projected.buffer(0)
    pa_projected["area_km2"] = pa_projected.geometry.area / 1e6
    pa_projected["diameter_km"] = 2 * np.sqrt(pa_projected["area_km2"] / np.pi)
    return pa, pa_projected


def read_refugia():
    refugia = {}
    for scenario in CONFIG["scenarios"]:
        path = CONFIG["fig4_dir"] / scenario / f"core_{scenario}.shp"
        gdf = gpd.read_file(path).to_crs(CONFIG["projected_crs"])
        gdf["geometry"] = gdf.buffer(0)
        refugia[scenario] = gdf
    return refugia


def velocity_path(scenario):
    return CONFIG["fig3_dir"] / scenario / "10km" / f"climate_velocity_{scenario}_10km.tif"


def compute_residence_time(pa_wgs84, pa_projected, scenario):
    """Climate residence time is equivalent diameter divided by median velocity."""
    with rasterio.open(velocity_path(scenario)) as src:
        velocity_crs = src.crs
        pa_for_mask = pa_wgs84.to_crs(velocity_crs)
        median_velocity = []

        for _, row in pa_for_mask.iterrows():
            try:
                array, _ = rasterio_mask(src, [mapping(row.geometry)], crop=True, nodata=np.nan)
                values = array[0].astype(np.float32)
                values[values <= 0] = np.nan
                values[values > 1e6] = np.nan
                median_velocity.append(np.nanmedian(values))
            except Exception:
                median_velocity.append(np.nan)

    result = pa_projected.copy()
    result["median_velocity"] = median_velocity
    result["residence_time"] = result["diameter_km"] / result["median_velocity"]
    result["scenario"] = scenario

    result = result[
        result["residence_time"].notna() &
        (result["residence_time"] > 0) &
        (result["residence_time"] < 10000)
    ].copy()
    return result


def compute_gap(refugia_gdf, pa_projected):
    pa_union = unary_union(pa_projected.geometry)
    gdf = refugia_gdf.copy()
    gdf["refugia_area_km2"] = gdf.geometry.area / 1e6

    covered = []
    for geom in gdf.geometry:
        try:
            covered.append(geom.intersection(pa_union).area / 1e6)
        except Exception:
            covered.append(0.0)

    gdf["covered_area_km2"] = covered
    gdf["uncovered_area_km2"] = (gdf["refugia_area_km2"] - gdf["covered_area_km2"]).clip(lower=0)
    gdf["coverage_ratio"] = (gdf["covered_area_km2"] / gdf["refugia_area_km2"]).clip(0, 1)

    total = gdf["refugia_area_km2"].sum()
    covered_total = gdf["covered_area_km2"].sum()

    return {
        "gdf": gdf,
        "total_area": total,
        "covered_area": covered_total,
        "uncovered_area": total - covered_total,
        "covered_pct": covered_total / total * 100,
    }


# ---------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------

def plot_main_figure(pa_wgs84, refugia, residence, gap):
    amazon_gdf = gpd.read_file(CONFIG["amazon_shape"]).to_crs("EPSG:4326")

    area_245 = gap["ssp245"]["total_area"] / 1e6
    area_585 = gap["ssp585"]["total_area"] / 1e6
    patch_245 = len(refugia["ssp245"])
    patch_585 = len(refugia["ssp585"])

    fig = plt.figure(figsize=(20, 16))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.30,
                           top=0.94, bottom=0.06, left=0.07, right=0.97,
                           height_ratios=[1.45, 1.0])

    # Panels a and b
    for col, scenario in enumerate(CONFIG["scenarios"]):
        ax = fig.add_subplot(gs[0, col])
        pa_wgs84.to_crs("EPSG:4326").plot(
            ax=ax, color=COLORS["pa_fill"], alpha=0.60,
            edgecolor=COLORS["pa_edge"], linewidth=0.25,
        )
        refugia[scenario].to_crs("EPSG:4326").plot(
            ax=ax, color=COLORS["refugia"], alpha=0.85, edgecolor="none",
        )
        amazon_gdf.plot(ax=ax, color="none", edgecolor="#111111", linewidth=1.0)

        ax.set_xlabel("Longitude (°)", fontsize=24)
        ax.set_ylabel("Latitude (°)" if col == 0 else "", fontsize=24)
        ax.tick_params(labelsize=22, direction="in", length=3)
        for spine in ax.spines.values():
            spine.set_visible(False)

        ax.text(0.02, 0.98, ["a", "b"][col], transform=ax.transAxes,
                fontsize=27, fontweight="bold", va="top")
        ax.set_title(f"Core Refugia — {CONFIG['scenario_labels'][scenario]}",
                     fontsize=26, fontweight="bold", pad=6)

        ax.text(0.03, 0.04,
                f"n = {format_int(len(refugia[scenario]))} patches\n"
                f"{gap[scenario]['total_area'] / 1e6:.2f} × 10⁶ km²",
                transform=ax.transAxes, fontsize=21, va="bottom", ha="left",
                bbox=dict(boxstyle="round,pad=0.35", fc="white",
                          ec="#cccccc", alpha=0.88, linewidth=0.6))

        if col == 0:
            legend = [
                mpatches.Patch(facecolor=COLORS["pa_fill"], edgecolor=COLORS["pa_edge"],
                               linewidth=0.8, label="Protected Areas"),
                mpatches.Patch(color=COLORS["refugia"], alpha=0.85, label="Core Refugia"),
                mpatches.Patch(color="none", ec="#111111", linewidth=1.0, label="Amazon boundary"),
            ]
            ax.legend(handles=legend, fontsize=21, loc="upper right",
                      framealpha=0.88, edgecolor="#cccccc", frameon=True)

    # Panel c
    ax_c = fig.add_subplot(gs[1, 0])
    data = [
        residence["ssp245"]["residence_time"].dropna().values,
        residence["ssp585"]["residence_time"].dropna().values,
    ]
    labels = [CONFIG["scenario_labels"]["ssp245"], CONFIG["scenario_labels"]["ssp585"]]
    colors = [COLORS["ssp245"], COLORS["ssp585"]]

    bp = ax_c.boxplot(data, patch_artist=True, showfliers=False, widths=0.5,
                      medianprops=dict(color="black", linewidth=2))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)

    for i, (values, color) in enumerate(zip(data, colors), start=1):
        median = np.nanmedian(values)
        ax_c.text(i, np.nanpercentile(values, 75) * 1.08,
                  f"median\n{median:.1f} yr",
                  ha="center", va="bottom", fontsize=20,
                  color=color, fontweight="bold")

    ax_c.set_xticklabels(labels, fontsize=24)
    ax_c.set_ylabel("Climate Residence Time (yr)", fontsize=24)
    ax_c.set_title("Climate Residence Time\nin Protected Areas",
                   fontsize=26, fontweight="bold", pad=5)
    ax_c.text(0.02, 0.98, "c", transform=ax_c.transAxes,
              fontsize=27, fontweight="bold", va="top")
    ax_c.grid(axis="y", alpha=0.25, linestyle="--", linewidth=0.7)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)

    # Panel d
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d2 = ax_d.twinx()

    x = np.array([0, 1])
    width = 0.35
    area_values = [area_245, area_585]
    patch_values = [patch_245, patch_585]

    area_bars = ax_d.bar(x - width / 2, area_values, width,
                         color=colors, alpha=0.85, label="Area")
    patch_bars = ax_d2.bar(x + width / 2, patch_values, width,
                           color=colors, alpha=0.45, hatch="///", label="Patch count")

    for bar, value in zip(area_bars, area_values):
        ax_d.text(bar.get_x() + bar.get_width() / 2, value + 0.02,
                  f"{value:.2f}", ha="center", va="bottom",
                  fontsize=21, fontweight="bold")
    for bar, value in zip(patch_bars, patch_values):
        ax_d2.text(bar.get_x() + bar.get_width() / 2, value + 10,
                   format_int(value), ha="center", va="bottom",
                   fontsize=21, fontweight="bold")

    ax_d.text(0.27, 0.90, f"{(area_585 - area_245) / area_245 * 100:+.1f}%",
              transform=ax_d.transAxes, fontsize=22,
              color=COLORS["ssp585"], fontweight="bold", ha="center")
    ax_d.text(0.73, 0.90, f"{(patch_585 - patch_245) / patch_245 * 100:+.1f}%",
              transform=ax_d.transAxes, fontsize=22,
              color=COLORS["ssp585"], fontweight="bold", ha="center")

    ax_d.set_xticks(x)
    ax_d.set_xticklabels(labels, fontsize=24)
    ax_d.set_ylabel("Core Refugia Area (million km²)", fontsize=22)
    ax_d2.set_ylabel("Number of Patches (≥50 km²)", fontsize=22, color="gray")
    ax_d.set_title("Core Refugia Area and Patch Count",
                   fontsize=26, fontweight="bold", pad=5)
    ax_d.text(0.02, 0.98, "d", transform=ax_d.transAxes,
              fontsize=27, fontweight="bold", va="top")
    ax_d.set_ylim(0, max(area_values) * 1.40)
    ax_d2.set_ylim(0, max(patch_values) * 1.40)
    ax_d.grid(axis="y", alpha=0.25, linestyle="--", linewidth=0.7)
    ax_d.spines["top"].set_visible(False)
    ax_d2.spines["top"].set_visible(False)

    fig.suptitle("Microclimate Refugia in the Amazon",
                 fontsize=27, fontweight="bold", y=0.985)

    CONFIG["output_dir"].mkdir(parents=True, exist_ok=True)
    out_pdf = CONFIG["output_dir"] / "Fig4_refugia_analysis.pdf"
    out_png = CONFIG["output_dir"] / "Fig4_refugia_analysis.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_pdf}")


def plot_gap_figure(gap):
    labels = [CONFIG["scenario_labels"]["ssp245"], CONFIG["scenario_labels"]["ssp585"]]
    covered = [gap["ssp245"]["covered_area"], gap["ssp585"]["covered_area"]]
    uncovered = [gap["ssp245"]["uncovered_area"], gap["ssp585"]["uncovered_area"]]
    covered_pct = [gap["ssp245"]["covered_pct"], gap["ssp585"]["covered_pct"]]

    fig, ax = plt.subplots(figsize=(8, 7))
    x = np.arange(2)
    width = 0.5

    ax.bar(x, [v / 1e6 for v in covered], width,
           color=COLORS["covered"], alpha=0.85, label="Within Protected Areas")
    ax.bar(x, [v / 1e6 for v in uncovered], width,
           bottom=[v / 1e6 for v in covered],
           color=COLORS["uncovered"], alpha=0.85, label="Outside Protected Areas")

    for i, (value, pct) in enumerate(zip(covered, covered_pct)):
        ax.text(i, (value / 1e6) / 2, f"{pct:.1f}%\ncovered",
                ha="center", va="center", fontsize=22,
                fontweight="bold", color="white")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=24)
    ax.set_ylabel("Refugia Area (million km²)", fontsize=24)
    ax.set_ylim(0, max([(c + u) / 1e6 for c, u in zip(covered, uncovered)]) * 1.3)
    ax.set_title("Gap Analysis — Refugia vs Protected Areas",
                 fontsize=26, fontweight="bold", pad=5)
    ax.legend(fontsize=21, framealpha=0.88, edgecolor="#cccccc", frameon=True, loc="upper right")
    ax.grid(axis="y", alpha=0.25, linestyle="--", linewidth=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out_pdf = CONFIG["output_dir"] / "Fig4_supplementary_gap.pdf"
    out_png = CONFIG["output_dir"] / "Fig4_supplementary_gap.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_pdf}")


# ---------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------

def main():
    CONFIG["output_dir"].mkdir(parents=True, exist_ok=True)

    pa_wgs84, pa_projected = read_protected_areas()
    refugia = read_refugia()

    residence = {}
    for scenario in CONFIG["scenarios"]:
        residence[scenario] = compute_residence_time(pa_wgs84, pa_projected, scenario)

    gap = {}
    for scenario in CONFIG["scenarios"]:
        gap[scenario] = compute_gap(refugia[scenario], pa_projected)
        gap[scenario]["gdf"].to_file(CONFIG["output_dir"] / f"refugia_gap_{scenario}.shp")

    residence_table = pd.concat(residence.values(), ignore_index=True)
    residence_table[["scenario", "area_km2", "diameter_km", "median_velocity", "residence_time"]].to_csv(
        CONFIG["output_dir"] / "residence_time_results.csv", index=False
    )

    plot_main_figure(pa_wgs84, refugia, residence, gap)
    plot_gap_figure(gap)


if __name__ == "__main__":
    main()
