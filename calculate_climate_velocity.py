"""
Climate velocity calculation workflow

This script computes:
1. Horizontal spatial temperature gradients and gradient directions from a baseline raster.
2. Temporal warming trends from a stack of decadal temperature rasters.
3. Climate-change velocity as the ratio between temporal trend and spatial gradient.

The script uses repository-relative paths and does not contain local machine-specific paths.
"""

from pathlib import Path
import math
import numpy as np
import rasterio
from scipy.stats import linregress


# =============================================================================
# Configuration
# =============================================================================

PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "outputs"

VELOCITY_INPUT_DIR = DATA_DIR / "velocity_inputs"
VELOCITY_OUTPUT_DIR = OUTPUT_DIR / "velocity_outputs"

VELOCITY_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)


# =============================================================================
# 1. Spatial gradient calculation
# =============================================================================

def calculate_average_max_gradient_and_direction(
    window,
    resolution_x,
    resolution_y,
    lat_center,
):
    """
    Compute horizontal spatial gradient magnitude and direction using a 3 x 3 window.

    Parameters
    ----------
    window : numpy.ndarray
        A 3 x 3 array containing the focal pixel and its neighbouring cells.
    resolution_x : float
        Pixel size in degrees along longitude.
    resolution_y : float
        Pixel size in degrees along latitude.
    lat_center : float
        Latitude of the focal pixel in degrees.

    Returns
    -------
    grad_mag : float
        Spatial gradient magnitude in units of input raster value per km.
    grad_dir : float
        Gradient direction in degrees.
    """
    center = window[1, 1]
    if np.isnan(center):
        return np.nan, np.nan

    lat_top = lat_center + resolution_y
    lat_mid = lat_center
    lat_bottom = lat_center - resolution_y

    width_top = resolution_x * 111.325 * math.cos(math.radians(lat_top))
    width_mid = resolution_x * 111.325 * math.cos(math.radians(lat_mid))
    width_bottom = resolution_x * 111.325 * math.cos(math.radians(lat_bottom))

    height = resolution_y * 111.352

    dx_top = (window[0, 2] - window[0, 0]) / (2 * width_top)
    dx_mid = (window[1, 2] - window[1, 0]) / width_mid
    dx_bottom = (window[2, 2] - window[2, 0]) / (2 * width_bottom)
    dx = (dx_top + dx_mid + dx_bottom) / 4.0

    dy_left = (window[2, 0] - window[0, 0]) / (2 * height)
    dy_mid = (window[2, 1] - window[0, 1]) / height
    dy_right = (window[2, 2] - window[0, 2]) / (2 * height)
    dy = (dy_left + dy_mid + dy_right) / 4.0

    grad_mag = np.sqrt(dx ** 2 + dy ** 2)
    grad_dir = np.degrees(np.arctan2(dy, dx))

    return grad_mag, grad_dir


def compute_spatial_gradient(
    input_file,
    output_gradient_file,
    output_direction_file,
    noise_range=(-0.05, 0.05),
):
    """
    Compute spatial gradient magnitude and direction from a baseline temperature raster.

    Parameters
    ----------
    input_file : str or pathlib.Path
        Input baseline temperature raster.
    output_gradient_file : str or pathlib.Path
        Output raster for spatial gradient magnitude.
    output_direction_file : str or pathlib.Path
        Output raster for spatial gradient direction.
    noise_range : tuple
        Small uniform noise range added to reduce artefacts in perfectly flat areas.
    """
    input_file = Path(input_file)
    output_gradient_file = Path(output_gradient_file)
    output_direction_file = Path(output_direction_file)

    output_gradient_file.parent.mkdir(parents=True, exist_ok=True)
    output_direction_file.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open(input_file) as src:
        data = src.read(1).astype(float)
        nodata = src.nodata
        if nodata is not None:
            data[data == nodata] = np.nan

        transform = src.transform
        resolution_x = abs(transform[0])
        resolution_y = abs(transform[4])
        lat_reference = src.bounds.top
        meta = src.meta.copy()

    data_with_noise = data + np.random.uniform(
        noise_range[0],
        noise_range[1],
        data.shape,
    )

    rows, cols = data.shape
    grad_mag = np.full_like(data, np.nan, dtype=float)
    grad_dir = np.full_like(data, np.nan, dtype=float)

    for row in range(1, rows - 1):
        lat_center = lat_reference - row * resolution_y

        for col in range(1, cols - 1):
            window = data_with_noise[row - 1: row + 2, col - 1: col + 2]

            if np.any(np.isnan(window)):
                continue

            g_h, d_h = calculate_average_max_gradient_and_direction(
                window=window,
                resolution_x=resolution_x,
                resolution_y=resolution_y,
                lat_center=lat_center,
            )

            grad_mag[row, col] = g_h
            grad_dir[row, col] = d_h

    meta.update(dtype="float32", nodata=np.nan, count=1)

    with rasterio.open(output_gradient_file, "w", **meta) as dst:
        dst.write(grad_mag.astype("float32"), 1)

    with rasterio.open(output_direction_file, "w", **meta) as dst:
        dst.write(grad_dir.astype("float32"), 1)

    print(f"Spatial gradient written to: {output_gradient_file}")
    print(f"Spatial direction written to: {output_direction_file}")


# =============================================================================
# 2. Temporal trend calculation
# =============================================================================

def compute_temporal_trend(
    input_folder,
    time_points,
    output_folder,
    file_suffix=".tif",
    nodata_value=-9999.0,
):
    """
    Compute per-pixel temporal linear trends from a stack of decadal rasters.

    Parameters
    ----------
    input_folder : str or pathlib.Path
        Folder containing decadal temperature rasters.
    time_points : list or numpy.ndarray
        Time coordinate for each raster, usually the central year of each decade.
    output_folder : str or pathlib.Path
        Output folder for slope, intercept and p-value rasters.
    file_suffix : str
        File suffix used to identify raster files.
    nodata_value : float
        No-data value written to output rasters.
    """
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    raster_files = sorted(input_folder.glob(f"*{file_suffix}"))

    if len(raster_files) != len(time_points):
        raise ValueError(
            f"Number of rasters ({len(raster_files)}) does not match "
            f"number of time points ({len(time_points)})."
        )

    raster_data = []
    meta = None

    for raster_path in raster_files:
        with rasterio.open(raster_path) as src:
            arr = src.read(1).astype(float)

            if src.nodata is not None:
                arr[arr == src.nodata] = np.nan

            raster_data.append(arr)

            if meta is None:
                meta = src.meta.copy()

    stack = np.stack(raster_data, axis=0)
    _, rows, cols = stack.shape

    slopes = np.full((rows, cols), np.nan, dtype=float)
    intercepts = np.full((rows, cols), np.nan, dtype=float)
    p_values = np.full((rows, cols), np.nan, dtype=float)

    time_points = np.asarray(time_points, dtype=float)
    skipped_pixels = 0

    for row in range(rows):
        for col in range(cols):
            y = stack[:, row, col]
            valid = np.isfinite(y)

            if valid.sum() < 2:
                skipped_pixels += 1
                continue

            slope, intercept, _, p_value, _ = linregress(
                time_points[valid],
                y[valid],
            )

            slopes[row, col] = slope
            intercepts[row, col] = intercept
            p_values[row, col] = p_value

    meta.update(dtype="float32", nodata=nodata_value, count=1)

    with rasterio.open(output_folder / "slopes.tif", "w", **meta) as dst:
        dst.write(slopes.astype("float32"), 1)

    with rasterio.open(output_folder / "intercepts.tif", "w", **meta) as dst:
        dst.write(intercepts.astype("float32"), 1)

    with rasterio.open(output_folder / "p_values.tif", "w", **meta) as dst:
        dst.write(p_values.astype("float32"), 1)

    print(f"Skipped {skipped_pixels} pixels with fewer than two valid observations.")
    print(f"Temporal trend rasters written to: {output_folder}")


# =============================================================================
# 3. Climate velocity calculation
# =============================================================================

def compute_climate_velocity(
    time_gradient_file,
    spatial_gradient_file,
    output_file,
    min_gradient_threshold=1e-6,
):
    """
    Compute climate-change velocity as temporal trend divided by spatial gradient.

    Parameters
    ----------
    time_gradient_file : str or pathlib.Path
        Raster containing temporal warming trend.
    spatial_gradient_file : str or pathlib.Path
        Raster containing spatial temperature gradient.
    output_file : str or pathlib.Path
        Output climate velocity raster.
    min_gradient_threshold : float
        Minimum spatial gradient threshold used to avoid division by very small values.
    """
    time_gradient_file = Path(time_gradient_file)
    spatial_gradient_file = Path(spatial_gradient_file)
    output_file = Path(output_file)

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open(time_gradient_file) as src_time:
        time_gradient = src_time.read(1).astype(float)
        time_meta = src_time.meta.copy()

        if src_time.nodata is not None:
            time_gradient[time_gradient == src_time.nodata] = np.nan

    with rasterio.open(spatial_gradient_file) as src_spatial:
        spatial_gradient = src_spatial.read(1).astype(float)

        if src_spatial.nodata is not None:
            spatial_gradient[spatial_gradient == src_spatial.nodata] = np.nan

    safe_spatial_gradient = np.where(
        spatial_gradient > min_gradient_threshold,
        spatial_gradient,
        np.nan,
    )

    climate_velocity = time_gradient / safe_spatial_gradient

    time_meta.update(dtype="float32", nodata=np.nan, count=1)

    with rasterio.open(output_file, "w", **time_meta) as dst:
        dst.write(climate_velocity.astype("float32"), 1)

    print(f"Climate velocity raster written to: {output_file}")


# =============================================================================
# 4. Example workflow
# =============================================================================

def main():
    """
    Example workflow using repository-relative paths.

    The expected input structure is:

    data/
    └── velocity_inputs/
        ├── mean_temperature_baseline.tif
        └── decadal_mean_temperature/
            ├── decade_2020s.tif
            ├── decade_2030s.tif
            └── ...
    """
    spatial_input = VELOCITY_INPUT_DIR / "mean_temperature_baseline.tif"

    spatial_gradient_output = VELOCITY_OUTPUT_DIR / "spatial_gradient.tif"
    spatial_direction_output = VELOCITY_OUTPUT_DIR / "spatial_direction.tif"

    temporal_input_folder = VELOCITY_INPUT_DIR / "decadal_mean_temperature"
    temporal_output_folder = VELOCITY_OUTPUT_DIR / "temporal_trend"

    velocity_output = VELOCITY_OUTPUT_DIR / "climate_velocity.tif"

    # Example central years for decadal future rasters.
    # Modify this list if your input rasters use different time periods.
    time_points = [2028, 2038, 2048, 2058, 2068, 2078, 2088, 2096]

    compute_spatial_gradient(
        input_file=spatial_input,
        output_gradient_file=spatial_gradient_output,
        output_direction_file=spatial_direction_output,
        noise_range=(-0.05, 0.05),
    )

    compute_temporal_trend(
        input_folder=temporal_input_folder,
        time_points=time_points,
        output_folder=temporal_output_folder,
        file_suffix=".tif",
        nodata_value=-9999.0,
    )

    compute_climate_velocity(
        time_gradient_file=temporal_output_folder / "slopes.tif",
        spatial_gradient_file=spatial_gradient_output,
        output_file=velocity_output,
        min_gradient_threshold=1e-6,
    )


if __name__ == "__main__":
    main()
