# LS2D
(LS)^2D: LES and SCM - Large Scale Dynamics

Python package to download ERA5 (using CDS or MARS), and calculate the initial conditions and large scale forcings, for large eddy simulation (LES) or single column model (SCM) experiments.

### Required Python modules
- numpy
- matplotlib
- netCDF4
- cdsapi (*)

(*): required if you want to use the Copernicus Data Store (CDS) Python API. See https://cds.climate.copernicus.eu/api-how-to for information on how to setup the API.

### Setup

In the future, LS2D will be uploaded to PyPI (allowing for a simple `pip install LS2D`). For now, download the LS2D source code, e.g. with `git`:

    git clone https://github.com/julietbravo/LS2D.git
    
In each script where you want to use LS2D, add the LS2D root directory to the Python path:

    import sys
    sys.path.append('/path/to/LS2D')
    
### Usage

An example script is provided in `examples/example_1.py`, which downloads the ERA5 data, calculates the initial conditions and large scale forcings, and creates an example plot.

### The `settings` dictionary

All settings for LS2D are wrapped in a dictionary:

- `central_lat`: central latitude of LES/SCM domain
- `central_lon`: central longitude of LES/SCM domain
- `area_size`: spatial size of ERA5 download (central lat/lon +/- `area_size` degrees)
- `era5_path`: storage location of ERA5 downloads/data
- `era5_expver`: ERA5 experiment version number (`1`=normal ERA5, `5`=near realtime). With CDS, only `1` works.
- `case_name`: experiment name, only used to create subdirectory in `era5_path`.
- `start_date`: Python `datetime` object with start date/time
- `end_date`: Python `datetime` object with end date/time
- `write_log`: Write ERA5 download to screen (`False`) or log file (`True`)
- `data_source`: Download method (`CDS` or `MARS`)
- `ntasks`: option for parallel downloads (only for `CDS`).
