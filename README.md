# (LS)<sup>2</sup>D: LES and SCM - Large Scale Dynamics

[![PyPI version](https://badge.fury.io/py/ls2d.svg)](https://pypi.org/project/ls2d/)
[![Static Badge](https://img.shields.io/badge/JAMES-10.1029%2F2023MS003750-blue?link=https%3A%2F%2Fdoi.org%2F10.1029%2F2023MS003750)](https://doi.org/10.1029/2023MS003750)

(LS)<sup>2</sup>D is a Python toolkit, developed to simplify all the steps required to downscale ERA5 with doubly-periodic large-eddy simulation (LES), or single-column models (SCMs). For the retrieval of ERA data, it relies on the Google Analysis-Ready, Cloud-Optimized (ARCO) ERA5 archive, the Copernicus Data Store (CDS), or the Meteorological Archival and Retrieval System (MARS) at ECMWF computer systems.

> [!IMPORTANT]
> Since v1.1.0, (LS)<sup>2</sup>D has two paths to download and process ERA5:
> - Google [ARCO ERA5](https://github.com/google-research/arco-era5) (recommended): `ls2d.download_era5_arco()`, `ls2d.read_era5_arco()` and `ls2d.create_column_input()`; see `examples/example_era5_arco.py`.
> - CDS or MARS: `ls2d.download_era5()` and `ls2d.Read_era5`; see `examples/example_era5.py`. CDS downloads are extremely slow, and further slowed down by the new CDS request size limits (~09/2026), but CDS is potentially more up-to-date than ARCO.
>
> Note that:
> - The ARCO path is not backward compatible: ERA5 files downloaded from CDS or MARS can't be read by `read_era5_arco()`. They can still be read with `Read_era5`.
> - Results from the ARCO and CDS paths are close but not identical, due to small changes in e.g. the interpolations.
> - ARCO stores global fields, and the smallest part that can be read is a band of ~91 latitudes x 1440 longitudes. Downloads are therefore much larger than the saved data: one day for a 1x1° area is ~7 MB on disk, but ~5-10 GB to download. See `src/ls2d/google/arco_tools.py` for details.

> [!IMPORTANT]
> In September 2024, Copernicus upgraded the Copernicus and Atmosphere Data Stores (CDS/ADS), introducing changes in the NetCDF file structure. The updated files now have different dimension names, time units, and even some reordered dimensions (!!). Unfortunately, these changes made the new NetCDF files incompatible with older versions or files extracted from MARS. Although (LS)<sup>2</sup>D briefly supported the new format, this led to several issues. Therefore, starting November 2024, (LS)<sup>2</sup>D will automatically patch new NetCDF files to ensure compatibility with the previous NetCDF file structure. On the user side, nothing should change:
> - Old NetCDF files from CDS or MARS remain compatible.
> - New NetCDF files from CDS/ADS are automatically patched after download.
> - Unpatched new NetCDF files are patched once prior to reading.
>   
> To upgrade to the new CDS and ADS, please follow the steps described [here](https://confluence.ecmwf.int/display/CKB/Please+read%3A+CDS+and+ADS+migrating+to+new+infrastructure%3A+Common+Data+Store+%28CDS%29+Engine). If (LS)<sup>2</sup>D was installed using `pip`, upgrade (LS)<sup>2</sup>D with `pip install --upgrade ls2d`.

### References

(LS)<sup>2</sup>D is described in:

B.J.H. van Stratum, C.C. van Heerwaarden, & J. Vilà-Guerau de Arellano (2023). *The benefits and challenges of downscaling a global reanalysis with doubly-periodic large-eddy simulations.* JAMES, https://doi.org/10.1029/2023MS003750

If you use (LS)<sup>2</sup>D, we kindly request citing this paper.

### Installation

If you want to use CDS to download the ERA5 data, then please start by following the steps explained at https://cds.climate.copernicus.eu/how-to-api .

#### PyPI

It is easiest to install (LS)<sup>2</sup>D from PyPI:

    pip install ls2d
    
By default, this excludes the `cdsapi` as a dependency. If you do want to install that as a dependency, use:
    
    pip install ls2d[cds]
   
#### Manual

For a manual installation, you can clone the package from Github, and install it in editable mode:

    git clone https://github.com/LS2D/LS2D.git
    cd LS2D
    pip install -e .

Or use `pip install -e .[cds]` to include the `cdsapi` dependency. Changes to the local source code are directly available without re-installing.
    
### Usage

Some examples are provided at https://github.com/LS2D/LS2D/tree/main/examples. The script `example_1.py` downloads the ERA5 data, calculates the initial conditions and large scale forcings, and creates an example plot.

The examples directory also contains example cases for MicroHH (https://github.com/microhh/microhh).

### Downloading ERA5 data

(LS)<sup>2</sup>D contains three methods to download ERA5: from the Google [ARCO ERA5](https://github.com/google-research/arco-era5) archive (recommended), through the [Copernicus Data Store](https://cds.climate.copernicus.eu) (open for everyone), or using MARS at ECMWF systems. 

The ERA5 model level data that (LS)<sup>2</sup>D requires is stored on tape archives, so downloads using CDS are extremely slow with long queueing times. For that reason, `ls2d.download_era5()` will stop the Python script once the download requests are submitted to CDS. On subsequent calls of `ls2d.download_era5()`, (LS)<sup>2</sup>D will check the status of the CDS request, and if the request is finished, download the ERA5 data. Due to CDS request size limits, the model level data is requested in blocks of 6 hours, which are merged after downloading.

### The `settings` dictionary

All settings for (LS)<sup>2</sup>D are wrapped in a dictionary:

- `central_lat`: central latitude of LES/SCM domain
- `central_lon`: central longitude of LES/SCM domain
- `area_size`: spatial size of ERA5 download (central lat/lon +/- `area_size` degrees)
- `era5_path`: storage location of ERA5 downloads/data
- `era5_expver`: ERA5 experiment version number (`1`=normal ERA5, `5`=near realtime). With CDS, only `1` works.
- `case_name`: experiment name, only used to create subdirectory in `era5_path`.
- `start_date`: Python `datetime` object with start date/time
- `end_date`: Python `datetime` object with end date/time
- `write_log`: Write ERA5 download to screen (`False`) or log file (`True`)
- `data_source`: Download method (`CDS` or `MARS`). `MARS` only works on e.g. the ECMWF supercomputer.

### Logging

By default, (LS)<sup>2</sup>D prints formatted log messages to the screen. To use your own logging setup instead, call `ls2d.disable_default_logging()`; messages are then passed on to your own handlers via the `ls2d` logger.

## Contributing guidelines

Contributions through pull requests are always appreciated. To keep the code style consistent, we stick to PEP8 with two modifications (see [`project.toml`](./pyproject.toml)):

- The maximum line length is set to 120 characters instead of 80.
- Single quotes are used for strings instead of double quotes.

To ensure consistent formatting, the code should be automatically formatted using [ruff](https://docs.astral.sh/ruff/). It can be installed via pip and can be run on the entire codebase from the project root (this directory) like

    ruff format .
