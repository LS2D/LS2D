# LS2D
(LS)<sup>2</sup>D: LES and SCM - Large Scale Dynamics

[![PyPI version](https://badge.fury.io/py/ls2d.svg)](https://pypi.org/project/ls2d/)

(LS)<sup>2</sup>D is a Python toolkit, developed to simplify all the steps required to downscale ERA5 with doubly-periodic large-eddy simulation (LES), or single-column models (SCMs). For the retrieval of ERA data, it relies on the Copernicus Data Store (CDS), or the Meteorological Archival and Retrieval System (MARS) at ECMWF computer systems.

### References

(LS)<sup>2</sup>D is described in:

B.J.H. van Stratum, C.C. van Heerwaarden, & J. Vilà-Guerau de Arellano (2023). *The benefits and challenges of downscaling a global reanalysis with doubly-periodic large-eddy simulations.* JAMES, https://doi.org/10.1029/2023MS003750

If you use (LS)<sup>2</sup>D, we kindly request citing this paper.

### Installation

If you want to use CDS to download the ERA5 data, then please start by following the steps explained at https://cds.climate.copernicus.eu/api-how-to .

#### PyPI

It is easiest to install (LS)<sup>2</sup>D from PyPI:

    pip install ls2d
    
By default, this excludes the `cdsapi` as a dependency. If you do want to install that as a dependency, use:
    
    pip install ls2d[cds]
   
#### Manual

For a manual installation, you can clone the package from Github:

    git clone https://github.com/LS2D/LS2D.git

In each script where you want to use (LS)<sup>2</sup>D, add the (LS)<sup>2</sup>D root directory to the Python path:

    import sys
    sys.path.append('/path/to/LS2D')
    
You will have to manually install the dependencies with `pip install numpy scipy netCDF4 matplotlib cdsapi`.
    
### Usage

Some examples are provided at https://github.com/LS2D/LS2D/tree/main/examples. The script `example_1.py` downloads the ERA5 data, calculates the initial conditions and large scale forcings, and creates an example plot.

The examples directory also contains example cases for MicroHH (https://github.com/microhh/microhh).

### Downloading ERA5 data

(LS)<sup>2</sup>D contains two methods to download ERA5: through the [Copernicus Data Store](https://cds.climate.copernicus.eu) (open for everyone), or using MARS at ECMWF systems. 

The ERA5 model level data that (LS)<sup>2</sup>D requires is stored on tape archives, so downloads using CDS tend to be slow with long queueing times. For that reason, `ls2d.download_era5()` will stop the Python script once the download requests are submitted to CDS. On subsequent calls of `ls2d.download_era5()`, (LS)<sup>2</sup>D will check the status of the CDS request, and if the request is finished, download the ERA5 data. 

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
