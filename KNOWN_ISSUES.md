# Known Issues

## [Issue 1: Longitude range](#issue-1-longitude-range)

Depending on the central longitude of the CAMS downloads, longitude values in the NetCDF files may be in the range `0 to 360` or `-180 to 180`. Unfortunately, the ADS system is inconsistent in this behavior (see https://github.com/LS2D/LS2D/issues/8#issuecomment-3149595007).

To ensure consistency, we patch any NetCDF files that use longitudes in the 0 to 360 range. This patching must be done immediately after download and before regridding. As a result, older NetCDF files that were downloaded but not patched are incompatible and must be re-downloaded.



