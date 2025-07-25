#! /bin/sh
microhh_path=/home/bart/meteo/models/microhh
rrtmgp_path=$microhh_path/rte-rrtmgp-cpp/rrtmgp-data/

# RTE-RRTMGP:
ln -sf $rrtmgp_path/rrtmgp-clouds-sw.nc cloud_coefficients_sw.nc
ln -sf $rrtmgp_path/rrtmgp-clouds-lw.nc cloud_coefficients_lw.nc
ln -sf $rrtmgp_path/rrtmgp-gas-sw-g112.nc coefficients_sw.nc
ln -sf $rrtmgp_path/rrtmgp-gas-lw-g128.nc coefficients_lw.nc

# Land-surface model:
ln -sf $microhh_path/misc/van_genuchten_parameters.nc .
