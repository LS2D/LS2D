#! /bin/sh
microhh_path=/home/bart/meteo/models/microhh
rrtmgp_path=$microhh_path/rte-rrtmgp-cpp/rte-rrtmgp/

# RTE-RRTMGP:
ln -sf $rrtmgp_path/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc cloud_coefficients_sw.nc
ln -sf $rrtmgp_path/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc cloud_coefficients_lw.nc
ln -sf $rrtmgp_path/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc coefficients_sw.nc
ln -sf $rrtmgp_path/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc coefficients_lw.nc

# Land-surface model:
ln -sf $microhh_path/misc/van_genuchten_parameters.nc .
