#
# This file is part of LS2D.
#
# Copyright (c) 2017-2026 Wageningen University & Research
# Author: Bart van Stratum (WUR)
#
# LS2D is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LS2D is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LS2D.  If not, see <http://www.gnu.org/licenses/>.
#

# Python modules
from datetime import datetime

import ls2d

settings = {
    'central_lat': 51.97,
    'central_lon': 4.93,
    'area_size': 0.25,
    'case_name': 'cabauw',
    'era5_path': '/home/scratch1/meteo_data/LS2D_ERA5_ARCO/',
    'start_date': datetime(year=2016, month=8, day=15, hour=6),
    'end_date': datetime(year=2016, month=8, day=15, hour=18),
}

# Download required ERA5 files from Google ARCO.
ls2d.download_era5_arco(settings)
