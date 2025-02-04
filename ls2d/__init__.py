#
# This file is part of LS2D.
#
# Copyright (c) 2017-2024 Wageningen University & Research
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

# Ban Python 2.x:
import sys
if sys.version_info.major < 3:
    from ls2d.src.messages import error
    error('(LS)2D requires Python 3.x')

# Make packages directly available as e.g.:
# ls2d.download_era5() instead of ls2d.ecmwf.download_era5()
from ls2d.ecmwf import download_era5
from ls2d.ecmwf import download_cams

from ls2d.ecmwf import Read_era5
from ls2d.ecmwf import Read_cams

from ls2d.src import grid
from ls2d.src import logger
