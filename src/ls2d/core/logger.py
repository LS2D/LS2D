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
import logging

# Third party modules
from colorlog import ColoredFormatter

logger = logging.getLogger('ls2d')
logger.setLevel(logging.DEBUG)

if not logger.handlers:
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(
        ColoredFormatter(
            '[%(asctime)s] [%(name)s] %(log_color)s[%(levelname)s] %(message)s%(reset)s',
            datefmt='%Y/%m/%d %H:%M:%S',
            log_colors={
                'DEBUG': 'fg_244',
                'INFO': '',
                'WARNING': 'fg_208',
                'ERROR': 'red',
                'CRITICAL': 'red',
            },
        )
    )
    logger.addHandler(handler)
    logger.propagate = False
