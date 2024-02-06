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

import os

#ls2d_root = '/Users/bart/meteo/models/LS2D/'
ls2d_root = '/home/bart/meteo/models/LS2D/'

old = 'Copyright (c) 2017-2023 Wageningen University & Research'
new = 'Copyright (c) 2017-2023 Wageningen University & Research'

for root, dirs, files in os.walk(ls2d_root):
    for file in files:
        bits = file.split('.')
        if len(bits) > 1 and bits[-1] == 'py':

            file = os.path.join(root, file)

            with open(file, 'r') as f:
                filedata = f.read()
            
            filedata = filedata.replace(old, new)
            
            with open(file, 'w') as f:
                f.write(filedata)
