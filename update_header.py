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

import glob

replace_in = '# Copyright (c) 2017-2024 Wageningen University & Research'
replace_out = '# Copyright (c) 2017-2024 Wageningen University & Research'

files = glob.iglob('**/*.py', recursive=True)
for file in files:
    with open(file, 'r') as f:
        lines = f.readlines()
    with open(file, 'w') as f:
        for l in lines:
            f.write(l.replace(replace_in, replace_out))
