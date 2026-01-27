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

import datetime
import sys

_opts = {
   'blue'   : '\033[94m',
   'green'  : '\033[92m',
   'purple' : '\033[95m',
   'red'    : '\033[91m',
   'bf'     : '\033[1m',
   'ul'     : '\033[4m',
   'end'    : '\033[0m'
}

def header(message, time=True):
    """
    Format of print statements indicating new main routine
    """
    if time:
        now = datetime.datetime.now()
        print('{}{}{}{} {}[{}]{}'.format(_opts['blue'], _opts['bf'], message, _opts['end'], _opts['green'], now.strftime('%d-%m: %H:%M'), _opts['end']))
    else:
        print('{}{}{}{}'.format(_opts['blue'], _opts['bf'], message, _opts['end']))

def message(message):
    """
    Format of print statements
    """
    print(' - {}'.format(message))

def warning(message):
    """
    Format of print warnings
    """
    print('{}{}WARNING:{} {}'.format(_opts['purple'], _opts['bf'], _opts['end'], message))

def error(message, exit=True):
    """
    Format of print errors
    """
    print('{}{}ERROR:{} {}'.format(_opts['red'], _opts['bf'], _opts['end'], message))
    if exit:
        sys.exit()

