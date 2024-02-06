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

def grad2(a, b, delta):
    """
    2nd order accurate gradient at location of X:

        delta
       <----->
    |  a  |  b  |
         [X]
    """

    return (b - a) / delta

def grad2c(a, b, delta):
    """
    2nd order accurate gradient at location of X:

           delta
          <----->
    |  a  |     |  b  |
            [X]
    """

    return (b - a) / (2*delta)

def grad4(a, b, c, d, delta):
    """
    4th order accurate gradient at location of X:

              delta
             <----->
    |  a  |  b  |  c  |  d  |
               [X]
    """

    return (a - 27*b + 27*c - d) / (24*delta)

def grad4c(a, b, c, d, delta):
    """
    4th order accurate gradient at location of X:

                 delta
                <----->
    |  a  |  b  |     |  c  |  d  |
                  [X]
    """

    return (a - 8*b + 8*c - d) / (12*delta)
