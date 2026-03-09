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

from ls2d import logger

# Demo of all log levels
logger.setLevel('DEBUG')

logger.debug('This is a debug() message')
logger.info('This is an info() message')
logger.warning('This is a warning() message')
logger.error('This is an error() message that does not raise exception')

# Suppress DEBUG messages (only INFO and above)
logger.setLevel('INFO')
logger.debug('This debug() message will not be shown..')
logger.info('This info() message will be shown')

# critical() raises RuntimeError in addition to logging
logger.critical('This is a critical() message that raises RuntimeError')