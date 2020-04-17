#
# This file is part of LS2D.
#
# Copyright (c) 2017-2018 Bart van Stratum
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

from messages import *

# One day datetime offset
one_day = datetime.timedelta(days=1)

def get_required_analysis(start, end):

    # Analysis start at 00 UTC, so first analysis = start day
    first_analysis = datetime.datetime(start.year, start.month, start.day)

    # If end time is after 23 UTC, include next day for the analysis files
    if end.hour == 23 and end.minute > 0:
        last_analysis = datetime.datetime(end.year, end.month, end.day) + one_day
    else:
        last_analysis = datetime.datetime(end.year, end.month, end.day)

    # Create list of datetime objects:
    dates = [first_analysis + i*one_day for i in range((last_analysis-first_analysis).days + 1)]

    return dates


def get_required_forecast(start, end):

    # Forecast runs through midnight, so last analysis = last day
    last_forecast = datetime.datetime(end.year, end.month, end.day)

    # If start time is before 06 UTC, include previous day for the forecast files
    if start.hour > 6:
        first_forecast = datetime.datetime(start.year, start.month, start.day)
    else:
        first_forecast = datetime.datetime(start.year, start.month, start.day) - one_day

    # Create list of datetime objects:
    dates = [first_forecast + i*one_day for i in range((last_forecast-first_forecast).days + 1)]

    return dates


def lower_to_hour(time):
    time_out = datetime.datetime(time.year, time.month, time.day, time.hour)
    if time.minute != 0 or time.second != 0:
        warning('changed date/time from {} to {}'.format(time, time_out))
    return time_out
