from datetime import datetime

# Add `src` subdirectory to Python path
import sys; sys.path.append('src/')

# Import the LS2D specific scripts
from download_ERA5 import download_ERA5_period
from read_ERA5 import Read_ERA

# ------------------------
# Settings
# ------------------------
central_lat  = 51.971
central_lon  = 4.927
area_size    = 2   # ERA5 area size (lat+/-size, lon+/-size degrees)

case_name    = 'cabauw'
#ERA5_path    = '/Users/bart/meteo/data/ERA5/LS2D/'
ERA5_path    = '/nobackup/users/stratum/ERA5/LS2D/'
model        = 'MicroHH'
LES_path     = '/Users/bart/meteo/models/MicroHH/testbed_cabauw/'

# Start and end date/time of experiment. For now, limited to full hours
start_date   = datetime(year=2016, month=5, day=1, hour=5)
end_date     = datetime(year=2016, month=5, day=2, hour=19, minute=23)
# ------------------------

# 1. Download the ERA5 data
download_ERA5_period(start_date, end_date, central_lat, central_lon, area_size, ERA5_path, case_name)

# 2. Read the ERA5 data
e5 = Read_ERA(start_date, end_date, central_lat, central_lon, ERA5_path, case_name)

# 3. Calculate the large-scale forcings from the ERA5 data
e5.calculate_forcings(n_av=1)

# 4. Write model specific output
if model == 'MicroHH':
    pass
elif model == 'DALES':
    pass
