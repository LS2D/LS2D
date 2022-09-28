from datetime import datetime, timedelta
import sys

"""
Lookup table, linking variable names to the full
names (needed by CDS/ADS) and the Grib IDs (needed by Mars).
"""
cams_lut = {
    # Meteorology:
    'ciwc': {'name': 'specific_cloud_ice_water_content', 'id': 247},
    'clwc': {'name': 'specific_cloud_liquid_water_content', 'id': 246},
    'q': {'name': 'specific_humidity', 'id': 133},
    'crwc': {'name': 'specific_rain_water_content', 'id': 75},
    'cswc': {'name': 'specific_snow_water_content', 'id': 76},
    'lnsp': {'name': 'surface_pressure', 'id': 152},
    't': {'name': 'temperature', 'id': 130},
    # CAMS chemistry:
    'ch3coch3': {'name': 'acetone', 'id': 217052},
    'aco2': {'name': 'acetone_product', 'id': 217053},
    'ald2': {'name': 'aldehydes', 'id': 217012},
    'nh2': {'name': 'amine', 'id': 217040},
    'nh3': {'name': 'ammonia', 'id': 217019},
    'nh4': {'name': 'ammonium', 'id': 217021},
    'co': {'name': 'carbon_monoxide', 'id': 210123},
    'dms': {'name': 'dimethyl_sulfide', 'id': 217018},
    'n2o5': {'name': 'dinitrogen_pentoxide', 'id': 217033},
    'c2h6': {'name': 'ethane', 'id': 217045},
    'c2h5oh': {'name': 'ethanol', 'id': 217046},
    'c2h4': {'name': 'ethene', 'id': 217010},
    'hcho': {'name': 'formaldehyde', 'id': 210124},
    'hcooh': {'name': 'formic_acid', 'id': 217043},
    'h2o2': {'name': 'hydrogen_peroxide', 'id': 217003},
    'ho2': {'name': 'hydroperoxy_radical', 'id': 217028},
    'oh': {'name': 'hydroxyl_radical', 'id': 217030},
    'c5h8': {'name': 'isoprene', 'id': 217016},
    'pb': {'name': 'lead', 'id': 217026},
    'ispd': {'name': 'methacrolein_mvk', 'id': 217050},
    'mcooh': {'name': 'methacrylic_acid', 'id': 217044},
    'ch4_c': {'name': 'methane_chemistry', 'id': 217004},
    'msa': {'name': 'methane_sulfonic_acid', 'id': 217022},
    'ch3oh': {'name': 'methanol', 'id': 217042},
    'ch3cocho': {'name': 'methyl_glyoxal', 'id': 217023},
    'ch3ooh': {'name': 'methyl_peroxide', 'id': 217007},
    'ch3o2': {'name': 'methylperoxy_radical', 'id': 217029},
    'no3_a': {'name': 'nitrate', 'id': 217051},
    'no3': {'name': 'nitrate_radical', 'id': 217032},
    'hno3': {'name': 'nitric_acid', 'id': 217006},
    'no2': {'name': 'nitrogen_dioxide', 'id': 210121},
    'no': {'name': 'nitrogen_monoxide', 'id': 217027},
    'ole': {'name': 'olefins', 'id': 217011},
    'ror': {'name': 'organic_ethers', 'id': 217036},
    'onit': {'name': 'organic_nitrates', 'id': 217015},
    'o3': {'name': 'ozone', 'id': 203},
    'par': {'name': 'paraffins', 'id': 217009},
    'ho2no2': {'name': 'pernitric_acid', 'id': 217034},
    'rooh': {'name': 'peroxides', 'id': 217014},
    'c2o3': {'name': 'peroxy_acetyl_radical', 'id': 217035},
    'pan': {'name': 'peroxyacetyl_nitrate', 'id': 217013},
    'c3h8': {'name': 'propane', 'id': 217047},
    'c3h6': {'name': 'propene', 'id': 217048},
    'ra': {'name': 'radon', 'id': 210181},
    'so2': {'name': 'sulphur_dioxide', 'id': 210122},
    'c10h16': {'name': 'terpenes', 'id': 217049},
    # CAMS aerosols:
    'aermr04': {'name': 'dust_aerosol_0.03-0.55um_mixing_ratio', 'id': 210004},
    'aermr05': {'name': 'dust_aerosol_0.55-0.9um_mixing_ratio', 'id': 210005},
    'aermr06': {'name': 'dust_aerosol_0.9-20um_mixing_ratio', 'id': 210006},
    'aermr09': {'name': 'hydrophilic_black_carbon_aerosol_mixing_ratio', 'id': 210009},
    'aermr07': {'name': 'hydrophilic_organic_matter_aerosol_mixing_ratio', 'id': 210007},
    'aermr10': {'name': 'hydrophobic_black_carbon_aerosol_mixing_ratio', 'id': 210010},
    'aermr08': {'name': 'hydrophobic_organic_matter_aerosol_mixing_ratio', 'id': 210008},
    'aermr01': {'name': 'sea_salt_aerosol_0.03-0.5um_mixing_ratio', 'id': 210001},
    'aermr02': {'name': 'sea_salt_aerosol_0.5-5um_mixing_ratio', 'id': 210002},
    'aermr03': {'name': 'sea_salt_aerosol_5-20um_mixing_ratio', 'id': 210003},
    'aermr11': {'name': 'sulphate_aerosol_mixing_ratio', 'id': 210011}
}

def download_cams(settings, variables):

    """
    1. Check if variables are in lookup table.
    """
    success = True
    for variable in variables:
        if variable not in cams_lut.keys():
            print('ERROR: invalid CAMS variable \"{}\"'.format(variable))
            success = False
    if not success:
        sys.exit('CAMS download failed!')

    """
    2. Create list of dates to download.
       Each day has data for 00, 03, 06, .., 21 UTC.
    """
    start = settings['start_date']
    end   = settings['end_date']

    # Strip start hour:
    date = datetime(start.year, start.month, start.day)

    dates = []
    while date <= end:
        dates.append(date)
        date += timedelta(hours=24)


    





if __name__ == '__main__':
    from datetime import datetime

    start = datetime(year=2016, month=8, day=1, hour=6)
    end   = datetime(year=2016, month=8, day=2, hour=22)

    settings = {
            'central_lat': 51,
            'central_lon': 4,
            'start_date': start,
            'end_date': end}

    species = ['no', 'no2', 'o3', 'oh']

    download_cams(settings, species)
    #cams = ls2d.read_cams(settings, species)
