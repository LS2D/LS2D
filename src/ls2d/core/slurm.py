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

# Pythom modules
import subprocess
import datetime
import glob
import sys


def execute(task):
    """
    Execute command line `task`
    """
    subprocess.call(task, shell=True, executable='/bin/bash')


def update_namelist(namelist_file, to_update):
    """
    Update `namelist_file` with dictionary `to_update`
    """
    with open(namelist_file, 'r') as f:
        lines = f.readlines()
    with open(namelist_file, 'w') as f:
        for l in lines:
            if len(l) > 0 and l[0] != '#' and '=' in l:
                name  = l.split('=')[0].strip()
                if name in to_update.keys():
                    f.write('{}={}\n'.format(name, to_update[name]))
                else:
                    f.write(l)
            else:
                f.write(l)


def read_namelist_value(namelist_file, variable):
    """
    Read single value from .ini file (first occurance only!)
    """
    with open(namelist_file, 'r') as f:
        for l in f.readlines():
            if len(l) > 0 and l[0] != '#' and '=' in l:
                name,value = l.split('=')
                if name == variable:
                    return value.strip()


def get_latest_restart_time(path='.'):
    """
    Get latest restart time from the `time.***` files
    """
    files = glob.glob('{}/time.*'.format(path))
    if len(files) == 0:
        print('Cant find any time.xxxx files, returning 0')
        return 0
    files.sort()
    return int(files[-1].split('.')[-1])


def seconds_to_time(seconds):
    """
    Convert seconds to d-hh:mm:ss format
    """
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)

    return '{0}-{1:02d}:{2:02d}:{3:02d}'.format(d,h,m,s)


def submit_case(
        case, total_time, max_time, wc_time, n_tasks,
        partition, work_dir, job_name, script_name,
        auto_submit, bypass_slurm=False):

    # Get time of latest restart file
    restart_time = get_latest_restart_time(work_dir)

    # Check if the experiment is finished, if not, re-submit
    if restart_time == total_time:
        print('Experiment finished!')
    else:
        # Switch between warm and cold start
        is_cold_start = restart_time==0

        if is_cold_start:
            print('Submitting cold start experiment')
        else:
            print('Re-submitting experiment')

        # Fix end time in case the wallclocklimit was hit,
        # to return to the normal `max_time` cycles.
        if restart_time % max_time != 0:
            end_time = (restart_time + max_time) // max_time * max_time
        else:
            end_time = restart_time + max_time

        # Update namelist
        to_update = {
                'wallclocklimit': wc_time/3600.,
                'starttime': restart_time,
                'endtime': end_time,
                'savetime': max_time}

        update_namelist('{}/{}.ini'.format(work_dir, case), to_update)

        # Wall-clock limit in `d-hh-mm-ss` format
        wc_time_str = seconds_to_time(int(wc_time))

        auto_submit_flag = '--auto_submit' if auto_submit else ''

        # Create Slurm run script
        with open('{}/{}'.format(work_dir, script_name), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH -p {}\n'.format(partition))
            f.write('#SBATCH -n {}\n'.format(n_tasks))
            f.write('#SBATCH -t {}\n'.format(wc_time_str))
            f.write('#SBATCH --job-name={}\n'.format(job_name))
            f.write('#SBATCH --output={}/mhh-%j.out\n'.format(work_dir))
            f.write('#SBATCH --error={}/mhh-%j.err\n'.format(work_dir))
            f.write('#SBATCH --constraint=haswell\n\n')

            f.write('cd {}\n\n'.format(work_dir))

            if bypass_slurm:
                if is_cold_start:
                    f.write('./microhh init {}\n'.format(case))
                f.write('./microhh run {}\n\n'.format(case))
            else:
                if is_cold_start:
                    if n_tasks == 1:
                        f.write('./microhh init {}\n'.format(case))
                    else:
                        f.write('srun ./microhh init {}\n'.format(case))

                if n_tasks == 1:
                    f.write('./microhh run {}\n\n'.format(case))
                else:
                    f.write('srun ./microhh run {}\n\n'.format(case))

            # Only re-submit if case exit was normal
            f.write('runjobval=$?\n')
            f.write('if [ "${runjobval}" -eq "0" ]; then\n')
            # Re-submit Python script
            f.write('  python3 {}/slurm.py -c {} -tt {} -mt {} -wc {} -p {} -n {} -w {} -j {} -s {} {}\n'.format(
                work_dir, case, int(total_time), int(max_time), int(wc_time), partition,
                int(n_tasks), work_dir, job_name, script_name, auto_submit_flag))
            f.write('fi\n')

        # Submit runscript
        if auto_submit:
            if bypass_slurm:
                execute('chmod +x {}/{}'.format(work_dir, script_name))
                execute('{}/{}'.format(work_dir, script_name))
            else:
                execute('sbatch {}/{}'.format(work_dir, script_name))


if __name__ == '__main__':
    import argparse

    #
    # Parse the input arguments
    #
    parser = argparse.ArgumentParser(
        description='Slurm launcher with restart capabilities')

    # Required arguments:
    parser.add_argument('-tt', '--total_time', required=True, type=int,
            help='Total integration time of experiment')
    parser.add_argument('-mt', '--max_time', required=True, type=int,
            help='Max integration time per submission')
    parser.add_argument('-c', '--case', required=True,
            help='MicroHH case name')
    parser.add_argument('-n', '--n_tasks', required=True, type=int,
            help='Total MPI tasks')
    parser.add_argument('-wc', '--wallclocklimit', required=True, type=int,
            help='Wall clock limit in seconds')

    # Optional arguments:
    parser.add_argument('--auto_submit', dest='auto_submit', action='store_true')

    parser.add_argument('-w', '--work_dir', default='.',
            help='Work directory')
    parser.add_argument('-j', '--job_name', default='mhh',
            help='Slurm job name')
    parser.add_argument('-s', '--script_name', default='run_restart.slurm',
            help='Slurm script name')
    parser.add_argument('-p', '--partition', default='normal',
            help='Slurm partition')

    args = parser.parse_args()

    #
    # Submit case
    #
    submit_case(
        args.case,
        args.total_time,
        args.max_time,
        args.wallclocklimit,
        args.n_tasks,
        args.partition,
        args.work_dir,
        args.job_name,
        args.script_name,
        args.auto_submit)
