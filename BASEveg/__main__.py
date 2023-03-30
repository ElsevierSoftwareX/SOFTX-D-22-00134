# Script to run the vegetation model in BM3
from functions.input_output import *
from vegetation_update import main_veg_module
from procedures import run_basement

import os
import csv
import sys
import re
import subprocess
from shutil import copyfile, copy
import platform
import argparse
import textwrap
import time

# Add to the python path the vegetation package
#sys.path.insert(0, "/path/to/your/package_or_module")

#
# argparse definition
#
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
    description=textwrap.dedent('''\
    This is BASEveg. It runs nuemrical simulations combining a vegetation growth module and BASEMENT.
    --------------------------------
    
    The folder structure of the project should be:

    - proj-folder
        --input
            --input veg (JSON)
            --(ideally all the other inputs, see the preprocessing step)
        --01_flood
        --02_flood
        --etc...

    '''))
parser.add_argument('-f', metavar = 'input_path',type=str, help='This is the path to the folder where input files are located.')
parser.add_argument('-s', type=str, metavar = 'BM_path',help='This is the path to the folder where BASEMENT executable are located.')
parser.add_argument('-n', type=int, default=1, metavar = 'ncores',help='This is the number of cores BASEMENT will use for running simulations. It has a default value (1).')
#parser.add_argument('--pre-processing', help='This is the optional agrument to run the pre-processing of the discharge (NOT ACTIVE).')
args=parser.parse_args()

if len(sys.argv) <= 1:
    print('Missing input parameters. See the help page')
    exit(1)

n_cores = args.n

n_cycles = 0
# loop over numeber of events (counting from the numeber of folders existing)
for name in os.listdir(args.f):
    match = re.search(r'cycle', name)
    if match:
        n_cycles += 1

print('-')
print(f'Simulating {n_cycles} cycles.')
print('-')

# get absoulute path
fpath = os.path.abspath(args.f)

input_to_BM_setup = '-f model.json -o setup.h5'
input_to_BM_simulation ='-f simulation.json -r setup.h5 -o results.h5'
input_to_BM_simulation_multi ='-f simulation.json -r setup.h5 -o results.h5 -n '+str(n_cores)
input_to_BM_results = '-f results.json -r results.h5 -o output'

BM_path_setup = os.path.join(args.s,'BMv3_BASEplane_setup')
BM_path_simulation = os.path.join(args.s,'BMv3_BASEplane_seq')
BM_path_simulation_multi = os.path.join(args.s,'BMv3_BASEplane_omp')
BM_path_results = os.path.join(args.s,'BMv3_BASEplane_results')

sim_code = 'all' # you can potentially run just some steps (e.g. to skip the 01_hq if you have done it already)

print('-')
print(f'Running {sim_code}')
print('-')

# set the time
start_time0 = time.time()

for event in range(0,n_cycles):

    start_time = time.time()

    # define the paths needed

    path_01 = os.path.join(fpath,'0'+str(event+1)+'_cycle','01_hq')
    path_02 = os.path.join(fpath,'0'+str(event+1)+'_cycle','02_veg_growth')
    path_03 = os.path.join(fpath,'0'+str(event+1)+'_cycle','03_init')
    path_04 = os.path.join(fpath,'0'+str(event+1)+'_cycle','04_flood')

    if event>0:
        copyfile('results.h5',os.path.join(path_01,'results.h5'))

    # change wd
    os.chdir(path_01)

    # Run BM_v3 with step-wise discharge
    if sim_code == 'all' or sim_code == 'hq':
        # run BM
        run_basement(BM_path_setup,BM_path_simulation,BM_path_simulation_multi,BM_path_results,
            input_to_BM_setup,input_to_BM_simulation,input_to_BM_simulation_multi,input_to_BM_results,n_cores)
    # copy .h5
    copyfile('results.h5',os.path.join(path_02,'results.h5'))
    # copy flow discharge
    #copyfile('inflow_discharge.txt',os.path.join(path_02,'inflow_discharge_steps.txt'))
    # copy results.json 
    copyfile('results.json',os.path.join(path_02,'results.json'))

    # change wd
    os.chdir(path_02)
    if sim_code == 'all' or sim_code == 'veg':
        # Run vegetation growth python script
        main_veg_module(path_02)
        # get the xdmf file
        subprocess.run(BM_path_results+' '+input_to_BM_results, shell=True)
    # copy .h5
    copyfile('results.h5',os.path.join(path_03,'results.h5'))

    # change wd
    os.chdir(path_03)
    # Run BM_v3 with constant discharge ( HYDRAULICS + VEGETATION )
    if sim_code == 'all' or sim_code == 'init':
        # run BM
        run_basement(BM_path_setup,BM_path_simulation,BM_path_simulation_multi,BM_path_results,
            input_to_BM_setup,input_to_BM_simulation,input_to_BM_simulation_multi,input_to_BM_results,n_cores)
    # copy .h5
    copyfile('results.h5',os.path.join(path_04,'results.h5'))

    # Run BM_v3 with unsteady discharge and bedload transport ( HYDRAULICS + VEGETATION + MORPHOLOGY )
    os.chdir(path_04)
    if sim_code == 'all' or sim_code == 'flood':
        # run BM
        run_basement(BM_path_setup,BM_path_simulation,BM_path_simulation_multi,BM_path_results,
            input_to_BM_setup,input_to_BM_simulation,input_to_BM_simulation_multi,input_to_BM_results,n_cores)
    
    print(' Cyce done in %s sec'% (time.time() - start_time))

print(' Simulation done in %s sec'% (time.time() - start_time0))