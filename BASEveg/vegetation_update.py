# Vegetation growth model, Caponi, 2020

from classes.mesh import *
from classes.discharges import TIMESERIE
from functions.input_output import *
from functions.math import *
from procedures import *

import os
import numpy as np
import csv
import time

#from numpy.polynomial import Polynomial as P
import warnings

def main_veg_module(path):
    """
    This is the main script for calculating the vegetation growth from BM3 results
    """

    # set the time
    start_time = time.time()
    start_time0 = start_time

    warnings.simplefilter('ignore', np.RankWarning)

    print('')
    print('Start vegetation growth calculation')
    print('')

    # get path to the input parameters (and output)
    #path = read_user_input(userinput,'f:')

    # Create subfolder for additional results printing
    pathout = os.path.join(path,'additional_output')
    if not os.path.exists(pathout):
        os.makedirs(pathout)

    # Read discharge time series - step-discharge - 
    qsim = read_float_values(os.path.join(path,'inflow_discharge_steps.txt'),1)
    qsim = list(dict.fromkeys(qsim)) # get unique values in order

    # Read input file including parameters
    parameters = read_json(os.path.join(path,'input_parameters_vegetation.json'))
    # --> access parameters like parameters['water_table'][0]['min_n_cells']

    # Read and allocate data from results.h5
    ntimesteps_toread = len(qsim) # I read only the timesteps in this cycle
    #print(f'ntimesteps to read is {ntimesteps_toread}')
    data = read_h5_and_get_data(os.path.join(path,'results.h5'),ntimesteps_toread)
    # --> list all keys of the dictionary print(list(data))

    # Read discharge time series file - dates -
    #formatTime = '%Y-%m-%d %H:%M:%S '
    #formattoMatch = '\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} '
    formatTime = '%d.%m.%Y %H:%M:%S'
    formattoMatch = '\d{2}.\d{2}.\d{4} \d{2}:\d{2}:\d{2}' 
    q = read_float_values(os.path.join(path,'inflow_discharge_date.txt'),2)
    dates = read_datetime(os.path.join(path,'inflow_discharge_date.txt'),formattoMatch,formatTime)
    duration = dates[-1]-dates[0]

    print(' 1. Read data ... %s sec'% (time.time() - start_time))
    #print("--- %s seconds ---" % (time.time() - start_time))
    start_time = time.time()

    # Assign NODE class
    nodes = []
    for ii in range(0,len(data['nodeCoords'])):
        nodes.append( NODE(data['nodeCoords'][ii][0],data['nodeCoords'][ii][1],0.0 ))

    # Assign to CELL class
    cells = []
    ii = 0
    for item in data['topology']:
        cells.append( CELL( nodes[item[0]], nodes[item[1]], nodes[item[2]]) )
        cells[-1].calcProperties("average")
        cells[-1].setElevation( data['cell_bottomEl'][ii] )
        ii += 1

    # Assign water table data to cells
    calculate_waterTable(cells,data,parameters,pathout,1)
    print(' 2. Calculation of the water table level ... %s sec'% (time.time() - start_time))
    start_time = time.time()

    # Get HQ relationship
    calculate_HQ(cells, parameters,qsim)
    print(' 3. Interpolation of the HQ relation ... %s sec'% (time.time() - start_time))
    start_time = time.time()

    if parameters['growth'][0]['type'] == 'root':
        # parameter calibration for roots
        calibration_root_model(cells,parameters,TIMESERIE(dates,q),1)
        # this also defines in which cells vegetation can grow (for the root model)
        print(' 4. Calibration of the root model.......DONE')
    else:
        print(' 4. No root model required.......')

    # Init vegetation on cells
    inizialize_vegetation(cells,data,parameters)
    print(' 5. Vegetation state initialization.......DONE')

    # Define where vegetation growth is not calculated, based on the minimum water table level
    calc_not_growing_cells(cells,parameters,TIMESERIE(dates,q))
    print(' 6. Computing cells where vegetation growth is not possible.......DONE')

    # exponential growth of the rooting depth -> needed as input for the root distribution
    advance_root_depth(cells, parameters,TIMESERIE(dates,q))
    print(' 7. Root depth advance in time.......DONE')

    if parameters['growth'][0]['type'] == 'root':
        # root distribution
        calc_root_distribution(cells,parameters,TIMESERIE(dates,q))
        print(' 8. Calc root distribution.......DONE')
    elif parameters['growth'][0]['type'] == 'no_root':
        # no_root model
        calc_veg_growth_parameters(cells,parameters,TIMESERIE(dates,q))
        print(' 8. Calc vegetation growth parameters.......DONE')
    else:
        # error 
        print('----Please enter the veg growth model - no_root - root')

    # Logistic growth or decay
    advance_biomass(cells,parameters,TIMESERIE(dates,q))
    print(' 9. Vegetation biomass advance in time........DONE')

    # write log and txt output
    #log_output(cells,parameters)
    write_txt_output(cells,parameters,TIMESERIE(dates,q),pathout)

    # modify data into the h5 file
    write_vegetation_on_h5(cells,os.path.join(path,'results.h5'))
    print(' 10. Vegetation data written........DONE')

    # modify data into the h5 file
    #write_watertable_on_h5(cells,os.path.join(path,'results.h5'))
    #print(' 11. Water table interpolation data written........DONE')

    print('')
    print('Vegetation growth calculation ends here. It took %s seconds in total.'%(time.time() - start_time0))
    print('')

if __name__ == '__main__':

    # parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        ----
        This runs (only) the vegetation growth module of BASEveg.
        ----
        It takes as input:
            - inflow_discharge_date
            - inflow_discharge_step
            - input parameters (JSON)
            - results.h5
        '''))
    parser.add_argument('-f',help='Folder paths where the input files are located.')
    args = parser.parse_args()

    if len(sys.argv) <= 1:
        print('Please see the help page, you are missing someinputs')
        exit(1)
    
    #path = read_user_inputs(sys.argv[1:],'f:',[]) 
    print('')
    print('#### PROGRESS ####')
    print('')
    main_veg_module(args.f)