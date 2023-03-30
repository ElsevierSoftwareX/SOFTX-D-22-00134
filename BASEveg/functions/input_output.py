import numpy as np
import json
import h5py
import sys
import getopt
import re
from datetime import datetime, timedelta, date
import subprocess
import os
from tabulate import tabulate

def read_float_values(filepath,ncolumn):
    """
    This reads a specific column of float values from a file
    """
    values = []
    fid = open(filepath,'r')
    for line in fid:
        values.append( float(line.split()[ncolumn]) )
    fid.close()
    return values

def read_datetime(filepath,form,formTime):
    """
    This searches and read for dates in specific format from a file
    """
    times = []
    fid = open(filepath,'r')
    for line in fid:
        matching = re.search(form, line)
        times.append( datetime.strptime(matching.group(), formTime) )
    fid.close()
    return times

def get_timeseries(filepath):
    """
    This function return a timeseries object
    """
    formatTime = '%d.%m.%Y %H:%M:%S '
    formattoMatch = '\d{2}.\d{2}.\d{4} \d{2}:\d{2}:\d{2} '
    fid = open(filepath,'r')
    for line in fid:
        matching = re.search(formattoMatch, line)	
        times.append( datetime.strptime(matching.group(), formatTime) )
    fid.close()
    return times

def read_user_inputs(argv,flags,modes):
    """
    This reads the user input with the optget module and return a dictionary. 
    It raises errors if flags (should be followed by :) are not given
    """
    try: 
        opts, args = getopt.getopt(argv, flags, modes) 
    except getopt.GetoptError as err:
        print(err)    
        #usage()
        sys.exit(2)
 
    #i = re.split(r':',flags)
    listofinputs = {}
    for opt, arg in opts: 
        if opt in ['-f']: 
            listofinputs['input'] = arg
        if opt in ['-s']:
            listofinputs['source'] = arg
        if opt in ['-o']:
            listofinputs['output'] = arg
        if opt in ['-n']:
            listofinputs['n_cores'] = arg
        if opt in ['--pre-processing']:
            listofinputs['pre-proc'] = True 
    if args:
        listofinputs['arg'] = args[0]
    return listofinputs

def read_user_input(argv):
    """
    This reads the user input with the optget module
    """
    try: 
        opts, args = getopt.getopt(argv, "f:") 
    except: 
        print("Error. Please insert with the falg -f the path to the input files.") 

    for opt, arg in opts: 
        if opt in ['-f']: 
            pathname = arg
    return pathname

def read_json(file_name):
    """
    This reads a .json file and give you a dict
    """
    with open(file_name,"r") as f:
        your_dict = json.load(f)
    return your_dict

def write_json(file_name,dict):
    """
    This writes a .json file given a dict and file_name
    """
    with open(file_name,"w") as f:
        json.dump(dict,f,indent=4, sort_keys=True)

def modify_json_TIME(file_name,times):
    """
    This modifies the TIME part in the .json file in simulation.json 
    """
    data = read_json(file_name)
    data['SIMULATION']['TIME']['start'] = times[0]
    data['SIMULATION']['TIME']['end'] = times[-1]
    data['SIMULATION']['TIME']['out'] = times[1]-times[0]
    write_json(file_name,data)

def modify_json_INITIAL(file_name,times):
    """
    This modifies the INITIAL part in the .json file in model.json 
    """
    data = read_json(file_name)
    if data['SETUP']['DOMAIN']['BASEPLANE_2D']['HYDRAULICS']['INITIAL']['type'] == 'dry':
        data['SETUP']['DOMAIN']['BASEPLANE_2D']['HYDRAULICS']['INITIAL']['type'] = 'continue'
        data['SETUP']['DOMAIN']['BASEPLANE_2D']['HYDRAULICS']['INITIAL']['file'] = 'results.h5'
        data['SETUP']['DOMAIN']['BASEPLANE_2D']['HYDRAULICS']['INITIAL']['time'] = times[0]

    if data['SETUP']['DOMAIN']['BASEPLANE_2D']['MORPHOLOGY']['INITIAL']['type'] == 'mesh':
        data['SETUP']['DOMAIN']['BASEPLANE_2D']['MORPHOLOGY']['INITIAL']['type'] = 'continue'
        data['SETUP']['DOMAIN']['BASEPLANE_2D']['MORPHOLOGY']['INITIAL']['file'] = 'results.h5'
        data['SETUP']['DOMAIN']['BASEPLANE_2D']['MORPHOLOGY']['INITIAL']['time'] = times[0]

    if data['SETUP']['DOMAIN']['BASEPLANE_2D']['VEGETATION']['INITIAL']['type'] == 'region_defined':
        data['SETUP']['DOMAIN']['BASEPLANE_2D']['VEGETATION']['INITIAL']['type'] = 'continue'
        data['SETUP']['DOMAIN']['BASEPLANE_2D']['VEGETATION']['INITIAL']['file'] = 'results.h5'
        data['SETUP']['DOMAIN']['BASEPLANE_2D']['VEGETATION']['INITIAL']['time'] = times[0]

    write_json(file_name,data)

def modify_json_BEDSTART(file_name,time):
    """
    This modifies the time at which bed load starts 
    """
    data = read_json(file_name)
    data['SETUP']['DOMAIN']['BASEPLANE_2D']['MORPHOLOGY']['PARAMETER']['morphodynamic_start'] = time
    write_json(file_name,data)

def modify_json_INFLOW(file_name,new_name):
    """
    This modifies the the file name of the inflow in model.json 
    """
    data = read_json(file_name)
    data['SETUP']['DOMAIN']['BASEPLANE_2D']['HYDRAULICS']['BOUNDARY']['STANDARD'][0]['discharge_file'] = new_name
    write_json(file_name,data)

def read_h5_and_get_data(file_name,ntimesteps_toread):
    """
    This reads the results.h5 file and returns needed data
    """
    with h5py.File(file_name, "r") as f:
        # get output times available
        timestep = np.array(f["RESULTS"]["CellsAll"]["HydState"])
        ntimestep = len(timestep) # total number of timesteps
        # number of cells
        data0 = np.array(f["RESULTS"]["CellsAll"]["HydState"][timestep[0]])
        ncells = len(data0)
        # min water depth
        min_waterDepth = np.array(f["Parameters"]["MinWaterDepth"])
        # get coordinates
        nodeCoords = np.array(f["NodesAll"]["Coordnts"])
        topology = np.array(f["CellsAll"]["Topology"])
        #get vegetation (at the last time step)
        vegdata = np.array(f["RESULTS"]["CellsAll"]["VegState"][timestep[-1]])
        # get water depth
        cell_wse = np.empty([ntimesteps_toread,ncells])
        cell_waterDepth = np.empty([ntimesteps_toread,ncells])
        #cell_bottomEl = np.array(f["CellsAll"]["BottomEl"])
        cell_bottomEl = np.array(f["RESULTS"]["CellsAll"]["BottomEl"][timestep[-1]])
        cont = 0
        for t in range(ntimestep - ntimesteps_toread,ntimestep):
            #print(ntimestep)
            #print(ntimesteps_toread)
            data = np.array(f["RESULTS"]["CellsAll"]["HydState"][timestep[t]])
            for n in range(0,ncells):
                cell_wse[cont][n] =  data[n][0]
                cell_waterDepth[cont][n] = cell_wse[cont][n] - cell_bottomEl[n]
            cont += 1
    return dict([('cell_bottomEl',cell_bottomEl),
    ('cell_wse',cell_wse),
    ('cell_waterDepth',cell_waterDepth),
    ('nodeCoords',nodeCoords),
    ('topology',topology),
    ('min_waterDepth',min_waterDepth),
    ('timestep',timestep),
    ('ntimestep',ntimestep),
    ('ncells',ncells),
    ('cell_vegstate',vegdata)]) 


def modify_h5(file_name,new_data,name_data):
    """
    This opens and modifes the results.h5 file to inlcude new calculated data. Data should be already in the correct form.
    Writes on the last timestep available in the .h5 file.
    """
    with h5py.File(file_name, "r+") as f:
        # get output times available
        timestep = np.array(f["RESULTS"]["CellsAll"]["HydState"])
        ntimestep = len(timestep)
        # number of cells
        data0 = np.array(f["RESULTS"]["CellsAll"]["HydState"][timestep[0]])
        ncells = len(data0)
        # write on file (on the last timestep)
        for n in range(0,ncells):
            f["RESULTS"]["CellsAll"][name_data][timestep[-1]][n] = new_data[n]

def write_vegetation_on_h5(cells,file_name):
    """
    This writes on .h5 the vegstate.
    """
    #provide data in the correct form 
    vegstate = np.empty([ len(cells),3 ])
    i = 0
    for item in cells:
        vegstate[i][0] = item.getBc()
        vegstate[i][1] = item.getBr()
        vegstate[i][2] = item.getD() 
        i += 1
    # modify the h5 file
    modify_h5(file_name,vegstate,'VegState')

def write_watertable_on_h5(cells,file_name):
    """
    This writes on .h5 the water table.
    """
    #provide data in the correct form 
    nwt = len(cells[0].getWaterTable())
    #print("Number of time the water table is interpolated is:")
    # modify the h5 file
    with h5py.File(file_name, "r+") as f:
        # get output times available
        timestep = np.array(f["RESULTS"]["CellsAll"]["WaterTEl"])
        # number of cells
        data0 = np.array(f["RESULTS"]["CellsAll"]["WaterTEl"][timestep[0]])
        ncells = len(data0)
        # write on file (on the last timestep)
        for t in range(1,nwt+1):
            #print(-t)
            for n in range(0,ncells):
                f["RESULTS"]["CellsAll"]["WaterTEl"][timestep[-t]][n] = cells[n].getWaterTable()[-t]

def write_txt_output(cells,parameters,timeseries,pathout):
    """
    Writes text output for checking the results
    """
    with open(os.path.join(pathout,"xyz.txt"),"w") as f:
        for item in cells:
            f.write(f"{item.getCenter()[0]:10.4f}\t{item.getCenter()[1]:10.4f}\t{item.getElevation()[0]:10.4f}\n")

    with open(os.path.join(pathout,"cell_features.txt"),"w") as f:
        #f.write(f'zb [m] \t zw_min [m] \t zw_mean [m] \t alpha [-] \t le [-] \n')
        table = np.zeros((len(cells),5))
        i = 0
        for item in cells:
            p = item.getPol()
            table[i][0] = item.getElevation()
            table[i][1] = p( timeseries.getMin() )
            table[i][2] = p( timeseries.getMean() )
            table[i][3] = item.getParameters()[0]
            table[i][4] = item.getParameters()[1]
            i += 1
        f.write(tabulate(table,headers=["elevation [m]","elevation min_wl [m]", "elevation mean_wl [m]", 'alpha [-]','le [-]']))
        #f.write(f"{item.getElevation()[0]:10.4f} \t {hmin:10.4f} \t {hmean:10.4f} \t {item.getParameters()[0]:=.2f} \t {item.getParameters()[1]:=.2f}\n")

    with open(os.path.join(pathout,"growth_biomass.txt"),"w") as f:
        table = np.zeros((len(cells),8))
        i = 0
        for item in cells:
            table[i][0] = item.getB0()
            table[i][1] = item.getB()
            table[i][2] = parameters['growth'][0]['deltat']
            table[i][3] = (item.getNETgrowth()*timeseries.calc_deltat())/(3600*24)
            table[i][4] = (item.getNETdecay()*timeseries.calc_deltat())/(3600*24)
            if parameters['growth'][0]['type'] == 'no_root':
                table[i][5] = parameters['growth'][0]['logistic'][0]['growth_rate']
            elif parameters['growth'][0]['type'] == 'root':
                table[i][5] = item.getNETgrowth()
            table[i][6] = parameters['growth'][0]['logistic'][0]['decay_rate']
            table[i][7] = parameters['growth'][0]['logistic'][0]['max_biomass']
            i += 1
        f.write(tabulate(table,headers=["B0 [-]","B [-]", "deltat [day]", 'time_growth [day]','time_decay [day]','growth_rate [1/day]','decay_rate [1/day]','Bmax [-]']))

    with open(os.path.join(pathout,"growth_root.txt"),"w") as f:
        table = np.zeros((len(cells),5))
        i = 0
        for item in cells:
            p = item.getPol()
            table[i][0] = item.getD0()
            table[i][1] = item.getD()
            table[i][2] = parameters['growth'][0]['deltat']
            table[i][3] = parameters['growth'][0]['exponential'][0]['growth_rate']
            table[i][4] = item.getElevation() - p( timeseries.getMin() )
            i += 1
        f.write(tabulate(table,headers=["D0 [m]","D [m]", "deltat [day]",'growth_rate [1/day]','Dmax [m]']))

    with open(os.path.join(pathout,"p_coeff.txt"),"w") as f:
        for item in cells:
            l = item.getPol()
            #print(l.c)
            for ii in range(0,len(l.c)):
                f.write(f'{l.c[ii]:=1.3e}\t')
            f.write('\n')

    with open(os.path.join(pathout,"root_profile.txt"),"w") as f:
        for item in cells:
            rootprofile = item.getRootDist()
            for r in rootprofile:
                f.write(f'{r:3.6f}\t')
            f.write('\n')

    with open(os.path.join(pathout,"vegstate.txt"),"w") as f:
        for item in cells:
            bc = item.getBc()
            br = item.getBr()
            d = item.getD()
            f.write(f'{bc:.4f}\t{br:.4f}\t{d:.4f}\n')

def write_txt_statistics(events,pathfile):
    """
    This writes stistics of the preprocessing discharge step to a table
    """
    with open(os.path.join(pathfile,'statistics_cycles.txt'),"w") as f:
        table = np.zeros((len(events),5))
        i = 0
        for ee in events:
            table[i][0] = np.min(ee.getValue())
            table[i][1] = np.max(ee.getValue())
            table[i][2] = np.mean(ee.getValue())
            table[i][3] = ee.calc_duration_days_effective()*24
            table[i][4] = ee.calc_duration_days()*24
            i += 1
        f.write(tabulate(table,headers=["q_min [m3/s]","q_max [m3/s]", "q_mean [m3/s]",'duration_effective [hours]','duration_total [hours]']))

def write_txt_time_series(timeseries,pathfile,namefile):
    """
    This writes on file timeseries to be simulated in BASEMENT. It takes the object TIMESERIES as input
    """
    with open(os.path.join(pathfile,namefile+".txt"),"w") as f:
        t = timeseries.getTime()
        q = timeseries.getValue()
        for i in range(0,len(q)):
            f.write(f'{t[i]*3600}\t{q[i]:.1f}\n')

def write_txt_timenum_series(timeseries,pathfile,namefile):
    """
    This writes on file timeseries as time (numeral)-discharge.
    """
    with open(os.path.join(pathfile,namefile+".txt"),"w") as f:
        t = timeseries.getTime()
        q = timeseries.getValue()
        for i in range(0,len(q)):
            f.write(f'{t[i]}\t{q[i]:.1f}\n')  

def write_txt_date_series(timeseries,pathfile,namefile):
    """
    This writes on file timeseries as datetime-discharge.
    """
    with open(os.path.join(pathfile,namefile+".txt"),"w") as f:
        t = timeseries.getTime()
        q = timeseries.getValue()
        for i in range(0,len(q)):
            f.write(f'{t[i].strftime("%d.%m.%Y %H:%M:%S")}\t{q[i]:.1f}\n') 

def log_output(cells,parameters):
    """
    This write on screen some output to check the results
    """
    ncells_withveg = 0
    cells_tot = len(cells)
    for item in cells:
        if item.getParameters()[0] == -1.0:
            pass
        else:
            ncells_withveg += 1

    #print(f'  -----> The percentage of cells with vegetation is {100*ncells_withveg/cells_tot:3.2f} %')

def terminal_output(message):
    """
    This writes in the terminal whatever is passed as message
    """
    print('--------')
    print('')
    for mess in message:
        print(mess)
        print('')
    print('')
    print('--------')
