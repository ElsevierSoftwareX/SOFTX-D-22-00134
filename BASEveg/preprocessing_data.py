from classes.discharges import *
from functions.input_output import *
from shutil import copy

import os
import numpy as np
import csv
import sys
import argparse
import textwrap
#import subprocess

def check_folders(simpath):
    """
    Check if folders and files exist
    """
    if not os.path.exists(simpath):
        print('No input files for discharge found.')

def creating_folders(outputpath,nevents):
    """
    This create and manage the folder structure within the simpath for a given number of events
    """
    # Create folders
    for i in range(1,nevents):
        if i>= 10:
            endpath = str(i)+'_cycle'
        else:   
            endpath = '0'+str(i)+'_cycle'
        #join the path
        name = os.path.join(outputpath,endpath)
        #create the folder
        if not os.path.exists(name):
            os.mkdir(name)
            # creating subfolders
            os.mkdir(os.path.join(name,'01_hq'))
            os.mkdir(os.path.join(name,'02_veg_growth'))
            os.mkdir(os.path.join(name,'03_init'))
            os.mkdir(os.path.join(name,'04_flood'))

def moving_files(inputpath,outputpath,nevent):
    """
    This function moves files from the input path to the destination path
    """
    cont = 1
    for i in range(1,nevent):
        
        # 01_hq
        copy( os.path.join(inputpath,'inflow_steps','inflow_discharge_steps'+str(i)+'.txt') ,
            os.path.join(outputpath,'0'+str(i)+'_cycle','01_hq','inflow_discharge.txt') )
        #os.rename('inflow_discharge_steps'+str(i)+'.txt','inflow_discharge_steps.txt')
        # 03_init
        copy( os.path.join(inputpath,'inflow_const','inflow_discharge_const'+str(i)+'.txt') ,
            os.path.join(outputpath,'0'+str(i)+'_cycle','03_init','inflow_discharge.txt') )
        # 04_flood
        copy( os.path.join(inputpath,'inflow_flood','inflow_discharge_flood'+str(i)+'.txt') ,
            os.path.join(outputpath,'0'+str(i)+'_cycle','04_flood','inflow_discharge.txt') )
        # 02_veg growth (json+inflows)
        copy( os.path.join(inputpath,'date_discharge','inflow_discharge_date'+str(cont)+'.txt') ,
            os.path.join(outputpath,'0'+str(i)+'_cycle','02_veg_growth','inflow_discharge_date.txt') )
        copy( os.path.join(inputpath,'inflow_steps','inflow_discharge_steps'+str(i)+'.txt') ,
            os.path.join(outputpath,'0'+str(i)+'_cycle','02_veg_growth','inflow_discharge_steps.txt') )
        copy( os.path.join(inputpath,'input_parameters_vegetation.json') ,
            os.path.join(outputpath,'0'+str(i)+'_cycle','02_veg_growth') ) 
        cont = cont + 2
        # basement files (mesh+json) 
        for file in os.listdir(inputpath):
            if file.endswith(".2dm"):
                meshfile = os.path.join(inputpath, file)
                break
        pathgroup = [os.path.join(outputpath,'0'+str(i)+'_cycle','01_hq'),os.path.join(outputpath,'0'+str(i)+'_cycle','03_init'),
            os.path.join(outputpath,'0'+str(i)+'_cycle','04_flood')]
        for path in pathgroup:
            copy( meshfile, path )
            copy( os.path.join(inputpath,'model.json'), path )
            copy( os.path.join(inputpath,'results.json'), path )
            copy( os.path.join(inputpath,'simulation.json'), path )

def check_model_JSON(outputpath,nevent):
    """
    This checks the model.json file, if it has all necessary block and filenames right for BASEveg sim
    """
    for event in range(1,nevent):

        # define the paths needed
        path_01 = os.path.join(outputpath,'0'+str(event)+'_cycle','01_hq')
        #path_02 = os.path.join(outputpath,'0'+str(event+1)+'_cycle','02_veg_growth')
        path_03 = os.path.join(outputpath,'0'+str(event)+'_cycle','03_init')
        path_04 = os.path.join(outputpath,'0'+str(event)+'_cycle','04_flood')

        # check filename inputs
        modify_json_INFLOW(os.path.join(path_01,'model.json'),'inflow_discharge.txt')
        modify_json_INFLOW(os.path.join(path_03,'model.json'),'inflow_discharge.txt')
        modify_json_INFLOW(os.path.join(path_04,'model.json'),'inflow_discharge.txt')
        """for file in os.listdir(path_04):
            if file.endswith(".txt"):
                namefileinflow = os.path.join(path_01, file)"""
        # 01
        times = read_float_values(os.path.join(path_01,'inflow_discharge.txt'),0)
        modify_json_TIME(os.path.join(path_01,'simulation.json'),times)
        modify_json_BEDSTART(os.path.join(path_01,'model.json'),times[-1]+1)
        if event>1:
            modify_json_INITIAL(os.path.join(path_01,'model.json'),times)
            modify_json_INITIAL(os.path.join(path_03,'model.json'),times)
            modify_json_INITIAL(os.path.join(path_04,'model.json'),times)
        # 03
        times = read_float_values(os.path.join(path_03,'inflow_discharge.txt'),0)
        modify_json_TIME(os.path.join(path_03,'simulation.json'),times)
        modify_json_BEDSTART(os.path.join(path_03,'model.json'),times[-1]+1)
        # 04
        times = read_float_values(os.path.join(path_04,'inflow_discharge.txt'),0)
        modify_json_TIME(os.path.join(path_04,'simulation.json'),times)
        


if __name__ == '__main__':
    # parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        ----
        This organizes the folder structure and the input parameters to run simulations with BASEveg.
        ----
        Within the input folder there should be:
        - JSON files for BM simulations
        - mesh
        - JSON file for vegetation
        - 4 folders created by the preprocessing step
            - date_discharge
            - inflow_const
            - inflow_flood
            - inflow_steps
        '''))
    parser.add_argument('-f',metavar=('input_path','destination_path'),nargs=2,help='Folder paths where the input files (from the preprocessing) are located and simulation files will be stored.')
    parser.add_argument('ncycle',nargs=1, type=int, help='Number of cycle events to simulate.')
    args = parser.parse_args()
    #print(args.f[0])  
    # create folder structure
    creating_folders(args.f[1],args.ncycle[0]+1)
    # move files
    moving_files(args.f[0],args.f[1],args.ncycle[0]+1)
    # small check
    check_model_JSON(args.f[1],args.ncycle[0]+1)