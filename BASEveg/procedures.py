# Collection of procedures (involves the use of functions) to be used in the vegetation growth model
# F. Caponi, 2020

import numpy as np
from functions.functions import *
from functions.input_output import terminal_output
from scipy.optimize import curve_fit
from scipy.special import gamma, factorial, gammainc
import subprocess
import os

def run_basement(BM_path_setup,BM_path_simulation,BM_path_simulation_multi,BM_path_results,
    input_to_BM_setup,input_to_BM_simulation,input_to_BM_simulation_multi,input_to_BM_results,n_cores):
    """
    This function runs BASEMENT in sinigle or multicores
    """
    subprocess.run(BM_path_setup+' '+input_to_BM_setup, shell=True)
    if n_cores==1:
        subprocess.run(BM_path_simulation+' '+input_to_BM_simulation, shell=True)
    else:
        subprocess.run(BM_path_simulation_multi+' '+input_to_BM_simulation_multi, shell=True)
    subprocess.run(BM_path_results+' '+input_to_BM_results, shell=True)

def inizialize_vegetation(cells,data,parameters):
    """
    This allocate the vegetation state on the CELL class
    """
    ii = 0
    for item in cells:
        # aboveground biomass
        item.setBc( data['cell_vegstate'][ii][0] )
        # belowground biomass
        item.setBr( data['cell_vegstate'][ii][1] )
        # initial value (to be stored)
        item.setB0( data['cell_vegstate'][ii][0]+data['cell_vegstate'][ii][1] )
        # rooting depth
        item.setD( data['cell_vegstate'][ii][2] )
        # initial value (to be stored)
        item.setD0( data['cell_vegstate'][ii][2] )
        ii += 1

def calc_not_growing_cells(cells,parameters,timeseries):
    """
    This calculates where vegetation cannot growth
    """
    for item in cells:
        # Define where vegetaiton can grow (in which cells)
        if parameters['growth'][0]['type'] == 'root':
            alpha = item.getParameters()[0]
            if alpha >= 0.0 and alpha <=1.0:
                item.setIsgrowing( True )
        elif parameters['growth'][0]['type'] == 'no_root':
            p = item.getPol()
            qmin = timeseries.getMin()
            wlmin = p(qmin)
            if item.getElevation() > wlmin: # exclude cells below the minimum water table level
                item.setIsgrowing( True )
        else:
            pass

def calc_veg_growth_parameters(cells,parameters,timeseries):
    """
    this calculates the vegetation growth parameters in case "no_root" model is selected
    """
    for item in cells:
        # get water level time series
        p = item.getPol()
        wlmin = p( timeseries.getMin() )
        wl = p( timeseries.getValue() ) #np array
        # otpimal box
        root_low_limit = item.getElevation() - item.getD()
        box_upper_limit = wlmin + parameters['growth'][0]['logistic'][0]['g2']
        box_lower_limit = wlmin + parameters['growth'][0]['logistic'][0]['g1']
        decay_upper_limit = wlmin + parameters['growth'][0]['logistic'][0]['d2']
        decay_lower_limit = wlmin + parameters['growth'][0]['logistic'][0]['d1']
        contg = 0
        contd = 0
        for ww in wl:
            if (root_low_limit < ww < item.getElevation() ) and (box_lower_limit < ww < box_upper_limit):
                contg += 1
            elif (ww < decay_lower_limit ) or (ww > decay_upper_limit ):
                contd += 1
        item.setNETgrowth(contg)
        # and how many in the decay zone
        item.setNETdecay(contd)


def calc_root_distribution(cells,parameters,timeseries):
    """
    This calculates the root distribution according to Tron's model, given the parameters
    h1 -> input as meters above the minimum water level (m)
    h2 -> minimum water level (m asl)
    L -> input as the width of the cappilary fringe (m)
    """
    # Create vector of soil depths where to calc the root profile
    #soil_ = np.arange(0,1,parameters['root_model'][0]['deltaz'])
    soil_ = np.linspace(0.0, 1.0, num=int(1/parameters['root_model'][0]['deltaz']))
    qmin = timeseries.getMin()
    for ee in cells:
        # only cells above h2
        p = ee.getPol()
        h2 = p( qmin )
        if ee.getIsgrowing():
            # get parameters in the right form
            h1 = h2 + parameters['root_model'][0]['h1']
            h1dim = dimensionless(h1,h2,ee.getElevation()[0])
            DSoil = ee.getElevation()[0]-h2
            L = parameters['root_model'][0]['L'] 
            Ldim = L/DSoil 
            #Ldim = h1dim - ((ee.getElevation()[0]-L)/DSoil)
            Ddim = ee.getD()/DSoil #dimensionless rooting depth
            alpha = ee.getParameters()[0]
            le = ee.getParameters()[1]

            # --Some corrections--
            # Approx of le, to avoid negative value
            if le<=1.5 and le>1:
                le = 1.0
                pass
            elif le>1.5 and le<2:
                le = 2.0

            if Ddim>h1dim:
                Ddim = h1dim# limit rooting depth to the maximum available, should be never the case
            elif Ldim>h1dim:
                Ldim=h1dim# maybe not necessary, should be never the case

            #print(f'h1dim {h1dim}')
            #print(f'h2 {h2}')
            #print(f'Ldim {Ldim}')
            #print(f'DSoil {DSoil}')

            k = []
            root_dist = []
            theta = []
            for s in soil_:
                # 1. define until which depth roots can grow
                if s > Ddim:
                    theta.append(0.0)
                else:
                    theta.append(1.0)
                # 2. calc k, Eq. 5, Tron et al., 2014
                if s < (h1dim-Ldim):
                    k.append( (gammainc((h1dim-s-Ldim)/alpha,le)-gammainc((h1dim-s)/alpha,le))/gamma(le) )
                elif (s < h1dim) and (s>h1dim-Ldim):
                    k.append( 1-(gammainc( (h1dim-s)/alpha,le)/gamma(le)) )
                else:
                    k.append(0.0)
                # 3. calc root profile, Eq. 10, Tron et al. 2014
                if (theta[-1]+theta[-1]*k[-1]+1-k[-1])==0.0:
                    print('Division by zero in the calculation of theta in the root model.')
                    root_dist.append( 0.0 )
                else:
                    root_dist.append( (2*theta[-1]*k[-1])/(theta[-1]+theta[-1]*k[-1]+1-k[-1]) )
            # save the root profile in the CELL
            ee.setRootDist(root_dist)
        else:
            fact = 1/parameters['root_model'][0]['deltaz']
            ee.setRootDist([0.0]*int(fact))

def calibration_water_table_pdf(le,zw_mean_dim,code):
    """
    This allows to calibrate the alpha and le values needed to calculate the root distribution
    """
    if code == 0:
        # calibration of the mean value, freqeuncy (le) constant
        tmp = (1- zw_mean_dim)/le
        return tmp
    else:
        # full model calibration -> to be done
        pass

def calculate_waterTable(cells,data,parameters,pathout,code):
    """
    This calculates the water table based on a step-discharge simulation (and save the data).
    """
    # run the IDW interpolator -> NEED TO RUN IN PARALLEL (to be done)
    ntimestep_tointerp = len(data['cell_waterDepth'])
    if code == 1:
        cells_coords = np.empty([data['ncells'],3])
        for ii in range(0,len(cells)):
            # x y z
            cells_coords[ii][0] = cells[ii].getCenter()[0] 
            cells_coords[ii][1] = cells[ii].getCenter()[1]
            cells_coords[ii][2] = cells[ii].getElevation()

        #startt = data['ntimestep']-ntimestep_tointerp
        # the first timestep (the initial one) is discharged here!!
        cells_waterTable = np.empty([ntimestep_tointerp,data['ncells']])
        for cont in range(0,ntimestep_tointerp):
            cells_waterTable[cont] = IDW_interpolator(
                cells_coords, 
                data['cell_waterDepth'][cont], 
                data['cell_wse'][cont], 
                data['min_waterDepth'][0], 
                parameters['water_table'][0]['min_n_cells'], 
                parameters['water_table'][0]['min_dist'], 
                parameters['water_table'][0]['exp_dist'], 
                parameters['water_table'][0]['increase_dist'])

        #np.savetxt(os.path.join(pathout,'water_table.txt'), cells_waterTable, fmt='%.2f')
        mat = np.matrix(cells_waterTable)
        with open(os.path.join(pathout,'water_table.txt'),'wb') as f:
            for line in mat:
                np.savetxt(f, line, fmt='%.2f')
    else:
        cells_waterTable = np.loadtxt(os.path.join(pathout,'water_table.txt'))
    # save data on CELL
    for i in range(0,len(cells)):
        tmplist = []
        for q in range(0,ntimestep_tointerp):
            tmplist.append(cells_waterTable[q][i])
        cells[i].setWaterTable(tmplist)

def calculate_HQ(cells, parameters,qsim):
    """
    This calculates the coefficients of the polynomial for fitting an HQ relation
    """
    for item in cells:
        tmpwt = item.getWaterTable()
        degree = parameters['discharges'][0]['poly_degree']
        if len(tmpwt)>degree:
            p = polyFit(qsim,tmpwt,parameters['discharges'][0]['poly_degree'])
            item.setPol(p)
        else:
            print('No H-Q calibration. The polynomial degree is greater than the actual water table data.')

def calibration_root_model(cells,parameters,timeseries,code):
    """
    This calibrates the parameters (alpha, and le) needed to compute the root profile according to the Tron et al., 2014 model
    code = 1 -> easy calibration (mean)
    code = 0 -> default, auto calibration
    """
    for item in cells:
        if item.getPol(): # in case I have a HQ relationship
            p = item.getPol()
            qmin = timeseries.getMin()
            qmean = timeseries.getMean()
            # calculation are performed only for cells that are above the mwl
            if (item.getElevation()[0] > p(qmin) ):
                if code:
                    # calibration (now is only if Q is const)
                    zw_mean_dim = dimensionless(qmean,qmin, item.getElevation()[0])
                    alpha = calibration_water_table_pdf(parameters['root_model'][0]['le'],zw_mean_dim,
                        parameters['root_model'][0]['calibration'])
                    if alpha > 1.0:
                        item.setParameters( [-1.0, -1.0] ) #I do not consider points below the mean level
                    else:
                        item.setParameters( [alpha, parameters['root_model'][0]['le']] )
                else:
                    pass

def advance_root_depth(cells,parameters,timeseries):
    """
    This advances in time the rooting depth through an exponential growth
    """
    qmin = timeseries.getMin()
    duration = timeseries.calc_duration_days()
    print(f' 7.1. Roots grow for {duration} days')
    for item in cells:
        #get parameters
        if item.getIsgrowing():
            p = item.getPol()
            h2 = p( qmin )
            max_rootdepth = item.getElevation() - h2
            sigma = parameters['growth'][0]['exponential'][0]['growth_rate']
            deltat = parameters['growth'][0]['deltat']
            #duration = parameters['growth'][0]['duration']
            #advance solution in time
            newvalue = advance_in_time(item.getD(), max_rootdepth, duration, deltat, 'exp', sigma)
            item.setD( newvalue[0] )

def advance_biomass(cells,parameters,timeseries):
    """
    This advances in time the biomass through a logistic growth
    """
    # duration of growth
    #if parameters['growth'][0]['type'] == 'root':
    duration = timeseries.calc_duration_days()
    deltat_timeseries = timeseries.calc_deltat()# seconds
    #get parameters
    lambda_bc = parameters['growth'][0]['logistic'][0]['prop_above']
    lambda_br = parameters['growth'][0]['logistic'][0]['prop_below']
    # check
    if lambda_bc+lambda_br>1:
        print('WARNING: the sum of the proportion of biomass above and below should be at max 1!')
    deltat = parameters['growth'][0]['deltat']
    max_biomass = parameters['growth'][0]['logistic'][0]['max_biomass']
    # output to terminal
    print(f' 8.1. Days available for vegetation growth = {duration} days.')
    #print(f'The deltat in the discharge time series (in seconds) is {deltat_timeseries}')
    for item in cells:
        # get actual values of biomass
        biomass_above = item.getBc()
        biomass_below = item.getBr()
        if item.getIsgrowing():
            # calc growth rate of the biomass
            if biomass_above+biomass_below == 0.0:
                # increase by a small fraction, otherwise logistic growth cannot be caluclated
                biomass_above = 0.001*lambda_bc
                biomass_below = 0.001*lambda_br
            if parameters['growth'][0]['type'] == 'no_root':
                duration = (item.getNETgrowth()*deltat_timeseries)/(3600*24) # days
                sigma = parameters['growth'][0]['logistic'][0]['growth_rate']
            elif parameters['growth'][0]['type'] == 'root':
                factor_sigma = parameters['growth'][0]['logistic'][0]['factor_growth_rate']
                mean_rootbiomass = calc_mean_root_biomass(item.getRootDist(),parameters['root_model'][0]['deltaz'])
                sigma = mean_rootbiomass*factor_sigma
                item.setNETgrowth(sigma)# in case of root model this is used as sigma
            else:
                sigma = 0.0

            # advance solution in time
            newvalue = advance_in_time(biomass_above+biomass_below, max_biomass, duration, deltat, 'logistic', sigma)
            deltabiomass = newvalue - (biomass_above+biomass_below)
            # allocate the solution
            item.setBc( biomass_above + deltabiomass*lambda_bc )
            item.setBr( biomass_below + deltabiomass*lambda_br )
            item.setB(item.getBc()+item.getBr())
        
        else: # decay
            if parameters['growth'][0]['type'] == 'no_root':
                if (item.getBc()+item.getBr() > 0.0):
                    duration = (item.getNETdecay()*deltat_timeseries)/(3600*24) # days
                    sigma = parameters['growth'][0]['logistic'][0]['decay_rate']
                    # advance solution in time
                    newvalue = advance_in_time(biomass_above+biomass_below, max_biomass, duration, deltat, 'exp_decay', sigma)
                    deltabiomass = newvalue - (biomass_above+biomass_below)
                    # allocate the solution (to be checked)
                    if ((biomass_above + deltabiomass*lambda_bc) > 0.0):
                        item.setBc( biomass_above + deltabiomass*lambda_bc )
                    else:
                        item.setBc( 0.0 )
                    if ((biomass_below + deltabiomass*lambda_br) > 0.0):
                        item.setBr( biomass_below + deltabiomass*lambda_br )
                    else:
                        item.setBr( 0.0 )
                    #item.setBr( biomass_below + deltabiomass*lambda_br )
                    item.setB(item.getBc()+item.getBr())
