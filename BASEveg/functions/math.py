# Collection of functions to be used in the vegetation growth model
# F. Caponi, 2020

import numpy as np
from scipy import integrate

def IDW_interpolator(
    cell_coords, cell_waterdepth, cell_wse, min_waterDepth, min_nCellsToInterpolate, min_distance, p_exp, increase_distance):
    """
    This algorithm finds the value of an unknwon variable by interpolating
    the values of that variable at known locations within a certain distance.
    To average the values, it uses a IDW (inverse distance weighted) method.
    """
    # find all wet cells
    cell_wet = cell_wse[cell_waterdepth > min_waterDepth]   
    cell_wet_coord = cell_coords[cell_waterdepth > min_waterDepth]
    # init some variables
    cell_interpolatedValue = np.empty([len(cell_wse),1])
    # loop over all cells
    print("Progress water table interpolation ", flush=True)
    deltaprogress = len(cell_wse)/10
    progress = deltaprogress
    ct = 10
    for ii in range(0,len(cell_wse)):
        # progress bar
        if ii > progress:
            print(f'-{ct}%', end="",flush=True)
            ct += 10
            progress += deltaprogress
        # check if cell is wet
        if cell_waterdepth[ii] > min_waterDepth:
            # water table = water surface elevation
            cell_interpolatedValue[ii] = cell_wse[ii]
        else:
            # coordinates of the ii cell (just to be clearer)
            x_target = cell_coords[ii][0]
            y_target = cell_coords[ii][1]
            min_distance_ = min_distance # just a copy
            # loop until it finds at least min_nCellsToInterpolate points to interpolate
            cont = 0
            while cont <= min_nCellsToInterpolate:
                # some init
                dist_inv_sum = 0
                dist_inv = np.empty([0])
                data_toInterpolate = np.empty([0])
                # loop over wet cells
                for jj in range(0,len(cell_wet)):
                    # calculate distance between target and wet cell
                    dist = np.sqrt( np.power((x_target-cell_wet_coord[jj][0]),2) + np.power((y_target-cell_wet_coord[jj][1]),2) )
                    # check if cell is within the radius min_distance
                    if dist <= min_distance_:
                        # calc inverse distance
                        dist_inv = np.append(dist_inv, np.power(dist, -p_exp) )
                        dist_inv_sum += dist_inv[-1]
                        data_toInterpolate = np.append(data_toInterpolate, cell_wet[jj] )
                        cont += 1
                if cont < min_nCellsToInterpolate:
                    # increase searching distance
                    min_distance_ += increase_distance
                    cont = 0

            # calculate weights
            lambda_value = np.empty([len(dist_inv)])
            lambda_value = (dist_inv * data_toInterpolate) / dist_inv_sum
            # assign interpolated value
            cell_interpolatedValue[ii] = np.sum(lambda_value)
    print("--DONE--")
    return np.reshape(cell_interpolatedValue,(1,len(cell_interpolatedValue))) 

def polyFit(list1,list2,degree):
    """
    Finds the coefficients of the fitting polynomial given two lists of values
    and the degree of the polynomial and create a class to store the parameters
    """
    y = np.polyfit(np.array(list1),np.array(list2),degree) # actual fitting
    p = np.poly1d(y) #storing parameters in a class
    return p

def logistic_growth(xi,xi_max,sigma,deltat):
    """
    This calculates the logistic growth given the parameters and initial state
    """
    delta_xi = deltat*sigma*xi*(1-xi/xi_max)
    return delta_xi

def f_Tron(x, le, alpha):
    """
    This calculates the pdf of the water table levels. All parameters (x, alpha) should be normalized by h2 (minimum water level)
    """
    return ( alpha**(-le)/gamma(le) )*np.exp( (x-1)/alpha )*(1-x)**(le-1) 

def exponential_growth(xi,xi_max,sigma,deltat):
    """
    This calculates the exponential growth given the parameters and initial state
    """
    delta_xi = deltat*sigma*(xi_max-xi)
    return delta_xi

def calc_mean_root_biomass(root_distribution,deltaz):
    """
    This calculates the mean root biomass from the root distrbution
    """
    return integrate.trapz(root_distribution,dx=deltaz)

def dimensionless(list1,min_x,max_x):
    """
    This allows to normalize a vector by their max and min values
    """
    if (max_x-min_x) <= 0.0:
        print("Division by zero! Please check min and max values for normalization.")
        print(f'ZB is {max_x}, MIN is {min_x} and QMEAN is {list1}')
    return 1-(( list1 - min_x )/( max_x-min_x ))

def advance_in_time(var0,varmax,duration,deltat,formula,growth_rate):
    """
    This advance the solution in time depending on the growth formula for the whole duration
    """
    time = 0.0
    var = var0 # initial value
    while time <= duration:
        if formula == 'exp':
            deltavar = exponential_growth(var,varmax,growth_rate,deltat)
        elif formula == 'logistic':
            deltavar = logistic_growth(var,varmax,growth_rate,deltat)
        else:
            print(f'Only formulas are ..exp.. and ..logistic..')
            break
        # advance value of var and time
        var += deltavar
        time += deltat #NB the last timestep should be adjusted
    return var

