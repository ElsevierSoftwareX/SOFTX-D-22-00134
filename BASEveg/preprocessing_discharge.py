# Script for preprocessing F. Caponi 2021

from classes.discharges import *
from functions.input_output import *

import os
import numpy as np
import csv
import sys
from datetime import datetime
import argparse
#import subprocess

def get_events_high(inputpar,q,dates):
    """
    This functions calculates the flood events in a timeseries
    """
    # Calculate number of events
    t = 0
    events = []
    while t < len(q):
        if q[t]>=inputpar['measured'][0]['Q_morpho']:
            tin = t+1
            qevent = []
            tevent = []
            while q[tin]>=inputpar['measured'][0]['Q_morpho']:
                qevent.append(q[tin])
                tevent.append(dates[tin])
                tin += 1
            t = tin-1
            if tevent:
                durationh = ((tevent[-1]-tevent[0]).total_seconds())/3600
                if durationh > inputpar['measured'][0]['T_morpho']:
                    events.append( TIMESERIE( tevent, qevent ))
        t += 1
    return events

def get_events_low(inputpar,q,dates):
    """
    This functions calculates the low flow events in a timeseries
    """
    events = []
    t = 0
    while t < len(dates):
        current_year = dates[t].year
        GSstart = datetime.strptime(inputpar['measured'][0]['growing_time'][0]+'.'+str(current_year), '%d.%m.%Y')
        GSend = datetime.strptime(inputpar['measured'][0]['growing_time'][1]+'.'+str(current_year), '%d.%m.%Y')
        if dates[t]>=GSstart and dates[t]<=GSend:
            tin = t+1
            qevent = []
            tevent = []
            while dates[tin]<=GSend:
                qevent.append(q[tin])
                tevent.append(dates[tin])
                tin += 1
            t = tin-1
            if tevent:
                events.append( TIMESERIE( tevent, qevent ))
        t += 1
    return events

def get_events_simulation(inputpar,events_high,events_low):
    """
    Returns a list of events (timeseries class) to be simulated
    """
    events = []
    order = []
    cont = 0
    for low in events_low:
        t0 = low.getTime()[0]
        t1 = low.getTime()[-1]
        idx0 = 0
        events_to_cycle = events_high[cont:]
        if cont == len(events_high):
            events_to_cycle = [events_high[-1]]
        for high in events_to_cycle:
            t = high.getTime()[0]
            if t0<t<t1:
                idx = low.getTime().index(t)
                tresult = low.getTime()[idx0:idx+1]
                qresults = low.getValue()[idx0:idx+1]
                # low flow - check duration
                if ((tresult[-1]-tresult[0]).total_seconds())/3600/24 >= inputpar['measured'][0]['T_baseflow']:
                    events.append(TIMESERIE(tresult,qresults))
                    events[-1].setType('low')
                    order.append('low')
                # flood
                events.append(high)
                events[-1].setType('high')
                order.append('high')
                #print(f"H start {events[-1].getTime()[0]} end  {events[-1].getTime()[-1]}")
                cont += 1
                try:
                    idx0 = low.getTime().index(high.getTime()[-1])
                except ValueError:
                    break
            elif t<t0:
                # flood
                events.append( high )
                order.append('high')
                events[-1].setType('high')
                #print(f"H start {events[-1].getTime()[0]} end  {events[-1].getTime()[-1]}")
                cont += 1
            elif t>t1:
                tresult = low.getTime()[idx0:]
                qresults = low.getValue()[idx0:]
                # low flow
                if ((tresult[-1]-tresult[0]).total_seconds())/(3600*24) >= inputpar['measured'][0]['T_baseflow']:
                    events.append(TIMESERIE(tresult,qresults))
                    #print(f"L start {events[-1].getTime()[0]} end  {events[-1].getTime()[-1]}")
                    order.append('low')
                    events[-1].setType('low')
                    break
    return events, order

def merge_events(inputpar,events,events_order):
    """
    This merges consecutive events that are equal (two floods or tow low flows) and return new list
    """ 
    newevents = []
    j = 0
    while j<len(events)-1:
        newt = events[j].getTime()
        newq = events[j].getValue()
        while events_order[j] == events_order[j+1]:
            newt.extend(events[j+1].getTime())
            newq.extend(events[j+1].getValue())
            j += 1
            if j>=len(events_order)-1:
                break
        # check first event
        if j==0 and events_order[j]=="high":
            print("Skip the first flood event. Simulate this before starting a simulation with vegetation.")
            break
        newevents.append(TIMESERIE(newt,newq))
        if events_order[j]=="low":
            newevents[-1].setType('low')
        else:
            newevents[-1].setType('high')
        j += 1
    # skip last event if it is a "low" event
    if newevents[-1].getType()=="low":
        newevents.pop()
    return newevents

def create_step_wise_discharge(inputpar, events):
    """
    This creates inflow boundary conditions for step-wise simulations (HQ relation)
    """
    stepwise_events = []
    for ee in events:
        if ee.getType() == 'low':
            # get some parameters
            qmin = min( ee.getValue() )
            qmax = max( ee.getValue() )
            tsteps = inputpar['measured'][0]['T_steps']
            nsteps = inputpar['measured'][0]['n_steps']
            # fill a list with values
            q = []
            t = []
            deltat2 = 1000
            deltaq = (qmax-qmin)/(nsteps-1)
            qstart = qmin
            tstart = 0
            for i in range(0,nsteps):
                # 1
                q.append(qstart)
                t.append(tstart)
                tstart = tstart+tsteps
                if i==0:
                    tsteps = tsteps - deltat2
                # 2
                q.append(qstart)
                t.append(tstart)
                # 3
                tstart = tstart + deltat2
                qstart = qstart + deltaq
                if i==nsteps-1:
                    qstart = qmax
            stepwise_events.append(TIMESERIE(t,q))
    return stepwise_events

def create_const_discharge(inputpar):
    """
    This create a TIMESERIE event with constant discharge
    """
    q = inputpar["measured"][0]["Q_morpho"]
    t = inputpar["measured"][0]["T_steps"]
    return TIMESERIE([0,t],[q,q])

def modify_incremental_times(events,events_steps,event_const):
    """
    This modifies the times in the given events to get incremental times
    """
    # length of steps and const discharge are always the same
    duration_steps = events_steps[0].getTime()[-1]
    duration_const = event_const.getTime()[-1]
    #print(duration_const)
    #print(duration_steps)

    nevents = len(events_steps)
    jj = 1 # I assume here that it is always low->high->low->high
    q_steps = []
    q_const = []
    q_flood = []
    increase = 0
    for ii in range(0,nevents):
        
        duration_flood = events[jj].calc_duration_seconds()
        duration_cycle = duration_steps + duration_const + duration_flood
        #print(duration_flood)

        q_steps.append(TIMESERIE([x+increase for x in events_steps[ii].getTime()],events_steps[ii].getValue()) )
        q_const.append(TIMESERIE([y+increase+duration_steps for y in event_const.getTime()],event_const.getValue()) )

        t = []
        t.append(q_const[-1].getTime()[-1])
        deltat = events[jj].calc_deltat()
        for p in range(0,len(events[jj].getTime())-1):
            t.append(t[p]+deltat)
        q_flood.append(TIMESERIE(t,events[jj].getValue()) )

        #print(q_steps[-1].getTime()[0])
        #print(q_flood[-1].getTime()[-1])
        #print('-')
        #print('încrease')
        #print(increase)
        #print('-')

        increase = q_flood[-1].getTime()[-1]#increase + duration_cycle
        jj += 2

    return q_steps, q_const, q_flood


def create_discharge_time_series(inputpath):
    """
    This creates all the inflow input for BM3 when using the vegetation and additional ouputs
    """
    # read json file with input
    inputpar = read_json(os.path.join(inputpath,'input_parameters_veg_pre.json'))

    # read txt file with time series of discharges
    formatTime = '%d.%m.%Y %H:%M:%S'
    formattoMatch = '\d{2}.\d{2}.\d{4} \d{2}:\d{2}:\d{2}'
    build_qt = False
    try:
        q = read_float_values(os.path.join(inputpath,inputpar['measured'][0]['file']),2)
        dates = read_datetime(os.path.join(inputpath,inputpar['measured'][0]['file']),formattoMatch,formatTime)
    except FileNotFoundError:
        # if none -> build a series from the input parameters (for artificial discharge series)
        #print('File with discharges not found. Building the discharge time series artificially from input parameters.')
        build_qt = True

    if build_qt:
        print('Building inflow time series from input parameters (ARTIFICIAL).')
        # get number of events (event = low flow + flood)
        n_events = len(inputpar['artificial'][0]['flood'][0]['Q_max'])
        print(f'...simulating {n_events} events.')
        
        # fill a list with TIMESERIES objects
        events = []
        tstart = 0
        for i in range(0,n_events):
            tend = tstart + inputpar['artificial'][0]['low_flow'][0]['T_steady']
            q = inputpar['artificial'][0]['low_flow'][0]['Q_mean'][i] # constant Q  
            events.append( TIMESERIE( [tstart, tend], [q, q] ))
            #
            tstart = tend
            #
            tmid =  tstart + inputpar['artificial'][0]['flood'][0]['T_rise'][i]
            tend = tmid + inputpar['artificial'][0]['flood'][0]['T_fall'][i]
            qbase = inputpar['artificial'][0]['flood'][0]['Q_base'][i]
            qmax = inputpar['artificial'][0]['flood'][0]['Q_max'][i]
            events.append( TIMESERIE( [tstart, tmid, tend], [qbase, qmax, qbase] ))
            #
            tstart = tend

    else:
        print('Building inflow time series from input parameters (MEASURED).')

        # Identify events depending on threshold defined
        events_high = get_events_high(inputpar,q,dates)
        events_low = get_events_low(inputpar,q,dates)

        # create a new list of events in temporal order 
        events_mod, events_order = get_events_simulation(inputpar,events_high,events_low)
        
        # merge consecutive low flow (or high flow) periods
        events_mod2 = merge_events(inputpar,events_mod,events_order)

        # build step-wise discharge input
        events_steps = create_step_wise_discharge(inputpar, events_mod2)

        # build constant discharge input (just once)
        event_const_q = create_const_discharge(inputpar)

        # modify times of events to get incremental...put in a dictonary the results
        q_steps, q_const, q_flood = modify_incremental_times(events_mod2,events_steps,event_const_q)

    print(f"The number of events is: {len(events_mod2)}")
    cont = 1
    for e in range(0,len(events_mod2)):
        datetoscreen1 = events_mod2[e].getTime()[0].strftime("%d.%m.%Y %H:%M:%S")
        datetoscreen2 = events_mod2[e].getTime()[-1].strftime("%d.%m.%Y %H:%M:%S")
        print(f"E{cont}: date {datetoscreen1} - {datetoscreen2}. Qpeak {max(events_mod2[e].getValue())}")
        cont += 1

    # Write to txt files
    if build_qt:
        cont = 1    
        for event in events:
            write_txt_time_series(event,inputpath,f"inflow_discharge{cont}")
            cont += 1
    else:
        # with real time series of discharge
        cont = 1    
        # create folder with results
        if not os.path.exists(os.path.join(inputpath,'date_discharge')):
            os.mkdir(os.path.join(inputpath,'date_discharge'))
        for event in events_mod2:
            write_txt_date_series(event,os.path.join(inputpath,'date_discharge'),f"inflow_discharge_date{cont}")
            cont += 1
        
        # step wise discharge
        cont = 1   
        if not os.path.exists(os.path.join(inputpath,'inflow_steps')): 
            os.mkdir(os.path.join(inputpath,'inflow_steps'))
        for event in q_steps:
            write_txt_timenum_series(event,os.path.join(inputpath,'inflow_steps'),f"inflow_discharge_steps{cont}")
            cont += 1

        # constant discharge
        cont = 1
        if not os.path.exists(os.path.join(inputpath,'inflow_const')): 
            os.mkdir(os.path.join(inputpath,'inflow_const'))
        for event in q_const:
            write_txt_timenum_series(event,os.path.join(inputpath,'inflow_const'),f"inflow_discharge_const{cont}")
            cont += 1

        # flood discharge
        cont = 1   
        if not os.path.exists(os.path.join(inputpath,'inflow_flood')): 
            os.mkdir(os.path.join(inputpath,'inflow_flood'))
        for event in q_flood:
            write_txt_timenum_series(event,os.path.join(inputpath,'inflow_flood'),f"inflow_discharge_flood{cont}")
            cont += 1

        cont = 1    
        # create statistics
        if not os.path.exists(os.path.join(inputpath,'statistics')):
            os.mkdir(os.path.join(inputpath,'statistics'))
        write_txt_statistics(events_mod2,os.path.join(inputpath,'statistics'))


if __name__ == '__main__':
    # parser
    parser = argparse.ArgumentParser(description='This preprocessing script reads and analyses discharge time series and produces input files for BASEveg.')
    parser.add_argument('-f', type=str, help='This is the path to the folder where input files are located. Input files are: discharge time series and input parameters (JSON).')
    args=parser.parse_args()
    #
    if len(sys.argv) <= 1:
        print('Please insert the input path with flag -f or see the help page.')
        exit(1)
    #path = args.f   
    create_discharge_time_series(args.f)
