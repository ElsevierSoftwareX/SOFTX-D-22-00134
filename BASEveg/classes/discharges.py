# Implementation of classes
# F. Caponi, 2021

import numpy as np
from statistics import mean

class TIMESERIE:
	"""
	Class for time series of discharges. Time and Discharges values are needed to initzialize.
	"""
	def __init__(self,t,q):
		self.__time = t
		self.__value = q
		self.__incremental_time = []
		self.__data = {}
		self.__N = []
		self.__deltaT = []
		self.__type = "none"

    # some get and set functions
	def getDeltaT(self):
	 	return self.__deltaT

	def getTime(self):
	 	return self.__time
	
	def getValue(self):
	 	return self.__value

	def getMin(self):
		return min(self.__value)
	
	def getMean(self):
		return mean(self.__value)
	
	def setType(self,value):
		self.__type = value 

	def getType(self):
		return self.__type

	def getIncrementalTime(self):
		return self.__incremental_time

	def setIncrementalTime(self, value):
		self.__incremental_time =  value

	#def calc_duration_days(self):
		#return (self.__time[-1]-self.__time[0]).days

	def calc_duration_seconds(self):
		return (self.__time[-1]-self.__time[0]).total_seconds()
	
	def calc_deltat(self):
		return (self.__time[1]-self.__time[0]).total_seconds()

	def calc_duration_days(self):
		time = self.__time
		duration = 0
		for i in range(0,len(time)-1):
			#print((time[i+1]-time[i]).days)
			duration = duration + (time[i+1]-time[i]).total_seconds()
		return duration/(3600*24)

	def calc_duration_days_effective(self): # it counts the number of days excluding gap in the series
		time = self.__time
		duration = 0
		deltat0 = (time[1]-time[0]).total_seconds() # I assume constant deltat in the series (= at the first one)
		for i in range(0,len(time)-1):
			deltat = (time[i+1]-time[i]).total_seconds()
			if deltat==deltat0:
				duration = duration + deltat
		return duration/(3600*24)
