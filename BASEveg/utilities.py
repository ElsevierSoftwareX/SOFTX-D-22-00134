# This is a collection of functions useful for defining some initial parameter for the model
# which can be used independently from the rest

from classes.discharges import *
from functions.input_output import *

import os
import numpy as np
import csv
import sys
from datetime import datetime


roots = advance_in_time(0.001,2,30,1,'exp',0.1)
biomass = advance_in_time(0.001,1,30,1,'logistic',0.1)


