# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 08:07:13 2018

@author: Andrea
"""

import interpolater_function as intFunc
import numpy as np
from UnitTesterSG import *


#importing the functions from UnitTesterSG module

#extracting the digit from the file name to use as prefix/suffix in check_results
def return_digit_from_filename():
    import os
    filename=os.path.basename(__file__)
    import re
    listOfNumbers=re.findall('\d+',filename)
    #the digit is the first element in list of numbers
    extractedDigit=listOfNumbers[0]
    return extractedDigit

#BELOW ARE THE LINES INTENDED TO BE CHANGED BY THE USER	
#1) import the function whose results need to be checked
    #They are all listed above
#2) getting the prefix (or suffix) arugument for check_results. This is just for the output filenames.
suffix=return_digit_from_filename()
#3) provide the input for the function you want to test (you can also import it from a pickled object, for example)
abscissa = np.array([10.0,16])

working_data = np.array(
    [[0,-5,-2],
     [3,2,6]])

#4) get the output of the function, which is what will typically be checked.

outputNewInterpolater=intFunc.Interpolater_new(working_data,abscissa,MaxAllowedDeltaYRatio=2.0, IgnorableDeltaYThreshold = 0.0001,cleanupSuperfluousInterpolatedRows=True)

resultObj= outputNewInterpolater  #You can alternatively populate resultObj with whatever you want, such as a list.
#5) A string is also typically provided, but is an optional argument. You can provide whatever string you want.
resultStr= str(resultObj)
#6) Checking the result of the function using check_results. In this case the result is sumList1 object. 
check_results(resultObj,resultStr,prefix='',suffix=suffix)
#7) Run UnitTesterSG (or this file). UnitTesterSG will loop through all files with this format: test1.py, test2.py, etc.,
# while running an individual testing py file will only run that test. Note that UnitTesterSG expects the numbers
#to be from 1 and going up: test1, test2, test3, etc. If you go "test1, test19" then test19 will (probably) not get executed by UnitTesterSG.


