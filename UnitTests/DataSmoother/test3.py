# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 16:14:43 2017

@author: tienhung2501
"""
#THE FOLLOWING LINES ARE MANDATORY FOR THE CODE
#importing the functions from UnitTesterSG module
from UnitTesterSG import *
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
from MSRESOLVE import DataSmoother
import time
import numpy as np
#2) getting the prefix (or suffix) arugument for check_results. This is just for the output filenames.
suffix=return_digit_from_filename()
#3) provide the input for the function you want to test (you can also import it from a pickled object, for example)
data = np.array([[0., 1.5, 3.4, 5., 6.5, 7.1, 7.5],
                 [0.1, 0.2, 0.45, .45, .55, .65, 0.69]])
data = data.transpose()
abscissa = np.array([1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
headers = [34, 35]
dataSmootherChoice = 'pointrange'
dataSmootherTimeRadius = 2
dataSmootherPointRadius = 2
headersToConfineTo = []
polynomialOrder = 2

#4) get the output of the function, which is what will typically be checked. 
output = DataSmoother(data,abscissa,headers,dataSmootherChoice,dataSmootherTimeRadius,dataSmootherPointRadius,headersToConfineTo,polynomialOrder)
resultObj= output  #You can alternatively populate resultObj with whatever you want, such as a list.
#5) A string is also typically provided, but is an optional argument. You can provide whatever string you want.
resultStr= str(resultObj)
#6) Checking the result of the function using check_results. In this case the result is sumList1 object. 
check_results(resultObj,resultStr,prefix='',suffix=suffix)
#7) Run UnitTesterSG (or this file). UnitTesterSG will loop through all files with this format: test1.py, test2.py, etc.,
# while running an individual testing py file will only run that test. Note that UnitTesterSG expects the numbers
#to be from 1 and going up: test1, test2, test3, etc. If you go "test1, test19" then test19 will (probably) not get executed by UnitTesterSG.
