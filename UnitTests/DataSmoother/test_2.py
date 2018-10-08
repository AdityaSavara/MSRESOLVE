# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 16:14:43 2017

@author: tienhung2501
"""
#THE FOLLOWING LINES ARE MANDATORY FOR THE CODE
import sys
import os
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#importing the functions from UnitTesterSG module
import UnitTesterSG as ut

#BELOW ARE THE LINES INTENDED TO BE CHANGED BY THE USER	
#1) import the function whose results need to be checked
from XYYYDataFunctionsSG import DataSmoother
import numpy as np
#2) getting the prefix (or suffix) arugument for check_results. This is just for the output filenames.
suffix=ut.returnDigitFromFilename(__file__)
#3) provide the input for the function you want to test (you can also import it from a pickled object, for example)
data = np.array([[0., 1.5, 3.4, 5., 6.5, 7.1, 7.5],
                 [0.1, 0.2, 0.45, .45, .55, .65, 0.69]])
data = data.transpose()
abscissa = np.array([1., 2., 3., 4., 5., 6., 7.])
headers = [34, 35]
dataSmootherChoice = 'timerange'
dataSmootherTimeRadius = 2
dataSmootherPointRadius = 2
headersToConfineTo = [34]
polynomialOrder = 2

#4) get the output of the function, which is what will typically be checked. 
output = DataSmoother(data,abscissa,headers,dataSmootherChoice,dataSmootherTimeRadius,dataSmootherPointRadius,headersToConfineTo,polynomialOrder)
resultObj= output  #You can alternatively populate resultObj with whatever you want, such as a list.
#5) A string is also typically provided, but is an optional argument. You can provide whatever string you want.
resultStr= str(resultObj)


relativeTolerance=1.0E-5
absoluteTolerance=1.0E-8

#run the Unit Tester
def test_Run(allowOverwrite = False):
    #if the user wants to be able to change what the saved outputs are
    if allowOverwrite:
        #This function call is used when this test is run solo as well as by UnitTesterSG
        ut.check_results(resultObj, resultStr, prefix = '', suffix=suffix, allowOverwrite = True, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
    #this option allows pytest to call the function
    if not allowOverwrite: 
        #this assert statement is required for the pytest module 
        assert ut.check_results(resultObj, resultStr, prefix = '', suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance) == True
    
    
if __name__ == "__main__":
   test_Run(allowOverwrite = True)
