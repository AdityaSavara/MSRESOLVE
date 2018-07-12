# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 16:14:43 2017

@author: tienhung2501
"""
#THE FOLLOWING LINES ARE MANDATORY FOR THE CODE
import sys
sys.path.insert(1, "..\\lib")
#importing the functions from UnitTesterSG module
import UnitTesterSG as ut
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
from XYYYDataFunctionsSG import DataSmoother
import numpy as np
#2) getting the prefix (or suffix) arugument for check_results. This is just for the output filenames.
suffix=return_digit_from_filename()
#3) provide the input for the function you want to test (you can also import it from a pickled object, for example)
data = np.array([[0., 1.5, 3.4, 5., 6.5, 7.1, 7.5],
                 [0.1, 0.2, 0.45, .45, .55, .65, 0.69]])
data = data.transpose()
abscissa = np.array([1., 2., 3., 4., 5., 6., 7.])
headers = [34, 35]
dataSmootherChoice = 'timerange'
dataSmootherTimeRadius = 2
dataSmootherPointRadius = 2
headersToConfineTo = []
polynomialOrder = 2

#4) get the output of the function, which is what will typically be checked. 
output = DataSmoother(data,abscissa,headers,dataSmootherChoice,dataSmootherTimeRadius,dataSmootherPointRadius,headersToConfineTo,polynomialOrder)
resultObj= output  #You can alternatively populate resultObj with whatever you want, such as a list.
#5) A string is also typically provided, but is an optional argument. You can provide whatever string you want.
resultStr= str(resultObj)

#run the Unit Tester
def test_Run(allowOverwrite = False):
    #if the user wants to be able to change what the saved outputs are
    if allowOverwrite:
        #This function call is used when this test is run solo as well as by UnitTesterSG
        ut.check_results(resultObj, resultStr, prefix = '', suffix=suffix)
    #this option allows pytest to call the function
    if not allowOverwrite: 
        #this assert statement is required for the pytest module 
        assert ut.check_results(resultObj, resultStr, prefix = '', suffix=suffix, allowOverwrite = False) == True
    
if __name__ == "__main__":
   test_Run(allowOverwrite = True)