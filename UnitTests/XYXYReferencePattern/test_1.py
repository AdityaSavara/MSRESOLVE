# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 15:20:19 2018

@author: Alex
"""

import sys
import os
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#This test file tests the ReferencePatternTimeChooser feature

import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)
import UnitTesterSG as ut
import pandas
import numpy as np

#Get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)

#import the test input files
import test_1_inputXYYY
import test_1_inputXYXY

##First Test input - First reference file
#This replaces the globals variables being pointed to in MSRESOLVE
MSRESOLVE.G = test_1_inputXYYY
#Run MSRESOLVE main
MSRESOLVE.main()
#Scaled concentrations are stored in ScaledConcentrations.csv
scaledConcentrationsA = pandas.read_csv('ScaledConcentrations.csv',header = 0)
scaledConcentrationsXYYYArray = np.array(scaledConcentrationsA)

##Second test input - Second reference file
#Replace global variables being pointed to in MSRESOLVE
MSRESOLVE.G = test_1_inputXYXY
#Run MSRESOLVE main
MSRESOLVE.main()
#Get scaledConcentrations
scaledConcentrationsB = pandas.read_csv('ScaledConcentrations.csv',header = 0)
scaledConcentrationsXYXYArray = np.array(scaledConcentrationsB)

#Set the expected results
expected_results = scaledConcentrationsXYYYArray
ut.set_expected_result(expected_results, str(expected_results), prefix = '', suffix=suffix)

#Set the output
output = scaledConcentrationsXYXYArray
resultObj = output
#String is provided
resultStr = str(resultObj)



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