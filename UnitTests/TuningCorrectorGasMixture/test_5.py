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

import UserInputGasMixture


MSRESOLVE.G = UserInputGasMixture

#now need to run MSRESOLVE which will change the coefficients...
MSRESOLVE.main()
print(MSRESOLVE.G.referenceCorrectionCoefficients)

expected_results = [2.7233E-07, 8.1371E-07, 1.83852E-06, 1.60614E-08]
ut.set_expected_result(expected_results, str(expected_results), prefix = '', suffix=suffix)

output = MSRESOLVE.resultsObjects['concentrationsScaledToCOarray'][0,1:]
resultObj = (output)
#String is provided
resultStr = str(resultObj)



#run the Unit Tester
def test_Run(allowOverwrite = False):
    #if the user wants to be able to change what the saved outputs are
    if allowOverwrite:
        #This function call is used when this test is run solo as well as by UnitTesterSG
        ut.check_results(resultObj, resultStr, prefix = '', suffix=suffix, relativeTolerance=0.05, absoluteTolerance=1.0E-5)
    #this option allows pytest to call the function
    if not allowOverwrite: 
        #this assert statement is required for the pytest module 
        assert ut.check_results(resultObj, resultStr, prefix = '', suffix=suffix, allowOverwrite = False, relativeTolerance=0.05, absoluteTolerance=1.0E-5) == True
    
if __name__ == "__main__":
   test_Run(allowOverwrite = True)