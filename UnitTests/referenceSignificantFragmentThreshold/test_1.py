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

#import the variables from the test input file
import test_1_UserInput as G

##First Test input - First reference file
MSRESOLVE.G = G
#Run MSRESOLVE main
MSRESOLVE.main()
#Get ExportedSLSUniqueMassFragmentsUsed
ExportedSLSUniqueMassFragmentsUsed = pandas.read_csv('ExportedSLSUniqueMassFragmentsUsed.csv')

#In this example, we know that mass 14 is going to have a value of 1, and all other values will be 0.
valuesAt14 = (np.array(ExportedSLSUniqueMassFragmentsUsed.loc[:,' 14.0'])) #.loc allows finding by name of column.

ExportedSLSUniqueMassFragmentsUsedMinus14 = ExportedSLSUniqueMassFragmentsUsed
del ExportedSLSUniqueMassFragmentsUsedMinus14[' 14.0']

valuesAtNot14 = np.array(ExportedSLSUniqueMassFragmentsUsedMinus14.iloc[:,2:]) # .iloc allows finding by indices.

meanAt14 = valuesAt14.mean()
meanAtNot14 = valuesAtNot14.mean()

expected_results = (1.0,0.0) #we expect 14 to be used, and others not.
ut.set_expected_result(expected_results, str(expected_results), prefix = '', suffix=suffix)

output = (meanAt14, meanAtNot14)
resultObj = (output)
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
