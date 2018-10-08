"""
Created on Wed Jun 13 08:07:13 2018

@author: Andrea
"""
import sys
import os
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#import the functions from UnitTesterSG
import UnitTesterSG as ut
#import your function
import XYYYDataFunctionsSG as dataFunctions
import numpy as np

#get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)

#input for the unit that will be tested
abscissa = np.array([10.0,16.0])

working_data = np.array(
    [[0,-5,-2],
     [3,2,6]])
    
concentrationBounds=np.array([[10,0.02,0.03,0.001,0.04,1,0.001],
                             [16,0.01,0.02,0.001,1,1.2,0.001]])

#These two lines can hardcode the expected results. They are not required. 
#expected_results = 6
#ut.set_expected_result(expected_results, str(expected_results), prefix = '', suffix=suffix)

#outputs with the function being tested using the input
outputmarginalChangeRestrictor=dataFunctions.marginalChangeRestrictor(working_data,abscissa,MaxAllowedDeltaYRatio=2.0, IgnorableDeltaYThreshold = 0.0001,extraOutput=True)
outputInterpolateAccompanyingArrays=dataFunctions.interpolateAccompanyingArrays(outputmarginalChangeRestrictor[1], concentrationBounds)

#places the object in a tuple
resultObj = (outputInterpolateAccompanyingArrays)

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
