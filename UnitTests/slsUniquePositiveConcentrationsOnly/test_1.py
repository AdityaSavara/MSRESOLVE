import os
import sys
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#import the functions from UnitTesterSG
import UnitTesterSG as ut
import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)
import test_1_Input as G, imp; imp.reload(G) #This is needed because we need the __var_list__
MSRESOLVE.G = G #This is because we need to overwrite whatever user input the user has with the default user input.
    
#get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)
#prefix. Make this '' if you do not want any prefix.
prefix = ''

import UnitTesterSG as ut
import numpy as np

#Get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)

#Run the main function in MSRESOLVE
MSRESOLVE.main()

#Now we will extract the concentrations and also set the expected ones.
output = np.genfromtxt("ScaledConcentrations.csv", skip_header = 1, delimiter =",")

#Expected output:
expected_output = np.genfromtxt("ScaledConcentrationsExpected_test_1.csv", skip_header = 1, delimiter =",")
ut.set_expected_result(expected_output,str(expected_output),prefix=prefix,suffix=suffix)

resultObj = output

#String is provided
resultStr = str(resultObj)

relativeTolerance = 1.0E-2
absoluteTolerance = 1.0E-9


#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
