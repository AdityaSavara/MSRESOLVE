import os
import sys
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#import the functions from UnitTesterSG
import UnitTesterSG as ut
import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)
import test_5_Input as G, importlib; importlib.reload(G) #This is needed because we need the __var_list__
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
output1 = np.genfromtxt("ScaledConcentrations.csv", skip_header = 1, delimiter =",")
output2 = np.genfromtxt("ScaledConcentrations_relative_uncertainties.csv", skip_header = 1, delimiter =",")

#Expected output:
try:
    import uncertainties
    from uncertainties import unumpy
    expected_output1 = np.genfromtxt("ScaledConcentrationsExpected_test_5.csv", skip_header = 1, delimiter =",", dtype = 'f8')
    expected_output2 = np.genfromtxt("ScaledConcentrations_relative_uncertainties_test_5.csv", skip_header = 1, delimiter =",", dtype = 'f8')
    print("THE UNCERTAINTIES MODULE ***IS*** PRESENT, RUNNING THE UNIT TEST ACCORDINGLY.")
except:
    print("THE UNCERTAINTIES MODULE ***IS NOT*** PRESENT, RUNNING THE UNIT TEST ACCORDINGLY.")
    expected_output1 = np.genfromtxt("ScaledConcentrationsExpected_test_5.csv", skip_header = 1, delimiter =",")
    expected_output2 = np.genfromtxt("ScaledConcentrations_relative_uncertainties_test_5_noUncertaintiesModule.csv", skip_header = 1, delimiter =",")

ut.set_expected_result((expected_output1,expected_output2) ,str((expected_output1,expected_output2)),prefix=prefix,suffix=suffix)

resultObj = (output1, output2)

#String is provided
resultStr = str(resultObj)

relativeTolerance = 1.0E-1
absoluteTolerance = 1.0E-9


#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
