import os
import sys
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#import the functions from UnitTesterSG
import UnitTesterSG as ut
    
#get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)
#prefix. Make this '' if you do not want any prefix.
prefix = ''

import numpy as np
import copy
import shutil
import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)

import test_2_Input
MSRESOLVE.G = test_2_Input
#this input file uses the applyRawSignalThresholds,which excludes molecules if they have significant mass fragments
# that are not present above a particular signal threshold.
# in this case, we set the signal threshold to be 0.02, and the significance to be at 1% for standardized reference patterns.
# the mass 70 signal is always below 0.02 and a significant for both crotyl alcohol and 2-butenal (crotonaldehyde)
# consequently, the solved concentrations for both of these molecules is zero at all times.
# Although we do not check it, it turns out that CO2 ends up as zero as well because after ethanol's contribution is subtracted, 
#   m45 goes below threshold, and CO2 has m45 above 1% in standardized signals. CO must have a similar reason for becoming set to zero.
MSRESOLVE.main()


ExportedSLSUniqueMassesUsedInSolvingMolecules = np.genfromtxt("ExportedSLSUniqueMassesUsedInSolvingMolecules.csv", skip_header = 0, delimiter=",", unpack=False)
# I know from previous experience that for this set of files Scaled Concentrations has this order:
# Time,
#Acetaldehyde Concentration Relative to CO,
#(E) 2-Butenal (Crotonaldehyde Concentration Relative to CO,
#CO Concentration Relative to CO,
#CO2 Concentration Relative to CO,
#Ethylene (Ethene) Concentration Relative to CO,
#Ethanol Concentration Relative to CO,
#Crotyl Alcohol Concentration Relative to CO,
#H2 Concentration Relative to CO,
#H2O Concentration Relative to CO
#H2 Concentration Relative to CO

#now, we will check that the values are expected in the columns.

resultObj = list(ExportedSLSUniqueMassesUsedInSolvingMolecules[0][2:]) #We discard the first 2 columns because those are time relatd.

expected_result = [57,70,56,54,44,58,59,45]



ut.set_expected_result(expected_result,expected_result_str=str(expected_result), prefix=prefix,suffix=suffix)

#String must be provided provided. Make it '' if you do not want to use a result string.
resultStr = str(resultObj)


#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False)
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True)