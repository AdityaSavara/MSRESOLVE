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

#We can use the default user input, which is already in MSRESOLVE.G.
#Need to change the "direct variable" version rather than the dictionary version.
MSRESOLVE.G.referencePatternsFileNamesList = ['AcetaldehydeMeasured.tsv']
MSRESOLVE.G.tuningCorrection = 'yes'
MSRESOLVE.G.createMixedTuningPattern = True
MSRESOLVE.G.referenceFileExistingTuningAndForm = ['ReferenceLiterature.tsv','xyyy']
MSRESOLVE.G.referenceFileDesiredTuningAndForm =[]

#apparently need to have dataAnalysis on to use this feature (should not need to, but as of Sep 2019, do need to.
MSRESOLVE.G.dataAnalysis ='yes'

#Turn off certain other settings to make the test faster.
MSRESOLVE.G.grapher ='no'
MSRESOLVE.G.dataSimulation ='no'
MSRESOLVE.G.timeRangeLimit = 'yes'
MSRESOLVE.G.timeRangeStart = 176  #start time (-int)
MSRESOLVE.G.timeRangeFinish = 178  #start time (-int)
MSRESOLVE.G.dataSmootherYorN = 'no'
MSRESOLVE.G.calculateUncertaintiesInConcentrations = False
#will use applyReferenceMassFragmentsThresholds (single value means it will be applied to all molecules)
MSRESOLVE.G.applyReferenceMassFragmentsThresholds= 'no'
MSRESOLVE.G.referenceMassFragmentFilterThreshold = [0.0]
MSRESOLVE.G.referenceSignificantFragmentThresholds = [0.0]
MSRESOLVE.G.solverChoice='inverse'
#Now to get started with the test itself...

#these are the tuning coefficients before doing anything (should be 0,0,1)
print(MSRESOLVE.G.referenceCorrectionCoefficients)
#now need to run MSRESOLVE which will change the coefficients...
MSRESOLVE.main()
print(MSRESOLVE.G.referenceCorrectionCoefficients)

#expected_results = (MSRESOLVE.G.referenceCorrectionCoefficients,MSRESOLVE.currentReferenceData.standardized_reference_patterns)
#ut.set_expected_result(expected_results, str(expected_results), prefix = '', suffix=suffix)

output =MSRESOLVE.currentReferenceData.standardized_reference_patterns
resultObj = (output)
#String is provided
resultStr = str(resultObj)



#run the Unit Tester
def test_Run(allowOverwrite = False):
    #if the user wants to be able to change what the saved outputs are
    if allowOverwrite:
        #This function call is used when this test is run solo as well as by UnitTesterSG
        ut.check_results(resultObj, resultStr, prefix = '', suffix=suffix, relativeTolerance=1.0e-5, absoluteTolerance=1.0E-5)
    #this option allows pytest to call the function
    if not allowOverwrite: 
        #this assert statement is required for the pytest module 
        assert ut.check_results(resultObj, resultStr, prefix = '', suffix=suffix, allowOverwrite = False, relativeTolerance=1.0e-5, absoluteTolerance=1.0E-5) == True
    
if __name__ == "__main__":
   test_Run(allowOverwrite = True)