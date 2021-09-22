# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 15:20:19 2018

@author: Alex
"""

import sys
import os
baseDir = os.getcwd()
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#This test file tests the ReferencePatternTimeChooser feature

import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)
import UnitTesterSG as ut

#Get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)

#import the test input file
import test_2_input

##First Test input - First reference file
#This replaces the globals variables being pointed to in MSRESOLVE
MSRESOLVE.G = test_2_input

#Run the main function in MSRESOLVE
MSRESOLVE.main()

#Create an MSReference object for the hand created Reference Data
[provided_reference_intensities, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form]=MSRESOLVE.readReferenceFile('HandInterpolatedReferenceData.csv','xyyy')
HandInterpolatedData = [MSRESOLVE.MSReference(provided_reference_intensities, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName=referenceFileName, form=form)]
#save each global variable into the class objects 
HandInterpolatedData[0].ExportAtEachStep = MSRESOLVE.G.ExportAtEachStep
HandInterpolatedData[0].iterationSuffix = MSRESOLVE.G.iterationSuffix
#Prepare the handInterpolatedData
HandInterpolatedData[0] = MSRESOLVE.PrepareReferenceObjectsAndCorrectionValues(HandInterpolatedData[0], MSRESOLVE.ExperimentData.mass_fragment_numbers, MSRESOLVE.ExperimentData,MSRESOLVE.G.extractReferencePatternFromDataOption,MSRESOLVE.G.rpcMoleculesToChange,MSRESOLVE.G.rpcMoleculesToChangeMF,MSRESOLVE.G.rpcTimeRanges,verbose=False)

#Get reference intensities
HandInterpolatedDataIntensities = HandInterpolatedData[0].standardized_reference_patterns

#Set the expected results
#Commented out due to rounding errors: HandInterpolatedData.csv gave three small rounding errors
#expected_result = HandInterpolatedDataIntensities
#ut.set_expected_result(expected_result, str(expected_result), prefix = '', suffix=suffix)

#Get the interpolated reference pattern
#currentReferenceData is the last referencePattern used by MSRESOLVE
interpolatedIntensities = MSRESOLVE.currentReferenceData.standardized_reference_patterns

#set output
output = interpolatedIntensities
#Places object in a tuple
resultObj = (output)

#String is provided
resultStr = str(resultObj)
print(resultObj)

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