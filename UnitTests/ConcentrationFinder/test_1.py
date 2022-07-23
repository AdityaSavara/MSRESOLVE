import os
import sys
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#import the functions from UnitTesterSG
import UnitTesterSG as ut
import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)
import DefaultUserInput as G, importlib; importlib.reload(G) #This is needed because we need the __var_list__
MSRESOLVE.G = G #This is because we need to overwrite whatever user input the user has with the default user input.
    
#get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)
#prefix. Make this '' if you do not want any prefix.
prefix = ''

import UnitTesterSG as ut
import numpy as np

#Get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)

#Replace the global variables being pointed to in MSRESOLVE
#This reference file has two molecules: Acetaldehyde and Acetaldehyde_easy_to_ionize
#They are identical with the exceptions that Acetaldehyde has an ionization factor of 1 and signal peaks at m29 and m29.2 are 9999 and 0, respectively, and Acetaldehyde_easy_to_ionize has an ionization factor of 2 and signal peaks at m29 and m29.2 are 0 and 9999, respectively.
MSRESOLVE.G.referencePatternsFileNamesList = ['AcetaldehydeNISTRefMixed2_test_1.csv']
#Use the truncated reference data that contains only m29 and m29.2 with signals of 1 for each time for both mass fragments
MSRESOLVE.G.dataToAnalyzeFileName = '2-CrotAcetExp#2Truncated.csv'
MSRESOLVE.G.concentrationFinder = 'yes' #Turn on concentrationFinder
MSRESOLVE.G.TSC_List_Type = 'SeparateMolecularFactors' #Use factors for separate molecules rather than numerous reference files
MSRESOLVE.G.moleculesTSC_List = ['Acetaldehyde'] #Using the 'known' concentration values for Acetaldehyde with SeparateMoleculesFactors, concentrationFinder will calculate a conversion factor based on this info and apply to it the scaled concentration of itself and Acetaldehyde_easy_to_ionize
MSRESOLVE.G.massNumberTSC_List = [29]
MSRESOLVE.G.moleculeSignalTSC_List = [1.66945]
MSRESOLVE.G.moleculeConcentrationTSC_List = [0.05]
MSRESOLVE.G.grapher = 'no'
MSRESOLVE.G.exportAtEachStep = 'no'
MSRESOLVE.G.timeRangeLimit = 'no'

#Run the main function in MSRESOLVE
MSRESOLVE.main()

ResolvedConcentrationsData = MSRESOLVE.resultsObjects['concentrationsarray'] #get the concentrations array from the global resultsObjects dictionary
ut.set_expected_result(2.0,str(2.0),prefix=prefix,suffix=suffix)

#set output
output = ResolvedConcentrationsData[0][1]/ResolvedConcentrationsData[0][2] #find the ratio of the second column to the third column.  The time value is the first column.  Use the first value in the column since each value in a particular column is the same
#Places object in a tuple
resultObj = output

#String is provided
resultStr = str(resultObj)

relativeTolerance = 1.0E-2
absoluteTolerance = 1.0E-2


#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
