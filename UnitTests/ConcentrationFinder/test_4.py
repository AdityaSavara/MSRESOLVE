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
MSRESOLVE.G.referencePatternsFileNamesList = ['AcetaldehydeNISTRefMixed2_test_2.csv','AcetaldehydeNISTRefMixed2_test_2.csv'] #List the reference file twice
MSRESOLVE.G.referencePatternTimeRanges = [[1,1],[8,8]] #Give time ranges for each reference file (Every point that is not at 1 or 8 will need to be interpolated)
#Use the truncated reference data that contains only 8 points with a signal of 1 for m29
MSRESOLVE.G.dataToAnalyzeFileName = '2-CrotAcetExp#2Truncated2.csv'
MSRESOLVE.G.concentrationFinder = 'yes' #Turn on concentrationFinder
MSRESOLVE.G.TSC_List_Type = 'MultipleReferencePatterns' #Use factors for numerous reference files
MSRESOLVE.G.moleculesTSC_List = ['Acetaldehyde','Acetaldehyde'] #Use the same molecule
MSRESOLVE.G.massNumberTSC_List = [29,29] #Same mass fragments
MSRESOLVE.G.moleculeSignalTSC_List = [1.66945,1.66945] #With the same signal
MSRESOLVE.G.moleculeConcentrationTSC_List = [0.05,0.1] #Make the known concentration for the second reference pattern to be 2x what the known concentration is for the first reference file
MSRESOLVE.G.grapher = 'no'
MSRESOLVE.G.exportAtEachStep = 'no'
MSRESOLVE.G.timeRangeLimit = 'no'

#Run the main function in MSRESOLVE
MSRESOLVE.main()

ResolvedConcentrationsData = MSRESOLVE.resultsObjects['concentrationsarray'] #get the concentrations array from the global resultsObjects dictionary
TimesList = ResolvedConcentrationsData[:,0] #Get the times
ResolvedConcentrations = ResolvedConcentrationsData[:,1:] #Get the resolved concentrations

ExpectedResolvedConcentrations = np.zeros(len(ResolvedConcentrations)) #initialize as an array of zeros the length of resolved concentrations

ExpectedResolvedConcentrations[0] = ResolvedConcentrations[0] #The first concentration should match the first concentration in ResolvedConcentrations


for concentrationIndex in range(1,len(ExpectedResolvedConcentrations)-1): #Loop from index 2 to the next to last index
    #Interpolate the concentration between the first concentration (at the first time) and the last concentration (at the last time) at the current concentration's time index
    ExpectedResolvedConcentrations[concentrationIndex] = MSRESOLVE.DataFunctions.analyticalLinearInterpolator(ResolvedConcentrations[0],ResolvedConcentrations[-1],TimesList[concentrationIndex],TimesList[0],TimesList[-1])

ExpectedResolvedConcentrations[-1] = ResolvedConcentrations[-1] #The last concentration should match the last concentration in ResolvedConcentrations

ut.set_expected_result(ExpectedResolvedConcentrations,str(ExpectedResolvedConcentrations),prefix=prefix,suffix=suffix) #set the expected result to be the ExpectedResolvedConcentrations list

#set output
output = ResolvedConcentrations.flatten() #the output is our resolved concentrations. Before Numpy 1.16, this flatten was not necessary. But for some reason, Numpy 1.16 changed array nesting somewhere in the MSRESOLVE process, so flattening is necessary here in order to get them to match (or, the expected concentrations would need to become more nested, which was less convenient to do). Flattening here also keeps the test backwards compatible w/ numpy 1.14.

#Places object in a tuple
resultObj = output

#String is provided
resultStr = str(resultObj)

relativeTolerance = 1.0E-5
absoluteTolerance = 1.0E-8


#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)