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
import numpy

#Get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)

##First Test input - First reference file
#This replaces the globals variables being pointed to in MSRESOLVE
MSRESOLVE.G.referencePatternsFileNamesList = ['AcetaldehydeNISTRefMatchingMolecule.csv'] #Overwrite with desired reference file
MSRESOLVE.G.dataToAnalyzeFileName = '2-CrotAcetExp#2Truncated.csv'
MSRESOLVE.G.ionizationDataFileName = 'ProvidedIonizationDataExample.csv'
MSRESOLVE.G.grapher = 'no'
MSRESOLVE.G.exportAtEachStep = 'no'
MSRESOLVE.G.timeRangeLimit = 'no'
#Run the main function in MSRESOLVE
MSRESOLVE.main()

#Get data from the reference file
ReferenceInfo = numpy.genfromtxt('AcetaldehydeNISTRefMatchingMolecule.csv',dtype=None,delimiter=',', encoding=None)

#Get a list the same length as the number of molecules
ionizationFactorsRN2 = numpy.zeros(len(ReferenceInfo[0][1:]))

#The only two molecules in the reference data and ProvidedIonizationDataExample.csv that match is Carbon Monoxide and Carbon Dioxide
#From the data: Carbon Monoxide has an RS_Value of 1.05
#From the data: Carbon Dioxide has an average RS_Value of 1.4
#Acetaldehyde is in the first column of the reference data and ethanol is in the sixth column
ionizationFactorsRN2[2] = 1.05 
ionizationFactorsRN2[3] = 1.4
ionizationFactorsOutput = ionizationFactorsRN2[2:4] #get an array of just the two ionization factors


#set the expected results to be the output array
ut.set_expected_result(ionizationFactorsOutput[2:4],expected_result_str=str(ionizationFactorsRN2[2:4]),prefix=prefix,suffix=suffix)

#set output
output = MSRESOLVE.ReferenceDataList[0].relativeIonizationEfficiencies[2:4] #The ionization factors list is a subobject to the MSReference object, indices 2 and 3 refer to CO and CO2
#Places object in a tuple
resultObj = (output)

#String is provided
resultStr = str(resultObj)

#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False)
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True)