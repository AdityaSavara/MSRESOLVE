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

#This replaces the globals variables being pointed to in MSRESOLVE
MSRESOLVE.G.referencePatternsFileNamesList = ['AcetaldehydeNISTRefKnownFactors.csv'] #Overwrite with desired reference file
MSRESOLVE.G.dataToAnalyzeFileName = '2-CrotAcetExp#2Truncated.csv'
MSRESOLVE.G.ionizationDataFileName = 'ProvidedIonizationDataExample.csv'
MSRESOLVE.G.grapher = 'no'
MSRESOLVE.G.exportAtEachStep = 'no'
MSRESOLVE.G.timeRangeLimit = 'no'
#Run the main function in MSRESOLVE
MSRESOLVE.main()

#Get data from the reference file
ReferenceInfo = numpy.genfromtxt('AcetaldehydeNISTRefKnownFactors.csv',dtype=None,delimiter=',', encoding=None)

#ionization factors are on the fifth row
knownIonizationFactorsRN2 = ReferenceInfo[4][1:].astype(float) #convert values to a float

#the feature uses known ionization factors if they are available
ut.set_expected_result(knownIonizationFactorsRN2,expected_result_str=str(knownIonizationFactorsRN2),prefix=prefix,suffix=suffix)

#set output
output = MSRESOLVE.ReferenceDataList[0].relativeIonizationEfficiencies #The ionization factors list is a subobject to the MSReference object
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