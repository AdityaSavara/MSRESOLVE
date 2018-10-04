import os
import sys
baseDir = os.getcwd()
sys.path.insert(1, os.path.join(baseDir, os.pardir))
sys.path.insert(1, os.path.join(baseDir, os.pardir, os.pardir))
#import the functions from UnitTesterSG
import UnitTesterSG as ut
import DefaultUserInput as G #This is needed because we need the __var_list__
MSRESOLVE_var_list = G.__var_list__ #need to store this to reassign in the new namespace.
    
#get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)
#prefix. Make this '' if you do not want any prefix.
prefix = ''

import MSRESOLVE
import UnitTesterSG as ut
import numpy

#Get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)

#import the test input file
import test_1_input

##First Test input - First reference file
#This replaces the globals variables being pointed to in MSRESOLVE
MSRESOLVE.G = test_1_input

#Run the main function in MSRESOLVE
MSRESOLVE.main()

#Get data from the reference file
ReferenceInfo = numpy.genfromtxt('AcetaldehydeNISTRefKnownFactors.csv',dtype=None,delimiter=',',encoding=None)

#ionization factors are on the fifth row
knownIonizationFactorsRN2 = ReferenceInfo[4][1:].astype(float) #convert values to a float

#the feature uses known ionization factors if they are available
ut.set_expected_result(knownIonizationFactorsRN2,expected_result_str=str(knownIonizationFactorsRN2),prefix=prefix,suffix=suffix)

#set output
output = MSRESOLVE.ReferenceDataList[0].ionizationEfficienciesList #The ionization factors list is a subobject to the MSReference object
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