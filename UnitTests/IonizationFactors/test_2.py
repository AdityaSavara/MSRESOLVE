import os
import sys
baseDir = os.getcwd()
sys.path.insert(1, os.path.join(baseDir, os.pardir))
sys.path.insert(1, os.path.join(baseDir, os.pardir, os.pardir))
#import the functions from UnitTesterSG
import UnitTesterSG as ut
import DefaultUserInput as G #This is needed because we need the __var_list__
    
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
import UserInputIonizationUnitTests

##First Test input - First reference file
#This replaces the globals variables being pointed to in MSRESOLVE
MSRESOLVE.G = UserInputIonizationUnitTests
MSRESOLVE.G.referenceFileNamesList = ['AcetaldehydeNISTRefMatchingMolecule.csv'] #Overwrite with desired reference file

#Run the main function in MSRESOLVE
MSRESOLVE.main()

#Get data from the reference file
ReferenceInfo = numpy.genfromtxt('AcetaldehydeNISTRefMatchingMolecule.csv',dtype=None,delimiter=',',encoding=None)

#ionization factors are on the fifth row
ionizationFactorsRN2 = ReferenceInfo[4][1:]

#The only two molecules in the reference data and _ProvidedIonizationData.csv that match is Acetaldehyde and Ethanol
#From the data: Acetaldehyde has an RS_Value of 2.6
#From the data: Ethanol has an average RS_Value of 3.25
#Acetaldehyde is in the first column of the reference data and ethanol is in the sixth column
ionizationFactorsRN2[0] = '2.6'
ionizationFactorsRN2[5] = '3.25'


#convert to float
ionizationFactorsRN2 = ionizationFactorsRN2.astype(float)


#the feature uses known ionization factors if they are available
ut.set_expected_result(ionizationFactorsRN2,expected_result_str=str(ionizationFactorsRN2),prefix=prefix,suffix=suffix)

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