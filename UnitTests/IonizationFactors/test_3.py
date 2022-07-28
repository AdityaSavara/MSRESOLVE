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
MSRESOLVE.G.referencePatternsFileNamesList = ['AcetaldehydeNISTRefKnownTypes.csv'] #Overwrite with desired reference file
MSRESOLVE.G.dataToAnalyzeFileName = '2-CrotAcetExp#2Truncated.csv'
MSRESOLVE.G.ionizationDataFileName = 'ProvidedIonizationDataExample.csv'
MSRESOLVE.G.grapher = 'no'
MSRESOLVE.G.exportAtEachStep = 'no'
MSRESOLVE.G.timeRangeLimit = 'no'
#Run the main function in MSRESOLVE
MSRESOLVE.main()

#Get data from the reference file
ReferenceInfo = numpy.genfromtxt('AcetaldehydeNISTRefKnownTypes.csv',dtype=None,delimiter=',', encoding=None)
#Electron Numbers are on the third row
ElectronNumbers = ReferenceInfo[2][1:].astype(float) #convert to floats
#Ionization Types are on the fourth row
IonizationTypes = ReferenceInfo[3][1:].astype(str) #conver to strings

#In the reference data, Acetaldehyde was renamed to Ethanal and Ethanol was renamed to EtOH, this way all nine molecules will be solved via a linear fit rather than molecule matching between reference data and ionization data

#The polynomial coefficients determined in Excel (LinearFits.xlsx)
InertsCoefficients = [0.0427, 0.0871]
MainGroupCoefficients = [0.1163, -0.7131]

#Convert to poly1dObjects
InertsPolyl1dObject = numpy.poly1d(InertsCoefficients)
MainGroupPoly1dObject = numpy.poly1d(MainGroupCoefficients)

#initalize a row of zeros to hold the ionization factors
ionizationFactorsRN2 = numpy.zeros(len(ElectronNumbers))

#Find the molecules type and calculate its ionization factor
for moleculeIndex in range(len(ionizationFactorsRN2)):
    if IonizationTypes[moleculeIndex] == 'Inert':
        ionizationFactorsRN2[moleculeIndex] = numpy.polyval(InertsPolyl1dObject,ElectronNumbers[moleculeIndex])
    if IonizationTypes[moleculeIndex] == 'Main Group':
        ionizationFactorsRN2[moleculeIndex] = numpy.polyval(MainGroupPoly1dObject,ElectronNumbers[moleculeIndex])
    else:
        ionizationFactorsRN2[moleculeIndex] = ElectronNumbers[moleculeIndex]*(0.6/14) + 0.4

       

#the feature uses known ionization factors if they are available
ut.set_expected_result(ionizationFactorsRN2,expected_result_str=str(ionizationFactorsRN2),prefix=prefix,suffix=suffix)


#set output
output = MSRESOLVE.ReferenceDataList[0].relativeIonizationEfficiencies #The ionization factors list is a subobject to the MSReference object
#Places object in a tuple
resultObj = (output)

#String is provided
resultStr = str(resultObj)

#set tolerances
relativeTolerance = 1E-2
absoluteTolerance = 1E-5


#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True, relativeTolerance=relativeTolerance,absoluteTolerance=absoluteTolerance)
