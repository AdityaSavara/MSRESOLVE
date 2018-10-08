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

##First Test input - First reference file
#This replaces the globals variables being pointed to in MSRESOLVE
MSRESOLVE.G.referenceFileNamesList = ['AcetaldehydeNISTRefKnownTypes.csv'] #Overwrite with desired reference file
MSRESOLVE.G.collectedFileName = '2-CrotAcetExp#2Truncated.csv'
MSRESOLVE.G.grapher = 'no'
MSRESOLVE.G.exportAtEachStep = 'no'

#Run the main function in MSRESOLVE
MSRESOLVE.main()

#Get data from the reference file
ReferenceInfo = numpy.genfromtxt('AcetaldehydeNISTRefKnownTypes.csv',dtype=None,delimiter=',',encoding=None)
#Electron Numbers are on the third row
ElectronNumbers = ReferenceInfo[2][1:].astype(float) #convert to floats
#Ionization Types are on the fourth row
IonizationTypes = ReferenceInfo[3][1:]

#In the reference data, Acetaldehyde was renamed to Ethanal and Ethanol was renamed to EtOH, this way all nine molecules will be solved via a linear fit rather than molecule matching between reference data and ionization data

#The polynomial coefficients determined in Excel (LinearFits.xlsx)
AldehydesCoefficients = numpy.array([0.1083, 0])
CarbonContainingNonmetalsCoefficients = numpy.array([0.1658, -1.5517])
AlcoholsCoefficients = numpy.array([0.0198, 2.6971])
HydrogenNonMetalIdesCoefficients = numpy.array([0.0954, 1.2171])

#Convert to poly1dObjects
AldehydesPoly1dObject = numpy.poly1d(AldehydesCoefficients)
CarbonContainingNonmetalsPoly1dObject = numpy.poly1d(CarbonContainingNonmetalsCoefficients)
AlcoholsPoly1dObject = numpy.poly1d(AlcoholsCoefficients)
HydrogenNonMetalIdesPoly1dObject = numpy.poly1d(HydrogenNonMetalIdesCoefficients)

#initalize a row of zeros to hold the ionization factors
ionizationFactorsRN2 = numpy.zeros(len(ElectronNumbers))

#Find the molecules type and calculate its ionization factor
for moleculeIndex in range(len(ionizationFactorsRN2)):
    if IonizationTypes[moleculeIndex] == 'Alcohol': #If alcohol evaluate using alcohol's poly1d object
        ionizationFactorsRN2[moleculeIndex] = numpy.polyval(AlcoholsPoly1dObject,ElectronNumbers[moleculeIndex])
    if IonizationTypes[moleculeIndex] == 'Aldehyde': #If aldehyde evaluate using aldehyde's poly1d object
        ionizationFactorsRN2[moleculeIndex] = numpy.polyval(AldehydesPoly1dObject,ElectronNumbers[moleculeIndex])
    if IonizationTypes[moleculeIndex] == 'Carbon containing nonmetals': #If carbon containing nonmetal evaluate using carbon containing nonmetal's poly1d object
        ionizationFactorsRN2[moleculeIndex] = numpy.polyval(CarbonContainingNonmetalsPoly1dObject,ElectronNumbers[moleculeIndex])
    if IonizationTypes[moleculeIndex] == 'Hydrogen non-metal-ides': #If hydrogen non-metal-ide then use hydrogen non-metal-ide's poly1d object
        ionizationFactorsRN2[moleculeIndex] = numpy.polyval(HydrogenNonMetalIdesPoly1dObject,ElectronNumbers[moleculeIndex])

       

#the feature uses known ionization factors if they are available
ut.set_expected_result(ionizationFactorsRN2,expected_result_str=str(ionizationFactorsRN2),prefix=prefix,suffix=suffix)


#set output
output = MSRESOLVE.ReferenceDataList[0].ionizationEfficienciesList #The ionization factors list is a subobject to the MSReference object
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