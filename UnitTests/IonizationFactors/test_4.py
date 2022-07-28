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
MSRESOLVE.G.referencePatternsFileNamesList = ['AcetaldehydeNISTRefDefault.csv'] #Overwrite with desired reference file
MSRESOLVE.G.dataToAnalyzeFileName = '2-CrotAcetExp#2Truncated.csv'
MSRESOLVE.G.ionizationDataFileName = 'ProvidedIonizationDataExample.csv'
MSRESOLVE.G.grapher = 'no'
MSRESOLVE.G.exportAtEachStep = 'no'
MSRESOLVE.G.timeRangeLimit = 'no'
#Run the main function in MSRESOLVE
MSRESOLVE.main()

#Get data from the reference file
ReferenceInfo = numpy.genfromtxt('AcetaldehydeNISTRefDefault.csv',dtype=None,delimiter=',', encoding=None)
#Electron Numbers are on the third row
ElectronNumbers = ReferenceInfo[2][1:].astype(float) #convert to floats

#In the reference data, Acetaldehyde was renamed to Ethanal and Ethanol was renamed to EtOH, this way all nine molecules will be solved via Madix and Ko than molecule matching between reference data and ionization data

ionizationFactorsRN2 = numpy.zeros(len(ElectronNumbers)) #initalize an array to store the ionization factors

#apply the Madix and Ko equation using the electron numbers from the reference file
for moleculeIndex in range(len(ionizationFactorsRN2)):
    ionizationFactorsRN2[moleculeIndex] = (0.6*ElectronNumbers[moleculeIndex]/14)+0.4
     

#the feature uses known ionization factors if they are available
ut.set_expected_result(ionizationFactorsRN2,expected_result_str=str(ionizationFactorsRN2),prefix=prefix,suffix=suffix)

#Round functions added so strings continue to match
#The exact outputs are below
#Expected Output is     [2.5992 4.1154 0.7695 2.0959 2.7435 3.2119 3.4891 1.4079 2.1711]
#Calculated Output was  [2.6    4.1167 0.7668 2.0952 2.7428 3.2128 3.4905 1.407  2.1703]

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