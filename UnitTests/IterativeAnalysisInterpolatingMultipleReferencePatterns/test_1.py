import os
import sys
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#import the functions from UnitTesterSG
import UnitTesterSG as ut
import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)
import DefaultUserInput as G, importlib; importlib.reload(G) #This is needed because we need the __var_list__
MSRESOLVE.G = G #This is because we need to overwrite whatever user input the user has with the default user input.
MSRESOLVE_var_list = G.__var_list__ #need to store this to reassign in the new namespace.
    
#get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)
#prefix. Make this '' if you do not want any prefix.
prefix = ''

import numpy as np
import copy
import shutil


#this helper function will be used to append the new settings in each iteration in this test.
def appendToIterFile(iterNumber, listOfStringsToAppend):
    iterNumber = str(iterNumber) #in case someone puts in an integer.
    os.chdir('_iter_'+iterNumber)
    with open("UserInput_iter_"+iterNumber+".py", "a") as iterInputFile: #the a is to append.
        iterInputFile.write("\n") #we need to first start a new line.
        for stringToWrite in listOfStringsToAppend:
            iterInputFile.write(stringToWrite + "\n")
    os.chdir("..") #now go back to regular directory, again.


#1st I need to delete the old iteration directories if they exist.
#First make a directory list.
listOfDirectoriesAndFiles = os.listdir(".")
directoryList = []
for elem in listOfDirectoriesAndFiles:
    if os.path.isdir(elem) == True:
        directoryList.append(elem)
#Then delete each one.
for directory in directoryList:
    shutil.rmtree(directory)
    
#Also delete the old "ScaledConcentrations.csv" and "TotalConcentrations.csv"
if os.path.isfile('ScaledConcentrations.csv'): #check if file is there
    os.remove('ScaledConcentrations.csv') #remove it if present
if os.path.isfile('TotalConcentrations.csv'):
    os.remove('TotalConcentrations.csv')

#NON ITERATIVE WAY FIRST.    
MSRESOLVE_var_list = G.__var_list__ #need to store this to reassign in the new namespace.
import test_1_initial_input_noniterative
MSRESOLVE.G = test_1_initial_input_noniterative
MSRESOLVE.G.__var_list__ = MSRESOLVE_var_list #need to repopulate var list since namespace was re-assigned.
MSRESOLVE.main()

#NOW DO THINGS THE ITERATIVE WAY.

#ITERATION 1
#now populate the globals with our first iteration's userinput.
MSRESOLVE_var_list = G.__var_list__ #need to store this to reassign in the new namespace.
import test_1_initial_input_iterative
MSRESOLVE.G = test_1_initial_input_iterative
MSRESOLVE.G.__var_list__ = MSRESOLVE_var_list #need to repopulate var list since namespace was re-assigned.
#now run MSRESOLVE.py for the first iteration.
MSRESOLVE.main()

#ITERATION 2
#now we need to change the choices for iteration 2. 
#we 1st need to make the strings we want to append into a list.
#Then call main function again.
listOfStringsToAppend = ["chosenMoleculesNames = ['(E) 2-Butenal (Crotonaldehyde']", "chosenMassFragments = [39]"]
appendToIterFile(2, listOfStringsToAppend)
MSRESOLVE.main()

#ITERATION 3
#now ready for iteration three.
listOfStringsToAppend = ["chosenMoleculesNames = ['Ethanol']", "chosenMassFragments = [31]"]
appendToIterFile(3, listOfStringsToAppend)
MSRESOLVE.main()

#ITERATION 4
#now ready for iteration three.
listOfStringsToAppend = ["chosenMoleculesNames = ['H2O']", "chosenMassFragments = [18]"]
appendToIterFile(4, listOfStringsToAppend)
MSRESOLVE.main()

#ITERATION 5
#now ready for iteration three.
listOfStringsToAppend = ["chosenMoleculesNames = ['Acetaldehyde']", "chosenMassFragments = [41]"]
appendToIterFile(5, listOfStringsToAppend)
MSRESOLVE.main()

#ITERATION 6
#now ready for iteration three.
listOfStringsToAppend = ["chosenMoleculesNames = ['Ethylene (Ethene)']", "chosenMassFragments = [27]"]
appendToIterFile(6, listOfStringsToAppend)
MSRESOLVE.main()

#ITERATION 7
#now ready for iteration three.
listOfStringsToAppend = ["chosenMoleculesNames = ['CO', 'CO2', 'H2']", "chosenMassFragments = [2, 28, 44]"]
appendToIterFile(7, listOfStringsToAppend)
MSRESOLVE.main()

#ITERATION 8 does not occur because there are no remaining molecules.

arrayReadFromNonIterativeAnalysis = np.genfromtxt("ScaledConcentrations.csv", skip_header = 1, delimiter=",", unpack=True)
arrayReadFromIterativeAnalysis = np.genfromtxt("TotalConcentrations.csv", skip_header = 1, delimiter=",", unpack=True)
# now I am going to rearrange the  columns  in a way that I know they're going to match.
# I know from previous experience that for this set of files Scaled Concentrations has this order:
# Time,
#Acetaldehyde Concentration Relative to CO,
#(E) 2-Butenal (Crotonaldehyde Concentration Relative to CO,
#CO Concentration Relative to CO,
#CO2 Concentration Relative to CO,
#Ethylene (Ethene) Concentration Relative to CO,
#Ethanol Concentration Relative to CO,
#Crotyl Alcohol Concentration Relative to CO,
#H2 Concentration Relative to CO,
#H2O Concentration Relative to CO
#
#In contrast, Total Concentrations has this order:
# Time,
#Crotyl Alcohol Concentration Relative to CO,
#(E) 2-Butenal (Crotonaldehyde Concentration Relative to CO,
#Ethanol Concentration Relative to CO,
#H2O Concentration Relative to CO,
#Acetaldehyde Concentration Relative to CO,
#Ethylene (Ethene) Concentration Relative to CO,
#CO Concentration Relative to CO,
#CO2 Concentration Relative to CO,
#H2 Concentration Relative to CO
#
# so to rearrange the columns, I'm going to make  a deep copy of arrayReadFromNonIterativeAnalysis, and populate it.
arrayToFillFromIterativeAnalysis = copy.deepcopy(arrayReadFromNonIterativeAnalysis)
arrayToFillFromIterativeAnalysis[0] = arrayReadFromIterativeAnalysis[0]
arrayToFillFromIterativeAnalysis[1] = arrayReadFromIterativeAnalysis[5]
arrayToFillFromIterativeAnalysis[2] = arrayReadFromIterativeAnalysis[2]
arrayToFillFromIterativeAnalysis[3] = arrayReadFromIterativeAnalysis[7]
arrayToFillFromIterativeAnalysis[4] = arrayReadFromIterativeAnalysis[8]
arrayToFillFromIterativeAnalysis[5] = arrayReadFromIterativeAnalysis[6]
arrayToFillFromIterativeAnalysis[6] = arrayReadFromIterativeAnalysis[3]
arrayToFillFromIterativeAnalysis[7] = arrayReadFromIterativeAnalysis[1]
arrayToFillFromIterativeAnalysis[8] = arrayReadFromIterativeAnalysis[9]
arrayToFillFromIterativeAnalysis[9] = arrayReadFromIterativeAnalysis[4]

#before comparing the arrays, we now need to round them.
roundedArrayReadFromNonIterativeAnalysis  = np.round(arrayReadFromNonIterativeAnalysis, decimals=4)
roundedArrayFilledFromIterativeAnalysis = np.round(arrayToFillFromIterativeAnalysis, decimals=4)

ut.set_expected_result(roundedArrayReadFromNonIterativeAnalysis,expected_result_str=str(roundedArrayReadFromNonIterativeAnalysis), prefix=prefix,suffix=suffix)

resultObj = roundedArrayFilledFromIterativeAnalysis

#String must be provided provided. Make it '' if you do not want to use a result string.
resultStr = str(roundedArrayFilledFromIterativeAnalysis)

relativeTolerance = 1E-4
absoluteTolerance = 1E-8


#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
