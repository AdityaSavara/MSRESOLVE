
"""
Created on Wed Aug  1 13:47:18 2018

@author: Andrea
"""

#THE FOLLOWING LINES ARE MANDATORY FOR THE CODE
#importing the functions from UnitTesterSG module
import sys
import os
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))

import UnitTesterSG as ut
import numpy

#BELOW ARE THE LINES INTENDED TO BE CHANGED BY THE USER	
#1) import the function whose results need to be checked
import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)

#2) getting the prefix (or suffix) arugument for check_results. This is just for the output filenames.
suffix= ut.returnDigitFromFilename(__file__)
prefix=''

#3) provide the input for the function you want to test (you can also import it from a pickled object, for example)
massFragCombinations=([1,2,3,4,5],[2,1,3,4,5],[3,4,5,6,7],[9,8,7,6,5],[3,4,5,1,2])

#Initialize list variable that store the objective function output and the mass
#frag combinations.
largestMagnitudeSigFactorSumsList=[]
topMassFragCombinationsList=[]

keep_N_ValuesInSignificanceFactorCheck=5
moleculesLikelihood=numpy.array([1,0.5,1,1])
#Create reference patterns for each mass fragment combination
chosenReferenceForMassFragComb1=numpy.array([[1,2,5,50,1],
                             [2,0,0,4,0],
                             [3,1,1,40,1],
                             [4,3,2,1,0],
                             [5,0,0,0,0]])
chosenReferenceForMassFragComb2=numpy.array([[2,2,5,2,1],
                             [1,50,0,4,0],
                             [3,1,1,30,1],
                             [4,100,2,1,0],
                             [5,0,0,0,0]])
chosenReferenceForMassFragComb3=numpy.array([[3,2,5,2,1],
                             [4,0,52,4,0],
                             [5,1,96,30,1],
                             [6,3,100,1,0],
                             [7,0,0,0,0]])
chosenReferenceForMassFragComb4=numpy.array([[9,2,5,2,1],
                             [8,0,0,4,25],
                             [7,1,1,30,50],
                             [6,3,2,1,100],
                             [5,0,0,0,0]])
chosenReferenceForMassFragComb5=numpy.array([[3,2,5,2,1],
                             [4,0,0,4,0],
                             [5,1,1,30,1],
                             [1,3,2,1,0],
                             [2,0,0,0,0]])
#Create a list that makes the reference patterns iterable
refPatternList=[chosenReferenceForMassFragComb1,chosenReferenceForMassFragComb2,chosenReferenceForMassFragComb3,chosenReferenceForMassFragComb4,chosenReferenceForMassFragComb5]

#4) get the output of the function, which is what will typically be checked.
#In it's intended use, the funciton is to append to a list 
#(topSignificanceFactorCheckList). This list was initialized and is to be 
#appended to in a loop
for massFragCombinationIndex, massFragCombination in enumerate(massFragCombinations):
    refIntensity=refPatternList[massFragCombinationIndex][:,1:]
    [largestMagnitudeSigFactorSumsList,topMassFragCombinationsList, valuesStoredInSFTopList]=MSRESOLVE.significanceFactorCheck(refIntensity,largestMagnitudeSigFactorSumsList,topMassFragCombinationsList, massFragCombination, keep_N_ValuesInSignificanceFactorCheck, moleculesLikelihood)

#The output result is the best mass fragment according to the objective function
resultObj= topMassFragCombinationsList[0] #, output[1], output[2]]  #You can alternatively populate resultObj with whatever you want, such as a list.

#5) A string is also typically provided, but is an optional argument. You can provide whatever string you want.
resultStr= str(resultObj)

#Set expected results to be the best mass fragment combination
ut.set_expected_result(massFragCombinations[3],expected_result_str=str(massFragCombinations[3]),suffix=suffix)

#6) Checking the result of the function using check_results. In this case the result is sumList1 object. 

#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=1E-5, absoluteTolerance=1E-8)
    
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True, relativeTolerance=1E-5, absoluteTolerance=1E-8)
