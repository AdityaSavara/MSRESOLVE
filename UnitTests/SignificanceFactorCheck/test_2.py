
"""
Created on Wed Aug  1 13:47:18 2018

@author: Andrea
"""

#THE FOLLOWING LINES ARE MANDATORY FOR THE CODE
#importing the functions from UnitTesterSG module
import sys
sys.path.insert(1, "..\\lib")
sys.path.insert(1, "..")
sys.path.insert(1, "..\..")
import UnitTesterSG as ut
import numpy

#BELOW ARE THE LINES INTENDED TO BE CHANGED BY THE USER	
#1) import the function whose results need to be checked
import MSRESOLVE
#2) getting the prefix (or suffix) arugument for check_results. This is just for the output filenames.
suffix= ut.returnDigitFromFilename(__file__)
prefix=''
#3) provide the input for the function you want to test (you can also import it from a pickled object, for example)
massFragCombinations=([1,2,3,4,5],[2,1,3,4,5],[3,4,5,6,7],[9,8,7,6,5],[3,4,5,1,2])

largestMagnitudeSigFactorSumsList=[]

keep_N_ValuesInSignificanceFactorCheck=3
moleculesLikelihood=numpy.array([1,0.5,1,1])
chosenReference=numpy.array([[1,2,5,2,1],
                             [2,0,0,4,0],
                             [3,1,1,30,1],
                             [4,3,2,1,0],
                             [5,0,0,0,0]])

#4) get the output of the function, which is what will typically be checked.
for massFragCombinationIndex, massFragCombination in enumerate(massFragCombinations):
    chosenReference[massFragCombinationIndex,massFragCombinationIndex]=54#This is to change the array so the significance factor sums will be different during each iteration. If a different array is not provided each time, the significance factor sums will all be the same.
    #calculates the significance factor for each element in the chosen reference array and sums all of the significane values for the whole array. It keeps the mass fragments that have the largest magnitude of significance sum.
    [largestMagnitudeSigFactorSumsList, valuesStoredInSFTopList]=MSRESOLVE.significanceFactorCheck(chosenReference[:,1:],largestMagnitudeSigFactorSumsList,keep_N_ValuesInSignificanceFactorCheck,massFragCombination, moleculesLikelihood)
#print(output)
resultObj= [largestMagnitudeSigFactorSumsList, valuesStoredInSFTopList] #, output[1], output[2]]  #You can alternatively populate resultObj with whatever you want, such as a list.
#5) A string is also typically provided, but is an optional argument. You can provide whatever string you want.
resultStr= str(resultObj)
#6) Checking the result of the function using check_results. In this case the result is sumList1 object. 

#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False)
    
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True)
