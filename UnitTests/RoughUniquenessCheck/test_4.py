
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

#BELOW ARE THE LINES INTENDED TO BE CHANGED BY THE USER	
#1) import the function whose results need to be checked
import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)
import numpy
#2) getting the prefix (or suffix) arugument for check_results. This is just for the output filenames.
suffix= ut.returnDigitFromFilename(__file__)
prefix=''
#3) provide the input for the function you want to test (you can also import it from a pickled object, for example)
massFragCombinations=([1,2,3,4,5],[2,1,3,4,5],[3,4,5,6,7],[9,8,7,6,5],[3,4,5,1,2])

#Create reference patterns as 1s and zeros for each mass fragment combination.
#The referecne patterns as ones and zeros shown alreadt pass the two rough checks
#for SLS solvability since the checks are used before the roughuniquness check
#in order to prevent unnecessary calculations.
chosenReferenceForMassFragComb1=numpy.array([[1,1,1,1,1],
                             [1,0,0,1,0],
                             [1,1,1,0,1],
                             [1,1,1,1,0],
                             [1,0,1,1,1]])
chosenReferenceForMassFragComb2=numpy.array([[1,1,1,1,1],
                             [1,0,0,1,0],
                             [1,1,1,0,1],
                             [1,1,1,1,0],
                             [1,0,1,1,1]])
chosenReferenceForMassFragComb3=numpy.array([[1,1,1,1,1],
                             [1,0,0,1,0],
                             [1,1,1,0,1],
                             [1,1,1,1,0],
                             [1,0,1,1,1]])
chosenReferenceForMassFragComb4=numpy.array([[1,1,1,1,1],
                             [1,0,0,1,0],
                             [1,1,1,0,1],
                             [1,1,1,1,0],
                             [1,0,1,1,1]])
chosenReferenceForMassFragComb5=numpy.array([[1,1,1,1,1],
                             [1,0,0,1,0],
                             [1,1,1,0,1],
                             [1,1,1,1,0],
                             [1,0,1,1,1]])
#Create a list that makes the reference patterns iterable
refPatternList=[chosenReferenceForMassFragComb1,chosenReferenceForMassFragComb2,chosenReferenceForMassFragComb3,chosenReferenceForMassFragComb4,chosenReferenceForMassFragComb5]

#Generate the sums across the mass fragments for each molecule.
#The slicing occurs because the first element in each row represents the mass fragment number. The fragment number should not be kept during calculations.
#The mass fragment numbers are kept to maintian consistency with the fuction they are added to.
rowSumsList=[]
for refPattern in refPatternList:
    refIntensity=refPattern[:,1:]
    rowSumsList.append(numpy.sum(refIntensity, axis=0))

topRoughUniquenessSumsList=[]
topMassFragCombinationsList=[]

keep_N_ValuesInRoughUniquenessCheck=4

#4) get the output of the function, which is what will typically be checked.
for massFragCombinationIndex, massFragCombination in enumerate(massFragCombinations):
   #calculates a sum that roughly expresses how unique the molecular mass fragments are to the different molecules, but this is a quick and not-rigrous method. Then, the value is stored *only* if it is in the top N of the values so far.
    [topRoughUniquenessSumsList,topMassFragCombinationsList,valueStoredInRUTopList] = MSRESOLVE.roughUniquenessCheck(rowSumsList[massFragCombinationIndex], topRoughUniquenessSumsList,topMassFragCombinationsList, keep_N_ValuesInRoughUniquenessCheck, massFragCombination)

resultObj= topMassFragCombinationsList[0] #, output[1], output[2]]  #You can alternatively populate resultObj with whatever you want, such as a list.
#5) A string is also typically provided, but is an optional argument. You can provide whatever string you want.
resultStr= str(resultObj)

#Set expected results to be the best mass fragment combination
ut.set_expected_result(massFragCombinations[0], expected_result_str=str(massFragCombinations[0]),suffix=suffix)

#6) Checking the result of the function using check_results. In this case the result is sumList1 object. 

#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False)
    
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True)
