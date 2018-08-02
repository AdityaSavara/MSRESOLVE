
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

#BELOW ARE THE LINES INTENDED TO BE CHANGED BY THE USER	
#1) import the function whose results need to be checked
import MSRESOLVE
#2) getting the prefix (or suffix) arugument for check_results. This is just for the output filenames.
suffix= ut.returnDigitFromFilename(__file__)
prefix=''
#3) provide the input for the function you want to test (you can also import it from a pickled object, for example)
massFragCombinations=([1,2,3,4,5],[2,1,3,4,5],[3,4,5,6,7],[9,8,7,6,5],[3,4,5,1,2])

rowSumsList=([1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1])

topRoughUniquenessSumsList=[]

keep_N_ValuesInRoughUniquenessCheck=4

#4) get the output of the function, which is what will typically be checked.
for counter, massFragCombination in enumerate(massFragCombinations):
    #calculates a sum that roughly expresses how unique the molecular mass fragments are to the different molecules, but this is a quick and not-rigrous method. Then, the value is stored *only* if it is in the top N of the values so far.
    [topRoughUniquenessSumsList,valueStoredInRUTopList] = MSRESOLVE.roughUniquenessCheck(rowSumsList[counter], topRoughUniquenessSumsList, keep_N_ValuesInRoughUniquenessCheck, massFragCombination)

#The output of the function is sightly unexpected, the last mass fragment combinaiton is stored instead of the fragments [9,8,7,6,5].
#This is due to the bisect search used. When the values are all the same, it starts looking in the fragment array. This is not a problem
#for the intended use of the objective funciton since fragment combinations will be evaluated in order.

resultObj= [topRoughUniquenessSumsList,valueStoredInRUTopList] #, output[1], output[2]]  #You can alternatively populate resultObj with whatever you want, such as a list.
#5) A string is also typically provided, but is an optional argument. You can provide whatever string you want.
resultStr= str(resultObj)
#6) Checking the result of the function using check_results. In this case the result is sumList1 object. 

#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False)
    
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = True)