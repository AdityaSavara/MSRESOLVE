
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
import UnitTesterSG 
import numpy

#BELOW ARE THE LINES INTENDED TO BE CHANGED BY THE USER	
#1) import the function whose results need to be checked
import roughUniquenessCheck as ruc
#2) getting the prefix (or suffix) arugument for check_results. This is just for the output filenames.
suffix= UnitTesterSG.returnDigitFromFilename(__file__)
#3) provide the input for the function you want to test (you can also import it from a pickled object, for example)
massFragCombinations=([1,2,3,4,5],[2,1,3,4,5],[3,4,5,6,7],[9,8,7,6,5],[3,4,5,1,2])

rowSumsList=([1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1])

topRoughUniquenessCheckList=[]

keep_N_ValuesInRoughUniquenessCheck=4

#4) get the output of the function, which is what will typically be checked.
for counter, massFragCombination in enumerate(massFragCombinations):
 [topRoughUniquenessCheckList,valueStoredInRUTopList] = ruc.roughUniquenessCheck(rowSumsList[counter], topRoughUniquenessCheckList, keep_N_ValuesInRoughUniquenessCheck, massFragCombination)
#print(output)
resultObj= [topRoughUniquenessCheckList,valueStoredInRUTopList] #, output[1], output[2]]  #You can alternatively populate resultObj with whatever you want, such as a list.
#5) A string is also typically provided, but is an optional argument. You can provide whatever string you want.
resultStr= str(resultObj)
#6) Checking the result of the function using check_results. In this case the result is sumList1 object. 

#run the Unit Tester
def test_Run(allowOverwrite = False):
    #if the user wants to be able to change what the saved outputs are
    if allowOverwrite:
        #This function call is used when this test is run solo as well as by UnitTesterSG
        UnitTesterSG.check_results(resultObj, resultStr, prefix = '', suffix=suffix)
    #this option allows pytest to call the function
    if not allowOverwrite: 
        #this assert statement is required for the pytest module 
        assert UnitTesterSG.check_results(resultObj, resultStr, prefix = '', suffix=suffix, allowOverwrite = False) == True
    
if __name__ == "__main__":
   test_Run(allowOverwrite = True)