# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 14:16:13 2018

@author: 3cw
"""
#THE FOLLOWING LINES ARE MANDATORY FOR THE CODE
#importing the functions from UnitTesterSG module
import sys
import os
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
import UnitTesterSG 
import numpy

#BELOW ARE THE LINES INTENDED TO BE CHANGED BY THE USER	
#1) import the function whose results need to be checked
from MSRESOLVE import DistinguishedArrayChooser
#2) getting the prefix (or suffix) arugument for check_results. This is just for the output filenames.
suffix= UnitTesterSG.returnDigitFromFilename(__file__)
#3) provide the input for the function you want to test (you can also import it from a pickled object, for example)
refMassFrags=numpy.array([[10, 3, 2],
                          [2, 1, 4],
                          [6, 11, 7],
                          [13, 5, 4]])
    
correctionValues=numpy.array([[1, 1, 1],
                              [2, 2, 2],
                              [3, 3, 3],
                              [4, 4, 4]])
    
rawSignals=numpy.array([[10],
                        [20],
                        [30],
                        [40]])
                  
moleculesLikelihood = [1,0.2]
    
sensitivityValues = [1,1]    

#4) get the output of the function, which is what will typically be checked. 
output = DistinguishedArrayChooser(refMassFrags,correctionValues,rawSignals, moleculesLikelihood, sensitivityValues)
#print(output)
resultObj= output #, output[1], output[2]]  #You can alternatively populate resultObj with whatever you want, such as a list.
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
