# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 11:40:09 2018

@author: 3cw
"""

import os
import filecmp
import shutil


def CompareTwoFiles(fileName1, fileName2, filePath):
    os.chdir(filePath)
    comparison = filecmp.cmp(fileName1,fileName2)
    return comparison

def test_Run(pytest = True):
    result = CompareTwoFiles('DefaultUserInput.py', 'UserInput.py', '../..' )
    if pytest:
        assert result == True
        
    if not pytest: 
        if result == True:
            print("DefaultUserInput.py and UserInput.py do match" )
        if result == False:
            overwriteOption=str(input('DefaultUserInput.py and UserInput.py do NOT match. Would you like to overwrite DefaultUserInput.py (Y or N)?'))
            if overwriteOption == 'Y':
                shutil.copy('UserInput.py', 'DefaultUserInput.py')
            elif overwriteOption != 'N':
                print("UserError: Only Y or N allowed. Please run program again.")
    
if __name__ == "__main__":
   test_Run(pytest = False)