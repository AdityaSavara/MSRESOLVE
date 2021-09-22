# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 16:52:12 2018

@author: Alex
"""



'''
listCast converts objects to lists
'''
def listCast(inputObj):
    if type(inputObj) == str or type(inputObj) == int or type(inputObj) == float: #If input object is a string, int, or float, make a list of length one
        inputObj = [inputObj]
    else: #otherwise try to cast input object as list
        inputObj = list(inputObj)
    return inputObj

'''
parallelVectorize takes a list or array and creates an expanded list or array of a desired length that contains the same element in each index
'''
def parallelVectorize(inputObj, desiredLength, desiredArrayType='list'):
    if type(inputObj) == str: #strings are not compatible with this function so use listCast() first
        print("Warning: parallelVectorize is not designed to work with string inputs.")
    if len(inputObj) == desiredLength: #If the list or array is already of desired length, then just return the list or array as is
        return inputObj
    elif len(inputObj) == 1: #If the length of the list or array is 1 then create a new list or array of the desired length and fill each element with the element given in input object
        if desiredArrayType == 'list':
            parallelizedList = [] #initialize an empty list
            for index in range(desiredLength): #Loop through parallelizedList the number of times as desiredLength
                parallelizedList.append(inputObj[0]) #append the element in the input object
            return parallelizedList #Return the new list
        #TODO add array feature in parallelVectorize
        if desiredArrayType == "array":
            print("The array feature of parallelVectorize has not been implemented yet, but would not be hard to implement.")
            #return parallelizedArray
    elif len(inputObj) == 0: #blank list
        return inputObj
    else:
        raise ValueError('An expected input type was not received by parallelVectorize')
        
'''
compareElementsBetweenLists takes in two lists and their respective names
It compares elements in the first list and sees if they are in the second list
If all elements in list1 are in list2 then the dependencies have been met and the function will return None
If there is an element in list1 that is not in list2, a value error is raised informing the user that there is an element in list1 that is not in list2
'''        
def compareElementsBetweenLists(list1,list2,nameOfList1,nameOfList2):
    for i in range(0,len(list1)):
        if list1[i] not in list2:
            print(str(list1[i]) + ' in ' + nameOfList1 + ' but not in ' + nameOfList2)
            raise ValueError(str(list1[i]) + ' in ' + nameOfList1 + ' but not in ' + nameOfList2)
    return None

'''
strCheck is used to see if a variable that is supposed to be a string actually is a string
Inputs are the variable that is supposed to be a string and the name of the variable as a string
The name is default 'variable' so if the user forgets to put the variable name there it just prints 'variable must be of type string'
'''
def strCheck(stringVar,nameOfStringVar='variable'):
    if not isinstance(stringVar,str):
        raise TypeError('%s must be of type str' %(nameOfStringVar))
    return None

'''
standardzieString standardizes the stringObject that is input
It removes all whitespace on the outsides, forces to lowercase, and then takes all whitespace on the inside and standardizes it to the same type of space ' '.
Was previously in stringCompare (moved out 181011)
'''
def standardizeString(stringObject):
    import re
    #First store the strings into a variable that will be standardized
    standardizedString = stringObject
    #Strip the strings of any whitespace on the outsides
    standardizedString = standardizedString.strip()
    #Make both strings lowercase
    standardizedString = standardizedString.lower()
    #Using regex, find any style of whitespace on the inside and replace it with a standardized space
    standardizedString = re.sub('\s+',' ',standardizedString)
    
    return standardizedString

'''
stringCompare takes in two strings and compares a standardized version of the two to see if they match
Added: 181008
Last modified: 181011
'''
def stringCompare(firstString,secondString):
    standardizedFirstString = standardizeString(firstString)
    standardizedSecondString = standardizeString(secondString)
    
    #If the standardized strings match return True
    if standardizedFirstString == standardizedSecondString:
        return True
    else: #Otherwise return false
        return False
