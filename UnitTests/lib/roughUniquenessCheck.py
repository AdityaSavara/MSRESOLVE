# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 13:54:25 2018

@author: Andrea
"""
import numpy
import bisect

#Store and Pop takes in a sorted list. Using a binary search method, it will 
#find the location where the value to insert will be inserted(if possible). The
#value will be inserted there and the last value removed from the list (if
#applicable). The max lengh of list is used to ensure the list does not go over
#a certain size. The bisect method used sorts the list from lowest to highest.
def storeAndPop(sortedList, valueToInsert, maxLengthOfList):
    #Find the insertion idex where the value will be inserted by using a binary
    #search
    insertionIndex=bisect.bisect(sortedList, valueToInsert)
    #Initialize a variable to keep track of if a value was inserted into the
    #list.
    valueStoredInList=False

    #If the list isn't yet filled, the value will inherently be in the top N
    #value in the list. This vlaue can just be inserted at the insertionIndex.
    if len(sortedList)<maxLengthOfList:    
        sortedList.insert(insertionIndex, valueToInsert)
        valueStoredInList=True
    #If the list already contains N elements, a new element could either be 
    #inserted in the list or at the end of the list. Because the list is 
    #already at its maximum length, nothing shouold be added to the end. This
    #check is to make sure nothing is going to be added to the end.
    elif insertionIndex<maxLengthOfList:
        #insert the value to insert in the location found through the binary
        #search
        sortedList.insert(insertionIndex, valueToInsert)
        valueStoredInList=True
        #delete the last element since somthing was added to the list
        del sortedList[-1]
    return sortedList, valueStoredInList

#The rough uniqueness check is a limiting check that takes the mass fragment
#combinations that pass the row sums check and builds a list of 
#keep_N_ValuesInRoughUniquenessCheck that contain the largest number of zeros.
#The largest number of zeros would be most likely to pass the SLS method.
def roughUniquenessCheck(rowSumsList, topRoughUniquenessCheckList, keep_N_ValuesInRoughUniquenessCheck, massFragCombination):
    
    #We want to save the smallest sum since that would contain the smallest 
    #number of zeros.
    roughUniqueness=numpy.sum(rowSumsList) 
    
    #Create a tuple that stores the rough uniqueness value and the 
    #massFragCombination
    roughUniquenessTuple=tuple([roughUniqueness, massFragCombination])
    
    #Use Store and Pop to add the tuple to the list of top rough uniquness 
    #combinations. This will only save a user specified number of tuples.
    [topRoughUniquenessCheckList, valueStoredInRUTopList]=storeAndPop(topRoughUniquenessCheckList, roughUniquenessTuple, keep_N_ValuesInRoughUniquenessCheck)
    return topRoughUniquenessCheckList, valueStoredInRUTopList