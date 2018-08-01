
"""
Created on Wed Aug  1 14:12:07 2018

@author: Andrea
"""
import MSRESOLVE
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

#The significance factor check is a limiting check the selects the mass 
#fragment combinations having the largest sum of significance factors. It
#calculates the significance factors for each element in the sub-reference 
#array (refernce array with zeros for all the data for mass fragments that 
#aren't needed). Is then takes the sum of all of the significance factors.
def significanceFactorCheck(chosenReference ,topSignificanceFactorCheckList, keep_N_ValuesInSignificanceFactorCheck, massFragCombination, moleculesLikelihood):
    
    #Initialize the sum of the significance factors
    sigFactorSum=0
    
    #The elemSignificanceCalculator finds the significance factors a column at 
    #a time. This loops through the columns(molecules)
    for columnCounter in range(len(chosenReference)-2):
        
        #Finds the significance factors for each element in the column. This 
        #does not include the first column since that is the mass fragment 
        #numbers
        significanceColumnDataList=MSRESOLVE.ElemSignificanceCalculator(chosenReference[:,1:], columnCounter, moleculesLikelihood)
        
        #Sums the significance factors across the column and subtracts from the
        #sum for the whole ref data array. The subtraction is used to make the
        #sum negative. The binary search used will order from lowest to highest
        #The largest magnitude will then be the smallest number and be kept
        #during the store and pop function
        sigFactorSum-=sum(significanceColumnDataList)
        
        ####Currently there is no real need to maintain a significance data 
        ####list for the whole array
        
    #Creates a tuple that stores the significane factor sum and the mass
    #fragment combination
    sigFactorTuple=tuple([sigFactorSum, massFragCombination])
    
    #Uses store and pop to maintian a list of the mass fragment with the
    #largest significance factors.
    [topSignificanceFactorCheckList,valueStoredInSFTopList]=storeAndPop(topSignificanceFactorCheckList,sigFactorTuple,keep_N_ValuesInSignificanceFactorCheck)
    
    return topSignificanceFactorCheckList, valueStoredInSFTopList