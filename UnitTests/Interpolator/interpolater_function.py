# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 11:22:59 2017

@author: tienhung2501
"""
import numpy
import numpy as np
import math


###OLD interpolater function. The new one is below.
#### The comments in this function are no longer completely accurate since we have added some absolute values and deltas where they were not before.
#### ( Look at the version of the program XXXXXX to see what this function looked like previously )
# this function makes sure that the data cannot change by more than a factor of two from one time to another
#this is done by finding those lines that have an issue and then getting the average of those two lines and 
#inserts this average row in between these two respective rows, then putting averages up until there is no jump 
#greater than 100%
def Interpolater(data,abscissa,dataRangeSpecifierYorN = [],datafromcsv = [],MaxAllowedDeltaYRatio=2.0, IgnorableDeltaYThreshold = 0.0001):
#### The comments in this function are no longer completely accurate since we have added some absolute values and deltas where they were not before.
#### ( Look at the version XXXXXX of the program to see what this function looked like previously )
#The current Interpolater function, while functional, suffers from two serious performance issues. 
#
#First, the use of repeated insert statements is slow, because it requires that the array be reconstructed for each insertion.
#This is because arrays are of fixed length. Solution: As the array is checked, maintain a new list of points that includes 
#both the original data and any added points using 'list'.append(). At the end of the program, convert this list to an array
#and save it under the name of the old array.
#
#Second, the current loops instruct the program to check every jump and insert an average if necessary. 
#However, this means that if there were, for example, a 40 fold jump, the program would have to iterate through 
#all points 5 times to resolve the problem. Solution: Resolve each jump completely in turn, by maintaining 
#a pointer to the base of the jump until enough pointers have been added. 
#In this way, the program could only iterate once through the data.


    #this first loop doesn't allow for jumps greater than 100%, the line saying it cannot equal more than the 
    #previous term*2 The place_holder object allows the function to know how many times the loop has added a row,
    #thereby knowing the size of the new array. the if statement doesn't allow the last row_counter to be checked 
    #since there is not going to be a row beyond the length of the array
    
    for column_counter in range(len(data[0,:])):#array-indexed for loop
        place_holder = 0 #reset after each column
        place_holder2 = 0
        for row_counter in range(len(data[:,0])):#array-indexed for loop
            if row_counter != len(data[:,0]) - (1+place_holder): #we don't want to go past the last row, because every row uses the average with the one after it
                if abs(data[row_counter+1+place_holder,column_counter]-data[row_counter+place_holder,column_counter]) > abs(data[row_counter+place_holder,column_counter]*MaxAllowedDeltaYRatio):#If one row has a changes of greater than 100% of the last, then an average is put between them


                    while abs(data[row_counter+1+place_holder,column_counter]-data[row_counter+place_holder,column_counter]) > abs(data[row_counter+place_holder,column_counter]*MaxAllowedDeltaYRatio):#If one row has a changes by a factor greater than 2, then an average is put between that and the next row, and averages keep getting placed there until there is no more jump that large
                        if abs(data[row_counter+1+place_holder,column_counter] - data[row_counter+place_holder,column_counter]) > IgnorableDeltaYThreshold:#the additions are only made if the change is greater then 0.0001, (no tiny changes)
                            data = numpy.insert(data,row_counter+place_holder+1,sum(data[(row_counter+place_holder):(row_counter+place_holder+2),:])/2,axis = 0)
                            abscissa = numpy.insert(abscissa,row_counter+place_holder+1,sum(abscissa[(row_counter+place_holder):(row_counter+place_holder+2)])/2)
                            if dataRangeSpecifierYorN == 'yes': #this is for the datarangespecifier function, to make the csv file as long as the input csv file
                                datafromcsv = numpy.insert(datafromcsv,row_counter+place_holder+1,sum(datafromcsv[(row_counter+place_holder):(row_counter+place_holder+2),:])/2,axis = 0)
                            place_holder2 = place_holder2 + 1
                            
                        else:
                            break
                    place_holder = place_holder2 #keeps the array in check since many rows are being added (for all rows added the index is changed)

                    
    #this loop is very similar to the first, but it makes sure that there is no decrease greater than 50 %, a factor
    #of 2 just as above, so that values cannot be less than half the value before, furthermore, the place_holder changes
    #a little bit- in this way- it changes the index of the two values that get checked, since the average can't be more 
    #than 50% below the higher value, you must check the average versus the lower value - which means changing the index
    #in the previous loop, this index grew but did not affect its own loop, but only the outer loops, in this nested loop,
    #the index actually affects the loop itself
    
    for column_counter in range(len(data[0,:])):#array-indexed for loop
        place_holder = 0
        place_holder2 = 0
        for row_counter in range(len(data[:,0])):#array-indexed for loop
            if row_counter != len(data[:,0]) - (1+place_holder):#same comments as above, symmetrical code
                if abs(data[row_counter+1+place_holder,column_counter]-data[row_counter+place_holder,column_counter]) > abs(data[row_counter+place_holder,column_counter]/MaxAllowedDeltaYRatio):
                    while abs(data[row_counter+1+place_holder,column_counter]-data[row_counter+place_holder,column_counter]) > abs(data[row_counter+place_holder,column_counter]/MaxAllowedDeltaYRatio):
                        if abs(data[row_counter+place_holder,column_counter] - data[row_counter+1+place_holder,column_counter]) > IgnorableDeltaYThreshold:
                            data = numpy.insert(data,row_counter+place_holder+1,sum(data[(row_counter+place_holder):(row_counter+place_holder+2),:])/2,axis = 0)
                            abscissa = numpy.insert(abscissa,row_counter+place_holder+1,sum(abscissa[(row_counter+place_holder):(row_counter+place_holder+2)])/2)
                            if dataRangeSpecifierYorN == 'yes':
                                datafromcsv = numpy.insert(datafromcsv,row_counter+place_holder+1,sum(datafromcsv[(row_counter+place_holder):(row_counter+place_holder+2),:])/2,axis = 0)
                            place_holder = place_holder + 1
                        else:
                            break
    return data,abscissa,datafromcsv


#Commented out because no longer in use (slower than analytical Linear interpolater)
#Uses a 1st degree polynomial fit to linearly interpolate values between 2 array-like objects or single values based on a point inserted between two other values. 
#def littleLinearInterpolator(array1,array2, xToInterpolateTo, interpolateBetweenValue1,interpolateBetweenValue2):
#    
#    def oneDegreeFit(yData,xData):
#        fitObject=numpy.polyfit(xData,yData,1)
#       return fitObject
#
#    #Puts the values into arrays to be used by polyfit()
#    xDataSingleArray=(interpolateBetweenValue1, interpolateBetweenValue2)
#    yDataListOfArrays=(array1,array2)
#    
#    #Generates the coeffecients created through a linear regression of the data provided to it
#    vectorOfFitObjects=oneDegreeFit(yDataListOfArrays,xDataSingleArray)
#    
#    #Generates an array or single value of the a y-point(s) that corelates with the x value to be interpolated to. 
#    vectorOfInterpolatedValues=numpy.polyval(vectorOfFitObjects,xToInterpolateTo)
#    
#    return vectorOfInterpolatedValues

def analyticalLinearInterpolator (ordinateValues1, ordinateValues2, abscissaValueToInterpolateTo, abscissaValue1, abscissaValue2):
        slopesVector=(ordinateValues2-ordinateValues1)/(abscissaValue2-abscissaValue1)
        interceptVector=ordinateValues1-slopesVector*abscissaValue1
        interpolatedOrdinateValues=slopesVector*abscissaValueToInterpolateTo+interceptVector
        return interpolatedOrdinateValues
    
#This is a helper function made to find the sign of a certain element and return a multiplication factor containing the sign of the element.
#It returns -1 if the value is negative and 1 if the value is positive.
def signOfElement(element):
    sign=1.0
    if (element<0):
        sign=-1.0
    return sign

#For any zero or near-zero value, this funciton inserts rows above and below the row containing the zero/near-zero contianing half of the ignorableDeltaYThreshold directly above
#or below the value in question to prevent the interpolator from interpolating between the two points.
def adjustForNearZerosAndZeros(data, abscissa, IgnorableDeltaYThreshold):
    
    #Step 3.a
    #The nearZeroIndicesArray is created to recod the location of any zeros in the original data. A one represents a non-zero greater than the ignorableDeltaYThreshold while a 
    #zero represents a zero/nearZero This is used to find the locations where rows need to be added in this pre-processing step. It is not used outside of this function.
    nearZeroIndicesArray=(abs(data)>=IgnorableDeltaYThreshold)/1.0
    
    
    #The following statements convert the arrays into lists that would be more effecient for inserting rows
    nearZeroIndicesList=list(nearZeroIndicesArray)
    dataList=list(data)
    abscissaList=list(abscissa)
    
    #Since for loops do not revaluate the len(dataList) everytime it iterates, a number of rows added has to be included that ensures the commands contianed within
    #the loop only operate on the original rows in the dataList.
    rowsAdded=0
    
    #Step 3.b.1
    #This for loop iterates throught the rows of the original dataList and adds rows above and below if a zero or near-zero value is included in the data. The locations
    #of the zeros are recorded in the zeroIndicesArray, while the locations of zeros and near-zeros are kept in the nearZeroIndicesArray. The rows that are isnerted
    #have half of the threshold value inserted above and below the zero/nearZero with the other values in the row being linearly interpolated.
    for rowCounter in range(len(dataList)):
        
        #Keeps track of the current_row and abscissa for the row that is being evaluated
        current_row=dataList[rowCounter+rowsAdded]
        current_row_abscissa=abscissaList[rowCounter+rowsAdded]
        
        #Lists of the interpolated abscissa values that need to be inserted above or below a row based on the inserted half-threshold value.
        abscissaToInsertBelow=[]
        abscissaToInsertAbove=[]
        
        #Array of locations of zeros/nearZeros in the current row that need rows inserted above and/or below
        nearZeroIndicesRow=nearZeroIndicesList[rowCounter]
        
        #This section inserts the rows above if it is not the first row.
        if rowCounter !=0:
            
            previous_row=dataList[rowCounter+rowsAdded-1]
            previous_row_abscissa=abscissaList[rowCounter+rowsAdded-1]
            
            #Step 3.b.i
            #Loops through the values in the nearZeroIndicesRow to determine if any zeros/ are present in the row
            for columnCounter, nearZeroIndexValue in enumerate (nearZeroIndicesRow):
            
                absDifferenceAbove=abs(current_row[columnCounter]-previous_row[columnCounter])
            
                #If a zero/nearZero is present AND the difference between the current row and the above row is significant(greater than the IgnorableDeltaYThreshold)
                #a new abscissa will be calculated based on the halfMinThreshold to be inserted.
                if nearZeroIndexValue==0 and absDifferenceAbove>IgnorableDeltaYThreshold:
                
                    #Makes the sign of the halfMinThreshold the same as the one above.
                    halfMinThresholdAbove=(IgnorableDeltaYThreshold/2)*signOfElement(previous_row[columnCounter])
                    #Interpolates the new abscissa value
                    interpolatedAbscissaValue=analyticalLinearInterpolator(previous_row_abscissa,current_row_abscissa, halfMinThresholdAbove, previous_row[columnCounter], current_row[columnCounter])
                    abscissaToInsertAbove.append(interpolatedAbscissaValue)
            
            #Step 3.b.ii
            #Sorts the abscissa values so they are isnerted in the right order (increasing) and the rest of the row can be interpolated based off of the abscissa value.
            abscissaToInsertAbove=sorted(abscissaToInsertAbove)
        
            
            #Iterate through the interpolated abscissa values that need to be inserted
            for abscissaToInsertValue in abscissaToInsertAbove:
                
                #Step 3.b.iii
                #Inserts the new row in the dataList containing the interpolated values that correspond to the value in the abscissa.
                interpolatedDataRow=analyticalLinearInterpolator(previous_row, current_row , abscissaToInsertValue, previous_row_abscissa, current_row_abscissa )
                #Step 3.b.iiv
                dataList.insert(rowCounter+rowsAdded, interpolatedDataRow )
                
                #Insertes the values in the abscissa list in the order determined by the sorted abscissa to insert
                abscissaList.insert(rowCounter+rowsAdded, abscissaToInsertValue)
                
                rowsAdded+=1
        
        #This section basically repeats the same process as above, but for inserting rows below if the current_row is not the last row.        
        if rowCounter !=len(abscissa)-1:
            
            next_row=dataList[rowCounter+rowsAdded+1]
            next_row_abscissa=abscissaList[rowCounter+rowsAdded+1]
            
            #Step 3.b.i
            #Loops through the values in the nearZeroIndicesRow to determine if any zeros/ are present in the row
            for columnCounter, nearZeroIndexValue in enumerate (nearZeroIndicesRow):
            
                absDifferenceBottom=abs(current_row[columnCounter]-next_row[columnCounter])
                
                #If a zero/nearZero is present AND the difference between the current row and the above row is significant(greater than the IgnorableDeltaYThreshold)
                #a new abscissa will be calculated based on the halfMinThreshold to be inserted.
                if nearZeroIndexValue==0 and absDifferenceBottom>IgnorableDeltaYThreshold:
                    
                    #Makes the sign of the halfMinThreshold the same as the one above.
                    halfMinThresholdBelow=(IgnorableDeltaYThreshold/2)*signOfElement(next_row[columnCounter])
                    #Interpolates the new abscissa value
                    interpolatedAbscissaValue=analyticalLinearInterpolator(current_row_abscissa,next_row_abscissa, halfMinThresholdBelow, current_row[columnCounter], next_row[columnCounter])
                    abscissaToInsertBelow.append(interpolatedAbscissaValue)
            
            #Step 3.b.ii
            #Sorts the abscissa values so they are isnerted in the right order (increasing) and the rest of the row can be interpolated based off of the abscissa value.
            abscissaToInsertBelow=sorted(abscissaToInsertBelow)
        
            #Iterate through the interpolated abscissa values that need to be inserted
            for counter, abscissaToInsertValue in enumerate (abscissaToInsertBelow):
                
                #Step 3.b.iii
                #Inserts the new row in the dataList containing the interpolated values that correspond to the value in the abscissa.
                interpolatedDataRow=analyticalLinearInterpolator(next_row, current_row , abscissaToInsertValue, next_row_abscissa, current_row_abscissa )
                #Step 3.b.iiv
                dataList.insert(rowCounter+1+rowsAdded, interpolatedDataRow )
                
                #Insertes the values in the abscissa list in the order determined by the sorted abscissa to insert
                abscissaList.insert(rowCounter+1+rowsAdded, abscissaToInsertValue)
                
                rowsAdded+=1

    #Makes the lists back into arrays
    abscissa=numpy.array(abscissaList)
    data=numpy.array(dataList)
    
    #Step 3.c
    #Remakes the zeroIndicesArray and nearZeroIndicesArray so they include the rows inserted
    zeroIndicesArray=(data!=0.0)/1.0
    
    #Step 3.d
    #sets all the zeros in the data set equal to 1/10th of the threshold since the interpolater cannot handle zeros in its row ratio.                
    data[data==0]=IgnorableDeltaYThreshold/10

    return data,abscissa,zeroIndicesArray

#Commented out because no longer in use.
##Returns true when a sign change occurs between two values
#def signChangePresent(value1,value2):
#    signChangePresent=False
#    if value1*value2<0:
#        signChangePresent=True
#        return signChangePresent

#AdujustForSignChange inserts a new row when a sign change is present and the threshold value is crossed. This row contains a zero at the location of the sign change.
#The rest of the values in the row are linearly interpolated.
def adjustForSignChange(data, abscissa, IgnorableDeltaYThreshold):
		
    dataList=list(data)
    abscissaList=list(abscissa)
    
    #Since for loops do not revaluate the len(dataList) everytime it iterates, a number of rows added has to be included that ensures the commands contianed within
    #the loop only operate on the original rows in the dataList.
    rowsAdded=0
    
    #Step 2.a
    #This for loop iterates throught the rows of the original dataList and adds a row below is there is a sign change occuring between any Y-Columns in the X-Row.
    #This row contains a zero in the location of the sign change and the rest of the values are linearly interpolated.
    for rowCounter in range(len(dataList)-1):
        
        #Keeps track of the current_row and abscissa for the row that is being evaluated
        current_row=dataList[rowCounter+rowsAdded]
        current_row_abscissa=abscissaList[rowCounter+rowsAdded]
        
        #Lists of the interpolated abscissa values that need to be inserted above or below a row based on the inserted 0.
        abscissaToInsert=[]
        
        
        next_row=dataList[rowCounter+rowsAdded+1]
        next_row_abscissa=abscissaList[rowCounter+rowsAdded+1]
        
        #Step 2.a.i
        #The below creates vectors of boolean variables. The signChangeVector will return true at a location of a sign chnage across the row. The 
        #absDifferenceSignificant vector will return true at the location of any differences across the row that are greater than the IgnorableDeltaYThreshold.
        signChangeVector=(current_row*next_row<0)
        absDifferenceVector=abs(current_row-next_row)
        absDifferenceSignificant=absDifferenceVector>IgnorableDeltaYThreshold
        
        if any(signChangeVector) and any(absDifferenceSignificant):
            #Loops through the current row to determine if a sign change is present and the difference larger than the threshold
            for columnCounter, current_rowValue in enumerate (current_row):
                
                #Step 2.a.ii
                #If a sign change is present AND the difference between the current row and the above row is significant(greater than the IgnorableDeltaYThreshold)
                #a new abscissa will be calculated based on the 0 to be inserted.
                
                 if signChangeVector[columnCounter] and absDifferenceSignificant[columnCounter]:   
                    #Interpolates the new abscissa value
                    interpolatedAbscissaValue=analyticalLinearInterpolator(current_row_abscissa,next_row_abscissa, 0, current_rowValue, next_row[columnCounter])
                    abscissaToInsert.append(interpolatedAbscissaValue)
            
            #Sorts the abscissa values so they are isnerted in the right order (increasing) and the rest of the row can be interpolated based off of the abscissa value.
            abscissaToInsert=sorted(abscissaToInsert)
            
            #Step 2.a.iii
            #Iterate through the interpolated abscissa values that need to be inserted
            for counter, abscissaToInsertValue in enumerate (abscissaToInsert):
                
                #Inserts the new row in the dataList containing the interpolated values that correspond to the value in the abscissa.
                interpolatedDataRow=analyticalLinearInterpolator(next_row, current_row , abscissaToInsertValue, next_row_abscissa, current_row_abscissa )
                dataList.insert(rowCounter+1+rowsAdded, interpolatedDataRow )
                
                #Insertes the values in the abscissa list in the order determined by the sorted abscissa to insert
                abscissaList.insert(rowCounter+1+rowsAdded, abscissaToInsertValue)
                
                rowsAdded+=1

    #Makes the lists back into arrays
    abscissa=numpy.array(abscissaList)
    data=numpy.array(dataList)
      
    #The analyticalLinearInterpolator will make values extremely close to zero, but not exactly zero, so this changes them to exactly zero.
    data[abs(data)<=(10**-12)]=0
               
    
    return data, abscissa

#The main interpolator inserts abscissa values for both the min and max log ratios. This creates the potential to insert superfluous interpolated rows.
#NOTE: rows refer to the YYY data at each abscissa value.
#removeSuperfluousInterpolatedRows removes the rows corresponding with the abscissa value if the current row and the row after the next row are within 
#the MaxAllowedDeltaYRatio and the IgnorableDeltaYThreshold. 
def superfluousInterpolatedRowCleanup(data, abscissa, insertion_indices, zeroIndicesArray, MaxAllowedDeltaYRatio, IgnorableDeltaYThreshold):
    
    #Make all the zeros in the YYY data equal to 1/10th of the ignorableDeltaYThreshold because a ratio cannot be taken with a zero in
    #any of the rows. This will be set back to zero at the end with the use of the zeroIndicesArray
    data[data==0]=IgnorableDeltaYThreshold/10
    
    #Make the data, abscissa, and zeroIndicesArray so rows can be effeciently deleted
    dataList=list(data)
    abscissaList=list(abscissa)
    zeroIndicesList=list(zeroIndicesArray)
    
    #Initialize rowsChecked variable that keeps track of how many rows were checked, but did not result in deletion
    rowsChecked=-1
    
    #Loops through the insertion indices so original data rows are not deleted
    for insertion_index in insertion_indices[:-2]:
        
        #Initializes the rows that need to be used to decide which rows to delete
        current_row=dataList[insertion_index+rowsChecked]
        next_next_row=dataList[insertion_index+rowsChecked+2]
        
        #Take the ratio of each element in the row in regards to the nex_next_row. This will determine if the two rows fall
        #within the required ratio without need for a row between them.
        row_ratio=next_next_row/current_row
        
        #Creates a boolean array. If an element ratio is between the MaxAllowedDeltaYRatio and the reciprocal of the MaxAllowedDeltaYRatio,
        #True is is put in the location of the element.
        ratioWithinLimits=[(ratio<=MaxAllowedDeltaYRatio and ratio>=MaxAllowedDeltaYRatio**-1.0) for
                           ratio in row_ratio]
        
        #Take the difference between the current_row and the next_next_row
        rowDifference=numpy.abs(current_row-next_next_row)
        
        #Creates a booolean array with True if the difference is significant(larger than the ignorableDeltaYThreshold)
        significantDifference=[(difference>IgnorableDeltaYThreshold) for difference in rowDifference]
        
  
        #If the current_row and the next_next_row are within the ratio and threshold limits and there is no peak present between the rows, delete the next_row from the 
        #data, abscissa, and zeroIndicies
        if all(ratioWithinLimits) and all(significantDifference):
            del dataList[insertion_index+rowsChecked+1]
            del abscissaList[insertion_index+rowsChecked+1]
            del zeroIndicesList[insertion_index+rowsChecked+1]
        else:
            #Since the row was checked, but no row deleted, add to the rowsChecked variable.
            rowsChecked+=1
    
    #Return the data. abscissa, and zeroIndicesArray to lists        
    data=numpy.array(dataList)
    abscissa=numpy.array(abscissaList)
    zeroIndicesArray=numpy.array(zeroIndicesList)
    
    #Remake the originally zero values in the data.
    data=data*zeroIndicesArray
    
    return data,abscissa

# This function is called once an unacceptable gap has been found.
# the limits of that gap are lower_point and upper_point, their
# corresponding abscissa points are the other parameters.
# The function performs the interpolation using Charles' algorithm for
# the data. Then it interpolates that to find the appropriate new abscissa points.
def form_new_abscissa(lower_point, upper_point,
                      abscissa_next_point, abscissa_current_point,
                      log_ratio, factor):

    # Now use Charles' algorithm to find the number of points needed
    # to fill the gap sufficiently.
    ## According to Charles if the largest ratio is exaclty an integer
    ## we should subtract one and it will work perfectly,
    ## else just round down to the nearest integer.
    if log_ratio.is_integer():
        number_needed = int(np.abs(log_ratio) - 1.0)
    else:
        number_needed = int(math.floor(np.abs(log_ratio)))
            
    # Construct a 1-d list that fills the gap
    # Include the data points from the lower and upper rows
    filled_column_gap = [lower_point]

    for new_point in range(number_needed):
        filled_column_gap.append(filled_column_gap[new_point]*factor)

    filled_column_gap.append(upper_point)
            
    # Now interpolate this data to find the corresponding points
    # for modified abscissa.
    ## NOTE: For some reason numpy interpolate can only function
    ## with increasing first argument, hence the if statements
    if lower_point < upper_point:
        filled_abscissa_gap = np.interp(filled_column_gap,
                                        [lower_point,upper_point],
                                        [abscissa_next_point,
                                         abscissa_current_point])
    elif lower_point >  upper_point:
        filled_abscissa_gap = np.interp(list(reversed(filled_column_gap)),
                                        [upper_point,lower_point],
                                        [abscissa_current_point,
                                         abscissa_next_point])
        filled_abscissa_gap = list(reversed(filled_abscissa_gap))

    # return only the points inside the gap not the bounding points
    # which we know from the original data
    return list(filled_abscissa_gap[1:-1])
        

##For Interpolater to work properly, the YYY data must be in numpy array form
#NOTE: rows refer to the YYY data corresponding to each abscissa value.
def Interpolater_new(data,abscissa,MaxAllowedDeltaYRatio=2.0, IgnorableDeltaYThreshold=0.0001, cleanupSuperfluousInterpolatedRows=True):
    
    # First probably should make sure that the abscissa and the data
    # have the same number of rows (abscissa should index the rows of data)
    if (not data.shape[0] == abscissa.shape[0]):
        raise ValueError("data and abscissa dimensions do not match" +
                         " in the interpolater function")
    
    # This set of statments runs the data through all of the necessary preprocessing
    # including casting all of the data to floats, removing sign changes with differences
    # greater than the IgnorableDeltaYThreshold, and removing any zeros from the data set.
    #Step 1: Cast data and abscissa to floats
    data=data.astype(float)
    abscissa=abscissa.astype(float)
    
    #Step 2: Remove Sign changes from Data
    signChangeData=adjustForSignChange(data,abscissa,IgnorableDeltaYThreshold)
    data,abscissa=signChangeData[0],signChangeData[1]
#    print(data,abscissa)
    
    #Step 3: Remove the Effect of zero/nearZero values through row insertion
    zeroAdjustment=adjustForNearZerosAndZeros(data,abscissa,IgnorableDeltaYThreshold)
    data,abscissa,zeroIndicesArray=zeroAdjustment[0],zeroAdjustment[1],zeroAdjustment[2]
#    print(data,abscissa)

    #Step 4: Go through the actual data interpolation

    #Step 4.a
    #Initialize list to be added to in the main interpolation
    # Array to eventually be returned with interpolated data
    interpolated_data = data[0,:]
    # list to eventually be the interpolated abscissa
    interpolated_abscissa = [abscissa[0]]
    # list to hold the filled gap chunks of data for
    # one final insert a the end
    data_gaps = []
    # insertion indices will store the row index
    # where the filled data should be inserted to the main data
    # array
    insertion_indices = []
    
    #Step 4.b.
    dataList=list(data)
    abscissaList=list(abscissa)

    #Step 4.c.
    # loop through abscissa, at each step check the data
    # if any gaps are too large determine the data that fills
    # the worst gap sufficiently
    # then interpolate those new data points to the 'abscissa_values_for_gap_combined'.
    # Finally use abscissa_values_for_gap_combined to interpolate to all of the data columns
    # Note: skip iteration of last point so we dont overrun abscissa
    for point_index,point in enumerate(abscissa[0:-1]):

        # get the relevant rows
        current_row = dataList[point_index]
        next_row = dataList[point_index + 1]
        current_zeros_row=zeroIndicesArray[point_index]
        next_zeros_row=zeroIndicesArray[point_index+1]
        
        # abscissa points
        abscissa_current = abscissaList[point_index]
        abscissa_next = abscissaList[point_index+1]

        #Step 4.c.i.
        # get the absolute difference
        # (by absolute I mean both absolute value and
        # arithmetic difference in contrast to the ratios below)
        row_absolute_difference = np.abs(next_row - current_row)

        # Determine which gaps exceed the absolute threshold
        # (likely most will/should)

        absolute_gaps = [(significantDifference > IgnorableDeltaYThreshold) for
                             significantDifference in row_absolute_difference]
        
        
        #Step 4.c.i.1
        #Create a row with the locations of any of differences between the 
        #next_row and the current_row that are less than the IgnorableDeltaYThreshold
        negligibleAbsoluteDifferencesArray=[]
        for significantDifference in absolute_gaps:
            if not significantDifference:
                negligibleAbsoluteDifferencesArray.append(0)
            else:
                negligibleAbsoluteDifferencesArray.append(1)
        
        #Step 4.c.ii
        #Take the ratio of the values in the next_row to the corresponding values
        #in the current_row
        row_ratio = next_row/current_row
        
        #This makes the row_ratio of sign changes within the threshold equal to
        #1 so the remainder of the function does not try to interpolate between
        #points having a change in sign.
        row_ratio[row_ratio<0]=1
        
        #Step 4.c.ii.1
        # take the log of the ratio in the desired base
        # Note: I dont think that there is a vectorized np log_n function
        ratio_logarithm = np.array(
            [math.log(data_point, MaxAllowedDeltaYRatio)
             for data_point in row_ratio]
            )
        
        #Step 4.c.iii
        #Multiplies the ratio_logarrithm by the current_zeros_row and the 
        #nearZeroDifference to prevent the interpolator from interpolating 
        #between points that are in the ignorableDeltaYThreshold
        ratio_logarithm=ratio_logarithm*current_zeros_row*negligibleAbsoluteDifferencesArray

        
        # is there any abs(log(ratio)) that exceed 1.0?
        # that would indicate that the gap in that column is too large
        # for the relative threshold
        relative_gaps = [(np.abs(ratio) > 1.0) for ratio in ratio_logarithm]

        # combine the relative and absolute gap
        gaps = [(sig_rel_gap_present and sig_abs_gap_present) for (sig_rel_gap_present,sig_abs_gap_present) in
                    zip(relative_gaps, absolute_gaps)]
        
        #Step 4.c.iii.1
        #returns the zero values in the current and next rows back to zero. Also alters the value in the current row in the actual dataList.
        current_row=current_row*current_zeros_row
        dataList[point_index]=dataList[point_index]*current_zeros_row
        next_row=next_row*next_zeros_row
        
        #Step 4.c.i.2 
        #If there are problem gaps that exceed the absolute threshold
        # and the relative threshold deal with them (i.e. add the points
        # to the modified_abscissa that we need)
        if any(gaps):
            # So now there is a gap (i.e. large gap in at least in one column).
            # If there is a problem gap in both directions (e.g. gap1 == [2,12]
            # and gap2 == [7,1.5]) then we will need to perform the interpolation
            # and backwards interpolation to form the modified_abscissa for both
            # gaps, then just combine the two in an ordered fashion and append
            # like normal
            
            #Step 4.c.iv
            # keep only the column indices (in this pair of rows) which had a gap
            gap_indices = [index for (index,gap)
                                     in enumerate(gaps) if gap]

            gap_log_ratios = [ratio_logarithm[gap] for
                                 gap in gap_indices]

            # Find the maximum and miniumum of the gap log(ratio)
            # dont care if they are unique, if there are two identical extrema
            # they will require identical interpolation anyway
            max_log_ratio = max(gap_log_ratios)
            min_log_ratio = min(gap_log_ratios)

            # Doesn't need to be unique
            max_log_ratio_index = gap_indices[gap_log_ratios.index(max_log_ratio)]
            min_log_ratio_index = gap_indices[gap_log_ratios.index(min_log_ratio)]

            abscissa_values_for_gap_from_rising_data = []
            abscissa_values_for_gap_from_falling_data = []
            
            #Step 4.c.v
            # if there is an increasing gap find the points to add to
            # modified abscissa
            if max_log_ratio > 0:
                abscissa_values_for_gap_from_rising_data = form_new_abscissa(current_row[max_log_ratio_index],
                                                   next_row[max_log_ratio_index],
                                                   abscissa[point_index],
                                                   abscissa[point_index + 1],
                                                   max_log_ratio,
                                                   MaxAllowedDeltaYRatio)
            # if there is a decreasing gap find the points to add to
            # modified abscissa
            if min_log_ratio < 0:
                abscissa_values_for_gap_from_falling_data = form_new_abscissa(current_row[min_log_ratio_index],
                                                   next_row[min_log_ratio_index],
                                                   abscissa[point_index],
                                                   abscissa[point_index + 1],
                                                   min_log_ratio,
                                                   MaxAllowedDeltaYRatio**-1.0)

            #Step 4.c.v.1
            # Combine the two if there is both and remove any duplicates
            abscissa_values_for_gap_combined = abscissa_values_for_gap_from_rising_data + abscissa_values_for_gap_from_falling_data
            abscissa_values_for_gap_combined = list(set(abscissa_values_for_gap_combined))
            
           
            # the abscissa needs to be ordered and increasing or decreasing
            # depending on how whether the boundary abscissa points are
            # increasing or decreasing.
            if abscissa[point_index] < abscissa[point_index + 1]:
                abscissa_values_for_gap_combined = sorted(abscissa_values_for_gap_combined)

            if abscissa[point_index] > abscissa[point_index + 1]:
                abscissa_values_for_gap_combined = list(reversed(sorted(abscissa_values_for_gap_combined)))

            
            #Step 4.c.v.2       
            
            #The following code vectorizes the linear row interpolation based on the abscissa values to be inserted.
            #The slopesVector hold the slopes for each column of the YYY data
            slopesVector=(next_row-current_row)/(abscissa_next-abscissa_current)
            #interceptVector hold the linear intercepts for each column of the YYY data
            interceptVector=current_row-slopesVector*abscissa_current
            
            #The duplicatedAbscissaArray is created to easily perform calculations for the entire list of abscissa values to be inserted at once.
            #There is a numpy array row for each YYY column and a column for each abscissa value to be inserted. The each column in the 
            #duplicatedAbscissaArray will have the same value in each of the rows. This value is the abscissa value to be inserted and the the 
            #columns of the YYY data to be interpolated to
            #Thes
            duplicatedAbscissaArray=numpy.ones((len(slopesVector),len(abscissa_values_for_gap_combined)))*abscissa_values_for_gap_combined
            
            
            #The interpolated rows of the YYY data are calculated with a y=mx+b equation using the slopes and intercepts found above. To calculate
            #correctly, the duplicatedAbscissaArray must be transposed.
            filled_data=numpy.multiply(numpy.transpose(duplicatedAbscissaArray),slopesVector)+interceptVector
            
            #Step 4.c.v.3
            # Add the filled data and abscissa points to the
            # interpolated data and abscissa (we handle the
            # boundary points elsewhere)
            #data_gaps.append(filled_data.reshape(len(abscissa_values_for_gap_combined),len(current_row)))
            for row in filled_data[:,]:
                data_gaps.append(row)
                insertion_indices.append(point_index+1)
            #interpolated_data = np.vstack((interpolated_data, filled_data))
            interpolated_abscissa.extend(list(abscissa_values_for_gap_combined))
            
        # END if gap block

        # Regardless of the existence of a gap in this row
        # make sure to add the upper row to the interpolated_data
        # and the upper point to abscissa
        interpolated_data = np.vstack((interpolated_data, next_row))
        interpolated_abscissa.extend([abscissa_next])
        
        
    #Else if there are x-rows that need to be inserted, insert the rows and remake the zeroIndicesArray before
    #returning the data
    if len(insertion_indices)!=0:
        #Step 4.d.
        #For the proper data_gaps to be inserted into the proper insertion points in the list,
        #the for loop has to be iterated backwards. This must have to do with the order that the
        #insertion_indices and data_gaps are filled.
        for counter in range(len(insertion_indices)-1,-1,-1):
            dataList.insert(insertion_indices[counter],data_gaps[counter])
            
        #The for loop does not return any originally zero values in the last row of the dataList back to zero.
        #This statement multiplies the last row in the dataList by the last row in the zeroIndiciesArray to return
        #any originally zero values in this row back to zero.
        dataList[-1]=dataList[-1]*zeroIndicesArray[-1]
        
        #return the lists to arrays
        data=numpy.array(dataList)
        abscissa=numpy.array(interpolated_abscissa)
        
        #Step 4.e
        #create a fresh zeroIndicesArray contianing the the zero locations with the inserted rows
        zeroIndicesArray=(data!=0.0)/1.0

        #Step 5
        
        #If the user wants to remove the superfluous rows created from the main interpolator, the cleanupSuperfluousInterpolatedRows will be defined as
        #True. The default value for this is set to True.
        #NOTE: rows refer to the YYY data corresponding to each abscissa value.
        if cleanupSuperfluousInterpolatedRows:
            
            dataWithsuperfluousRowsRemoved=superfluousInterpolatedRowCleanup(data, abscissa, insertion_indices, zeroIndicesArray, MaxAllowedDeltaYRatio, IgnorableDeltaYThreshold)
            data=dataWithsuperfluousRowsRemoved[0]
            abscissa=dataWithsuperfluousRowsRemoved[1]
    

    #Else no rows needed to be interpolated above, and there would be no rows that need to be inserted so the 
    #data does not need to be altered and can be returned.
    else:
        data = numpy.multiply(data,zeroIndicesArray)


    #Step 6
    return data, abscissa, insertion_indices
