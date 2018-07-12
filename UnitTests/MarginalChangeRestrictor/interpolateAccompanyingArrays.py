# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 11:45:11 2018

@author: Andrea
"""
import numpy 

def analyticalLinearInterpolator (ordinateValues1, ordinateValues2, abscissaValueToInterpolateTo, abscissaValue1, abscissaValue2):
        slopesVector=(ordinateValues2-ordinateValues1)/(abscissaValue2-abscissaValue1)
        interceptVector=ordinateValues1-slopesVector*abscissaValue1
        interpolatedOrdinateValues=slopesVector*abscissaValueToInterpolateTo+interceptVector
        return interpolatedOrdinateValues

#The dataRangeSpecifierMarginalChangeRestrictor is a preprocessing step that changes the abscissa of the accompanying concentration bounds file to 
#match that of the interpolated data. When the times do not match that of the marginalChangeRestricted_abscissa, it will add a new row full of 
#linearInterpolations across the YYY columns of the bounds and other data contained in the accompanyingArray file.

#When used, the marginalChangeRestricted_abscissa must be ordered. The abscissa in the accompanyingArray must match that of the actual data 
#pre-marginalChangeRestrictor. The first and last number in the csv file must match that of the abscissa. The accompanying array cannot contain a 
#header.

def interpolateAccompanyingArrays(marginalChangeRestricted_abscissa, accompanyingArray):
    
    
    #The dataRangeSpecifierInterpolator will progress if points were added in the interpolator. The easiest way to do this is to
    #compare the lengths of the interpolated abscissa to that of the datafromcsv
    if len(marginalChangeRestricted_abscissa)!=len(accompanyingArray):
        accompanyingList=list(accompanyingArray)

        #Loops through the marginalChangeRestricted_abscissa.
        for rowCounter, marginalChangeRestricted_abscissa_value in enumerate(marginalChangeRestricted_abscissa[:-1]):
            
            #Check if the marginalChangeRestricted_abscissa_value matches the abscissa of the same row in the datafromcsvfile
            if marginalChangeRestricted_abscissa_value !=accompanyingList[rowCounter][0]:
                
                #Set the previous and current row to be used to interpolate between them based on the marginalChangeRestricted_abscissa_value that did not
                #match the value in the file
                previous_row=numpy.array(accompanyingList[rowCounter-1])
                current_row=numpy.array(accompanyingList[rowCounter])
                
                #Use the analyticalLinearInterpolator to interpolate 
                interpolatedExpandedRow=analyticalLinearInterpolator(previous_row,current_row,marginalChangeRestricted_abscissa_value,previous_row[0],current_row[0])
                accompanyingList.insert(rowCounter,interpolatedExpandedRow)

        
        expandedArray=numpy.array(accompanyingList)        
                  
    return expandedArray
