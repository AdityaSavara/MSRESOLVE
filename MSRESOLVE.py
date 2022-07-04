import bisect
import copy 
import ParsingFunctions as parse
import numpy
import csv
import time
import timeit
import math
import sys
import pandas
import XYYYDataFunctionsSG as DataFunctions
import os
import shutil
import importlib
from numpy import genfromtxt
import export_import as ei
#G stands for Global, and is used to draw data from the UserInput File, and to store data during processing.
import UserInput as G; importlib.reload(G) #import the user input and reload the module to get rid of any unwanted variables in the namespace

############################################################################################################################################
#########################################################Best Mass Fragment Chooser#########################################################
############################################################################################################################################
#This function is made to conduct prelimiary checks. The checks are designed to
#fail only if there is not chance of the the chosen reference data passing the
#SLS method. There are 2 cases where this occurs :1) a molecule does not 
#contain reference intensities for any mass fragments in the combination or 2)
#The reference data for the mass fragments does not contain any zeros.
def passesRowsSumChecks(rowSumsList, massFragCombination, allOverlappingPatterns, numberOfMassFragments=None, numberOfRows=None):
    if numberOfMassFragments==None:
        numberOfMassFragments=len(massFragCombination)##It is slightly faster to pass in the length variable to make it shorter      
    passesRowsSumChecks=True#Initialize return variable as true
    if 0 in rowSumsList: #Check if any row is entirely full of zeros
        passesRowsSumChecks=False
    elif numpy.array(rowSumsList).all() == numberOfMassFragments: #Check if the array passed to it is full. If it is full, it appends to the allOverlapping patterns list
        allOverlappingPatterns.append(massFragCombination)###May be an error with the use of allOverlapping patterns
        #FIXME: #TODO: should be printed to a file.
        #print(allOverlappingPatterns)
        passesRowsSumChecks=False

    return passesRowsSumChecks, allOverlappingPatterns #Return true if the rowsSumsList passes the check

#The function maintains two lists: 1)that contains the objective function values that 
#need to be kept in a sorted order 2) a parallel list where a value needs to be inserted 
#according to the sorting of the first one. It also takes in an integer value, N, that limits 
#the the length of the two lists to N values. Using a binary search method, it will 
#find the location where the value to insert will be inserted(if possible). The
#value will be inserted there and the last value removed from the list (if
#applicable).  The bisect method used requires the list be presorted in ascending order.
#Thus, the algorithm here inserts values in a way to create an ascending ordered list.

#By default storeAndPop is used to keep the best N values based on minimizing
#the objective function values. If instead it is desirable to retain values 
#with the objective function maximized, the optional argument 'optimumType'
#should be set to `Maximum`.

#Alternatively, the objective function values could be multiplied by -1.
#The function supports multidimensional objective functions in nested objects
#such as tuples or lists, e.g., (A,B,C) in which case it will be sorted
#by A, then B, then C.  For multidimensional objective functions, it is
#necessary to have the dimensions as either all "maximum" or all "minimum" type.
#Using a -1 factor can be helpful in this regard (e.g., passing in (A,-B,C) etc.

#If the values in the sample space for parallelList are not unique it is 
#possible that this repeated calls to this function could lead to 
#a parallelList of a particular value repeated many times. If repeated values
#are undesired, then excludeDuplicates can be set to False.
def storeAndPop(objectiveFunctionValuesList, objectiveFunctionValueToInsert, 
                parallelList, valueToInsertInParallelList, maxItemsAllowed,
                optimumType="Minimum", excludeDuplicates=True):
    
    #Find the insertion index where the value will be inserted by using a binary
    #search
    insertionIndex=bisect.bisect(objectiveFunctionValuesList,
                                 objectiveFunctionValueToInsert)

    #Initialize a variable to keep track of if a value was inserted into the
    #list.
    valueStoredInList=False
    
    #If it is a duplicate exit now without checking everything else
    #Note that we only check the value in parallel list to the left of
    #the insertion index. This is because bisect.bisect() will specify 
    #and insertion index such that a duplicate is inserted
    #to the right of the original.
    if (excludeDuplicates and len(parallelList) and 
        (parallelList[insertionIndex-1] == valueToInsertInParallelList)):
        
        #This value is a duplicate. Return the original lists.
        return objectiveFunctionValuesList, parallelList, valueStoredInList    
    
    #If the list isn't yet filled, the value will inherently be in the top N
    #value in the list. This value can just be inserted at the insertionIndex.
    if len(objectiveFunctionValuesList)<maxItemsAllowed:    
        objectiveFunctionValuesList.insert(insertionIndex,
                                           objectiveFunctionValueToInsert)
        parallelList.insert(insertionIndex, valueToInsertInParallelList)
        valueStoredInList=True
    #If the list already contains N elements, a new element could either be 
    #inserted in the list or at the end of the list. For optimumType == Minimum,
    #an item at the end would be worse and thus nothing should be added. 
    #However, for optimumType == Maximum an object at the end would be the best
    #and thus should be added.
    elif (len(objectiveFunctionValuesList) == maxItemsAllowed and 
            (insertionIndex<maxItemsAllowed or optimumType=="Maximum")):
        #insert the value to insert in the location found through the binary
        #search
        objectiveFunctionValuesList.insert(insertionIndex,
                                           objectiveFunctionValueToInsert)
        parallelList.insert(insertionIndex, valueToInsertInParallelList)
        valueStoredInList=True
        
        if optimumType == 'Minimum':
            #delete the last element since something was added to the list
            del objectiveFunctionValuesList[-1]
            del parallelList[-1]
        elif optimumType == 'Maximum':
            #delete the first element since something was added to the list
            del objectiveFunctionValuesList[0]
            del parallelList[0]   
        else:
            raise  ValueError("optimumType must be either 'Maximum' " +
                              "or 'Minimum'")

    return objectiveFunctionValuesList, parallelList, valueStoredInList

#The rough uniqueness check is a limiting check that takes the mass fragment
#combinations that pass the row sums check and builds a lists of 
#keep_N_ValuesInRoughUniquenessCheck that contain the largest number of zeros
#and the corresponding mass fragment combination for each of the values in the
#first list.
#The largest number of zeros would be most likely to pass the SLS method.
#It calculates a sum that roughly expresses how unique the molecular mass 
#fragments are to the different molecules, but this is a quick and not-rigrous 
#method. Then, the value is stored *only* if it is in the top N of the values 
#so far.
def roughUniquenessCheck(rowSumsList, smallestRowsSumsList,topMassFragCombinationsList, keep_N_ValuesInRoughUniquenessCheck, massFragCombination):

    #We want to save the smallest sum since that would contain the smallest 
    #number of zeros.
    roughUniqueness=numpy.sum(rowSumsList) 

    #Use Store and Pop to add the tuple to the list of top rough uniquness 
    #combinations. This will only save a user specified number of tuples.
    [smallestRowsSumsList,topMassFragCombinationsList,valueStoredInRUTopList]=storeAndPop(smallestRowsSumsList, roughUniqueness,topMassFragCombinationsList, massFragCombination, keep_N_ValuesInRoughUniquenessCheck)
    return smallestRowsSumsList,topMassFragCombinationsList, valueStoredInRUTopList

#The significance factor check is a limiting check the selects the mass 
#fragment combinations having the largest sum of significance factors. It
#calculates the significance factors for each element in the sub-reference 
#array (refernce array with zeros for all the data for mass fragments that 
#aren't needed). Is then takes the sum of all of the significance factors.
#The top N mass fragment combinations are then stored in a list in decending 
#order according to the magnitude of their significance sum. 
#Basically, it calculates the significance factor for each element in the
#chosen reference array and sums all of the significane values for the whole 
#array. It keeps the mass fragments that have the largest magnitude of 
#significance sum.
def significanceFactorCheck(chosenReferenceIntensities,largestMagnitudeSigFactorSumsList,topMassFragCombinationsList, massFragCombination, keep_N_ValuesInSignificanceFactorCheck, moleculesLikelihood):
    
    #Initialize the sum of the significance factors
    sigFactorSum=0
    
    #The elemSignificanceCalculator finds the significance factors a column at 
    #a time. This loops through the columns(molecules)
    for columnCounter in range(len(chosenReferenceIntensities[0])):
        
        #Finds the significance factors for each element in the column. This 
        #does not include the first column since that is the mass fragment 
        #numbers
        significanceColumnDataList=ElemSignificanceCalculator(chosenReferenceIntensities, columnCounter, moleculesLikelihood)
        
	#Sums the significance factors across the column and to the
        #sum for the whole ref data array. The larger in magnitude this is, the 'better'.
        sigFactorSum+=sum(significanceColumnDataList)
        
        ####Currently there is no real need to maintain a significance data 
        ####list for the whole array
    
    #The subtraction is used to make the sum negative. The binary search used 
    #will order from lowest to highest. A larger magnitude for this value then means
    #the most negative which will be kept during the store and pop function.
    #Note that we use the wording of largest magnitude so that our phrasing remains the same
    #when talking about the negative of the sum.
    negativeOfSigFactorSum=-1*sigFactorSum    
    
    #Uses store and pop to maintian a list of the mass fragment with the
    #largest significance factors.
    #The below line only keeps the combinations with the largest magnitude (most negative) of the negativeOfSigFactorSums.
    [largestMagnitudeSigFactorSumsList,topMassFragCombinationsList, valueStoredInSFTopList]=storeAndPop(largestMagnitudeSigFactorSumsList,negativeOfSigFactorSum,topMassFragCombinationsList, massFragCombination, keep_N_ValuesInSignificanceFactorCheck)
    return largestMagnitudeSigFactorSumsList, topMassFragCombinationsList, valueStoredInSFTopList
############################################################################################################################################
################################################Algorithm Part 1: Pre-Processing the data###################################################
############################################################################################################################################
# this function will take in the mass fragments that are getting changed, the slope of the line
#of which they are being changed by and finally the intercept of that line which is being subtracted
#from all of their values, so much of this function is going to be prompting and getting values from
#the user
def SlopeEliminator (ExperimentData,backgroundMassFragment,backgroundSlopes,backgroundIntercepts): 
    for mf_counter in range(len(backgroundMassFragment)):#array-indexed for loop
        for mass_fragment_numbers_counter in range(len(ExperimentData.mass_fragment_numbers)):#array-indexed for loop
            for times_counter in range(len(ExperimentData.times)): #array-indexed for loop
                if ExperimentData.mass_fragment_numbers[mass_fragment_numbers_counter] == backgroundMassFragment[mf_counter]:#if the right mass fragment number is selected then its slope is eliminated
                    subtraction_value =  ExperimentData.times[times_counter] * float(backgroundSlopes[mf_counter]) + float(backgroundIntercepts[mf_counter]) 
                    ExperimentData.workingData[times_counter,mass_fragment_numbers_counter] = ExperimentData.workingData[times_counter, mass_fragment_numbers_counter] - subtraction_value
    #return collected #no return is necessary. this is implied. the function takes a "pointer/reference" to the collected array and then modifies it directly.
        
#this function works by drawing the input values from the input file which includes a baselineType, a time range, and the mass fragment
#that is being altered. the time range is taken, and using a for loop, all the values between those two times (including those two times)
#are used in finding the average or polyfit, the collected signals are also picked up here, in order to use them for finding the polyfit. 
#So then, depending of the baselineType, the mass fragments row will be found using a nested for loop with an if statement, which then
# subtracts what ever was found- the polyfit or the average from each data point in the collected data set. 
def LinearBaselineCorrectorSemiAutomatic(ExperimentData,baselineType,massesToBackgroundCorrect,earlyBaselineTimes,lateBaselineTimes):
    
    #This section of code is a patch to permit easier usage of the Baseline Correction
    # it Applies background correction to all fragments if none are listed,
    # but corrector is turned on
    if len(massesToBackgroundCorrect) == 0:
        massesToBackgroundCorrect = ExperimentData.mass_fragment_numbers
    # Applies timeRange to all fragments if only one time is given
    if (len(earlyBaselineTimes) == 1 and len(massesToBackgroundCorrect) != 1):
        newTimeList = earlyBaselineTimes
        for i in range(len(massesToBackgroundCorrect)-1):
            newTimeList = numpy.vstack((newTimeList,earlyBaselineTimes))
        earlyBaselineTimes = newTimeList
    # Applies timeRange2 to all fragments if only one time is given
    if (len(lateBaselineTimes) == 1 and len(massesToBackgroundCorrect) != 1):
        newTimeList = lateBaselineTimes
        for i in range(len(massesToBackgroundCorrect)-1):
            newTimeList = numpy.vstack((newTimeList,lateBaselineTimes))
        lateBaselineTimes = newTimeList
        
    #the mass_fragment_numbers list has the list of all mass fragments that were measured.
    #massesToBackgroundCorrect is the list of mass fragments to do the correction on (which can be a smaller set of masses).
    selectedtimes = []
    selectedsignals = []
    slopelist = []
    interceptlist = []
    averagelist = []
    starttimelist = []
    endtimelist = []
    starttimelist2 = []
    endtimelist2 = []
    #since the format of the user input ranges change, this loop puts that information back into the form
    #originally used (so that the function will work), with this loop here
    for earlyBaselineIndex in range(len(earlyBaselineTimes)):#array-indexed for loop
        starttimelist.append(earlyBaselineTimes[earlyBaselineIndex][0])
        endtimelist.append(earlyBaselineTimes[earlyBaselineIndex][1])
        if len(lateBaselineTimes) != 0:#if there is only one time range, this part of the loop doesn't happen
            starttimelist2.append(lateBaselineTimes[earlyBaselineIndex][0])
            endtimelist2.append(lateBaselineTimes[earlyBaselineIndex][1])

    #this for loop goes through all the data, getting lists of the times and signals, which are then used
    #to make lists of slopes, intercepts and averages, which can be used to alter the collected data
    for massesToCorrectIndex in range(len(massesToBackgroundCorrect)):#array-indexed for loop
        for measured_masses_counter in range(len(ExperimentData.mass_fragment_numbers)):#array-indexed for loop
            if massesToBackgroundCorrect[massesToCorrectIndex] == ExperimentData.mass_fragment_numbers[measured_masses_counter]:#if the right mass fragment is found, we now have that index
                for timescounter in range(len(ExperimentData.times)):#array-indexed for loop
                    if ExperimentData.times[timescounter] >= starttimelist[massesToCorrectIndex]:# once the loop gets past the start time it starts appending all the times and signals
                        selectedtimes.append(ExperimentData.times[timescounter])
                        selectedsignals.append(ExperimentData.workingData[timescounter,measured_masses_counter])
                    if ExperimentData.times[timescounter] > endtimelist[massesToCorrectIndex]:# once it gets past the end time, it deletes all the values that are being appended
                        selectedtimes.pop()
                        selectedsignals.pop()
                    if len(starttimelist2) != 0:# same as above, it might not exist
                        if ExperimentData.times[timescounter] >= starttimelist2[massesToCorrectIndex]:
                            selectedtimes.append(ExperimentData.times[timescounter])
                            selectedsignals.append(ExperimentData.workingData[timescounter,measured_masses_counter])
                        if ExperimentData.times[timescounter] > endtimelist2[massesToCorrectIndex]:
                            selectedtimes.pop()
                            selectedsignals.pop()
                    if timescounter == len(ExperimentData.times)-1:#once the loop is finished getting all that, the loop saves all the data for each molecules before going to the next molecule                      
                        [slope,intercept] = numpy.polyfit(numpy.array(selectedtimes),numpy.array(selectedsignals),1)
                        average = numpy.average(selectedsignals)
                        slopelist.append(slope)
                        interceptlist.append(intercept)
                        averagelist.append(average)
                        selectedtimes = []
                        selectedsignals = []

        #these are the if statements that choose what happens based on user baselineType in the input data file
    if len(baselineType) == 1:
        baselineTypeholder = []
        for length in range(len(massesToBackgroundCorrect)):
            baselineTypeholder.append(baselineType[0])
        baselineType = baselineTypeholder
    for MassFragmentIndex in range(len(massesToBackgroundCorrect)):#array-indexed for loop
        if baselineType[MassFragmentIndex] == 'flat': #the different baselineTypes subtract different things
            for measured_masses_counter in range(len(ExperimentData.mass_fragment_numbers)):#array-indexed for loop
                if massesToBackgroundCorrect[MassFragmentIndex] == ExperimentData.mass_fragment_numbers[measured_masses_counter]:#when the index for the mass fragment is obtained
                    try:
                        ExperimentData.workingData[:,measured_masses_counter] = ExperimentData.workingData[:,measured_masses_counter] - averagelist[MassFragmentIndex]
                    except IndexError:
                       print("Warning: LinearBaselineCorrectorSemiAutomatic has failed.\
                       It is possible that you have entered an observed mass fragement in massesToBackgroundCorrect that \
                       does not appear in the reference data. If this is the case, \
                       remove that mass fragment and rerun the program.")
                       sys.exit()
        if baselineType[MassFragmentIndex] == 'linear':#other option
            for measured_masses_counter in range(len(ExperimentData.mass_fragment_numbers)):#array-indexed for loop
                if massesToBackgroundCorrect[MassFragmentIndex] == ExperimentData.mass_fragment_numbers[measured_masses_counter]:#same as above
                    try:
                        ExperimentData.workingData[:,measured_masses_counter] = ExperimentData.workingData[:,measured_masses_counter] - (slopelist[MassFragmentIndex]*ExperimentData.times + interceptlist[MassFragmentIndex])
                    except IndexError:
                       print("Warning: LinearBaselineCorrectorSemiAutomatic has failed.\
                       It is possible that you have entered an observed mass fragement in massesToBackgroundCorrect that \
                       does not appear in the reference data. If this is the case, \
                       remove that mass fragment and rerun the program.")
                       sys.exit()
     #no return is necessary. this is implied. the function takes a "pointer/reference" to the collected array and then modifies it directly.
        
        
#this function is going to set a certain threshold based on the data, and eliminate any values below this,
#setting them to 0, this serves to eliminate negatives and insignificant values - this means that the rest of the
#script can solve for the signals more easily and eliminate a good range of possibilities. There are two tests:
#either an absolute threshold, or one that is a percentage of the max signal. The signal should pass both checks.
def LowerBoundThresholdFilter (ExperimentData,massesToLowerBoundThresholdFilter,lowerBoundThresholdPercentage,lowerBoundThresholdAbsolute):
    #test case for whether all masses should be filtered
    #if no masses are listed
    if len(massesToLowerBoundThresholdFilter) == 0:
        #set all mass to be filtered
        massesToLowerBoundThresholdFilter = ExperimentData.mass_fragment_numbers
    #For lower bound absolute value, if somebody has set only one value, then all masses will receive that. Otherwise, only the chosen masses will.
    #We also need to convert this into a float if it's not already a float.
    if len(lowerBoundThresholdAbsolute) == 1:
        ThresholdTemporaryVariable = lowerBoundThresholdAbsolute[0]
        lowerBoundThresholdAbsolute = []#making list blank again so I can append to it as many times as needed.
        for mass in massesToLowerBoundThresholdFilter:
            lowerBoundThresholdAbsolute.append(ThresholdTemporaryVariable) #making float and appending.        
    if len(lowerBoundThresholdAbsolute)!=0:
        #This makes a numpy array at the end.
        lowerboundByAbsoluteArray = numpy.array(lowerBoundThresholdAbsolute)
    
    #if somebody has set only one value for the percentage or the absolute threshold, then all masses will receive that. Otherwise, only the chosen masses will.
    if len(lowerBoundThresholdPercentage) ==1:
        percentageTemporaryVariable = lowerBoundThresholdPercentage[0]
        lowerBoundThresholdPercentage = []#making list blank again so I can append to it as many times as needed.
        for mass in massesToLowerBoundThresholdFilter:
            lowerBoundThresholdPercentage.append(float(percentageTemporaryVariable)) #making float and appending.
    if len(lowerBoundThresholdPercentage)!=0:
        #finding maximum of each column so that I can do the percentage case:
        MaximumIntensitiesOfMassesToLowerBoundThresholdFilter = [ ]
        for measured_masses_counter in range(len(ExperimentData.mass_fragment_numbers)): #array-indexed for loop to loop across all masses
            for mtbc in range(len(massesToLowerBoundThresholdFilter)):#array-indexed for loop that loops acrosses masses to be filtered
                if ExperimentData.mass_fragment_numbers[measured_masses_counter] == massesToLowerBoundThresholdFilter[mtbc]: #This passes *only* for the indices of masses to be filtered.
                    MaximumIntensitiesOfMassesToLowerBoundThresholdFilter.append(max(ExperimentData.workingData[:,measured_masses_counter]))
        #This converts the percentages to absolute values and makes a numpy arry also.
        lowerboundByPercentageArray = numpy.array(lowerBoundThresholdPercentage) * numpy.array(MaximumIntensitiesOfMassesToLowerBoundThresholdFilter)
        
    #Now we will apply the actual filter to each point for each mass that needs to be applied to.
    #We use the two loops again.
    for measured_masses_counter in range(len(ExperimentData.mass_fragment_numbers)): #array-indexed for loop to loop across all masses
        for mtbc in range(len(massesToLowerBoundThresholdFilter)):#array-indexed for loop that loops acrosses masses to be filtered
            if ExperimentData.mass_fragment_numbers[measured_masses_counter] == massesToLowerBoundThresholdFilter[mtbc]: #This passes *only* for the indices of masses to be filtered.
                #Below looks at each point for a given mass, and then applies the appropriate filter.
                for point in enumerate(ExperimentData.workingData[:,measured_masses_counter]): #note that "point" is actually the a tuple with both the value and the index in form (index,value).
                    if len(lowerBoundThresholdAbsolute)!=0: #Check if length of lowerbound threshold list is not equal to zero before applying it
                        if point[1] < lowerboundByAbsoluteArray[mtbc]: #check if the value of the point is lower than the absolute threshold.
                            ExperimentData.workingData[:,measured_masses_counter][point[0]] = 0.0  #this changes the value to 0.0 at index given by point[0] within the array of masses.
                    if len(lowerBoundThresholdPercentage)!=0: #Check if length of lowerbound threshold by percentage list is not equal to zero before applying it
                        if point[1] < lowerboundByPercentageArray[mtbc]: #check if the value of the point is lower than the absolute threshold.
                            ExperimentData.workingData[:,measured_masses_counter][point[0]] = 0.0 #this changes the value to 0.0 at index given by point[0] within the array of masses.
    #return collected #no return is necessary. this is implied. the function takes a "pointer/reference" to the collected array and then modifies it directly.    
    
       
'''
This function simply searchs two arrays and returns an array with any and all 
overlapping numbers
'''
def MatchingNumbers(Array1, Array2):

    return list(set(Array1) & set(Array2))
'''
THis function determines and returns the ABC correction coefficients for Tuning Correction based off 
of a literature (NIST) reference pattern and a measured (Experimental)reference pattern
'''
def ABCDetermination(ReferencePatternExistingTuning_FileNameAndForm, ReferencePatternDesiredTuning_FileNameAndForm, exportCoefficients=True):
    #The two arguments are actually lists containing strings of the name and form type: (["FileNameOne.csv","xyyy"],["FileNameTwo.csv", "xyyy"])
    '''
    Step 1: Read in all neccessary information from fragment patterns
    Step 2: populate all variables and standardize signals by 100 

    '''
    
    if G.minimalReferenceValue !='yes':
        print("Warning: The ABCDetermination will occur without threshold filtering, since that setting is off.")
        
    if G.extractReferencePatternFromDataOption == 'yes':
        print("Warning: MSRESOLVE is set to do a tuning correction and is also set to extract the reference pattern from the data.  The tuning correction will be applied to the extracted pattern.  This is not a typical procedure. Typically, the extract pattern from data feature is off when the tuning corrector is on.")
       
    #For simplicity, we will put the items into temporary items, then into dictionaries that we can then access.
    ReferencePatternExistingTuningDict = {}
    [provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form] = readReferenceFile(*ReferencePatternExistingTuning_FileNameAndForm)
    ReferencePatternExistingTuningDict['molecules']=molecules
    ReferencePatternExistingTuningDict['provided_reference_patterns'] = provided_reference_patterns
    ReferencePatternExistingTuningDict['provided_reference_patterns'] = StandardizeReferencePattern(ReferencePatternExistingTuningDict['provided_reference_patterns'],len(molecules)) #this does have the molecular weight as the first column.
    if G.minimalReferenceValue =='yes':
        ReferencePatternExistingTuningDict['provided_reference_patterns'] = ReferenceThresholdFilter(ReferencePatternExistingTuningDict['provided_reference_patterns'],G.referenceValueThreshold)
    
    ReferencePatternDesiredTuningDict = {}
    [provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form] = readReferenceFile(*ReferencePatternDesiredTuning_FileNameAndForm)
    ReferencePatternDesiredTuningDict['molecules']=molecules
    ReferencePatternDesiredTuningDict['provided_reference_patterns'] = provided_reference_patterns
    ReferencePatternDesiredTuningDict['provided_reference_patterns'] = StandardizeReferencePattern(ReferencePatternDesiredTuningDict['provided_reference_patterns'],len(molecules)) #this does have the molecular weight as the first column.
    if G.minimalReferenceValue =='yes':
        ReferencePatternDesiredTuningDict['provided_reference_patterns'] = ReferenceThresholdFilter(ReferencePatternDesiredTuningDict['provided_reference_patterns'],G.referenceValueThreshold)

    '''
    Step 3a: Truncate to the molecules which match.
    '''
    OverlappingMolecules = numpy.intersect1d(ReferencePatternExistingTuningDict['molecules'],ReferencePatternDesiredTuningDict['molecules'] )
    if len(OverlappingMolecules) == 0:
        return [0,0,1], [[0,0,0],[0,0,0],[0,0,0]]    
    OverlappingFragments = numpy.intersect1d(ReferencePatternExistingTuningDict['provided_reference_patterns'][:,0],ReferencePatternDesiredTuningDict['provided_reference_patterns'][:,0])    
    if len(OverlappingFragments) == 0:
        return [0,0,1], [[0,0,0],[0,0,0],[0,0,0]]                                          
    [OverlappingMolecules,ReferencePatternExistingTuningDict['OverlappingIndices'], ReferencePatternDesiredTuningDict['OverlappingIndices']]  = numpy.intersect1d(ReferencePatternExistingTuningDict['molecules'],ReferencePatternDesiredTuningDict['molecules'], return_indices=True)
    ReferencePatternExistingTuningDict['OverlappingMolecules'] = OverlappingMolecules
    ReferencePatternDesiredTuningDict['OverlappingMolecules'] = OverlappingMolecules
    
    ReferencePatternExistingTuningDict['overlapping_provided_reference_patterns'] = []
    #Now need to get only the overlapping indices values. The numpy take function will work if we do it one row at a time:
    for rowIndex,row in enumerate(ReferencePatternExistingTuningDict['provided_reference_patterns']):
            massFragment = row[0]
            truncatedIntensities = numpy.take(row[1:],ReferencePatternExistingTuningDict['OverlappingIndices'])    
            rowTruncated = numpy.insert(truncatedIntensities, 0, massFragment)
            ReferencePatternExistingTuningDict['overlapping_provided_reference_patterns'].append(rowTruncated)    
    
    ReferencePatternDesiredTuningDict['overlapping_provided_reference_patterns'] = []
    #Now need to get only the overlapping indices values. The numpy take function will work if we do it one row at a time:
    for rowIndex,row in enumerate(ReferencePatternDesiredTuningDict['provided_reference_patterns']):
            massFragment = row[0]
            truncatedIntensities = numpy.take(row[1:],ReferencePatternDesiredTuningDict['OverlappingIndices'])   
            rowTruncated = numpy.insert(truncatedIntensities, 0, massFragment)
            ReferencePatternDesiredTuningDict['overlapping_provided_reference_patterns'].append(rowTruncated)    
    #Note that the columns of overlapping_provided_reference_patterns are not just truncated but also rearranged. 
    #Thus, going forward, OverlappingMolecules must be used for indices rather than using molecules.
  
   
    '''
    Step 3b: Determine Matching Mass fragments
    '''    
    #Only append rows that have an overlapping mass fragment.
    matchingFragmentsOnlyPatterns = []
    for row in ReferencePatternExistingTuningDict['overlapping_provided_reference_patterns']:
        if row[0] in OverlappingFragments:
            matchingFragmentsOnlyPatterns.append(row)
    #now swap in the new list.
    ReferencePatternExistingTuningDict['overlapping_provided_reference_patterns'] = matchingFragmentsOnlyPatterns    

    #Only append rows that have an overlapping mass fragment.
    matchingFragmentsOnlyPatterns = []
    for row in ReferencePatternDesiredTuningDict['overlapping_provided_reference_patterns']:
        if row[0] in OverlappingFragments:
            matchingFragmentsOnlyPatterns.append(row)
    #now swap in the new list.
    ReferencePatternDesiredTuningDict['overlapping_provided_reference_patterns'] = matchingFragmentsOnlyPatterns

    #Convert the two lists to numpy arrays so that they can be divided....
    ReferencePatternExistingTuningDict['overlapping_provided_reference_patterns'] = numpy.array(ReferencePatternExistingTuningDict['overlapping_provided_reference_patterns'])
    ReferencePatternDesiredTuningDict['overlapping_provided_reference_patterns'] = numpy.array(ReferencePatternDesiredTuningDict['overlapping_provided_reference_patterns'])
    #Need to standardize again, because the ratios during division below will be wrong if the largest peak was not among the matching fragments.
    ReferencePatternExistingTuningDict['overlapping_provided_reference_patterns'] = StandardizeReferencePattern(ReferencePatternExistingTuningDict['overlapping_provided_reference_patterns'],len(OverlappingMolecules)) #this does have the molecular weight as the first column.
    ReferencePatternDesiredTuningDict['overlapping_provided_reference_patterns'] = StandardizeReferencePattern(ReferencePatternDesiredTuningDict['overlapping_provided_reference_patterns'],len(OverlappingMolecules)) #this does have the molecular weight as the first column.

    '''
    Find a,b,c:
    '''
    numpy.seterr(divide='ignore', invalid='ignore') #This is to prevent some unexpected (but not at the moment concerning) warnings in the below lines from the nan values.
    RatioOfPatterns = ReferencePatternDesiredTuningDict['overlapping_provided_reference_patterns'][:,1:]/ReferencePatternExistingTuningDict['overlapping_provided_reference_patterns'][:,1:]    
    #For each row, excluding the Take the mean excluding nan values.
    meanRatioPerMassFragment = numpy.nanmean(RatioOfPatterns,1) #the 1 is for over axis 1, not over axix 0.    
    numpy.seterr(divide='warn', invalid='warn') #This is to [it tje warnings back how they normally are.
    from numpy import inf
    #Due to divide by zero errors, it is useful to get rid of any infinity type values.
        #Below 3 lines were first attempt to deal with infinte values, but were not necessary. Can skip to the fitting.
    meanRatioPerMassFragment[meanRatioPerMassFragment==0] = 'nan' #zeros are not useful.
        # meanRatioPerMassFragment[meanRatioPerMassFragment==inf] = 'nan'
        # meanRatioPerMassFragment[meanRatioPerMassFragment==-inf] = 'nan'
    #Following this post.. https://stackoverflow.com/questions/28647172/numpy-polyfit-doesnt-handle-nan-values
    FiniteValueIndices = numpy.isfinite(OverlappingFragments) & numpy.isfinite(meanRatioPerMassFragment)
    #But changing the next lines to add in temporary variables.
    FiniteOverlappingFragments = OverlappingFragments[FiniteValueIndices]
    if len(FiniteOverlappingFragments) == 0:
        return [0,0,1], [[0,0,0],[0,0,0],[0,0,0]]                                       
    FiniteMeanRatioPerMassFragment = meanRatioPerMassFragment[FiniteValueIndices]
    abcCoefficients, abcCoefficients_cov =numpy.polyfit(FiniteOverlappingFragments,FiniteMeanRatioPerMassFragment,2, cov=True) #The two is for 2nd degree polynomial.
    FiniteMeanReverseRatioPerMassFragment = 1/FiniteMeanRatioPerMassFragment
    abcCoefficients_reverse, abcCoefficients_reverse_cov =numpy.polyfit(FiniteOverlappingFragments,FiniteMeanReverseRatioPerMassFragment,2, cov=True) #The two is for 2nd degree polynomial.
    
    if exportCoefficients == True:
        numpy.savetxt('TuningCorrectorPattern_FiniteOverlappingFragments.csv', FiniteOverlappingFragments, delimiter=",")
        numpy.savetxt('TuningCorrectorPattern_FiniteMeanRatioPerMassFragment.csv', FiniteMeanRatioPerMassFragment, delimiter=",")
        #Export the TuningCorrectorCoefficients as a list.
        with open('TuningCorrectorCoefficients.txt', 'w') as the_file:
            the_file.write(str(list(abcCoefficients))) 
        try:         #Export the TuningCorrectorCoefficients_cov to a csv.
            numpy.savetxt('TuningCorrectorCoefficients_cov.csv', abcCoefficients_cov, delimiter=",")
        except:
            print("Could not save to file TuningCorrectorCoefficients_cov")
        #Do the same for the reverse relation.
        numpy.savetxt('TuningCorrectorPattern_FiniteMeanReverseRatioPerMassFragment.csv', FiniteMeanReverseRatioPerMassFragment, delimiter=",")
        with open('TuningCorrectorCoefficientsReverse.txt', 'w') as the_file:
            the_file.write(str(list(abcCoefficients_reverse))) 
        #Export the TuningCorrectorCoefficients_cov to a csv.
        try:
            numpy.savetxt('TuningCorrectorCoefficientsReverse_cov.csv', abcCoefficients_reverse_cov, delimiter=",")
        except:
            print("Could not save to file TuningCorrectorCoefficientsReverse_cov")


    
    #Factor = A*X^2 + B*X + C, so C=1.0 means the factor is 1.0 and independent of molecular weight.
    #To use this with mixed patterns (meaning, some patterns taken from Literature reference like NIST) you will need to first determine the A,B,C coefficients, then you'll have to divide the Literature reference pattern by this factor. This will compensate for when the code multiplies by the factor, thus putting the mixed patterns into the same tuning.
    return abcCoefficients, abcCoefficients_cov
 
#this function either creates or gets the three coefficients for the polynomial correction (Tuning Correction) and calculates
#the correction factor for the relative intensities of each mass fragment, outputting a corrected set
#of relative intensities
def TuningCorrector(referenceDataArrayWithAbscissa,referenceCorrectionCoefficients, referenceCorrectionCoefficients_cov, referenceFileExistingTuningAndForm=None,referenceFileDesiredTuningAndForm=None,measuredReferenceYorN="no"):
    #Tuning corrector is designed to work with standardized_reference_patterns, so first we make sure standardize the data.
    referenceDataArrayWithAbscissa=StandardizeReferencePattern(referenceDataArrayWithAbscissa)    
    if type(referenceCorrectionCoefficients) == type({}):#check if it's a dictionary. If it is, we need to make it a list.
        referenceCorrectionCoefficients = [referenceCorrectionCoefficients['A'],referenceCorrectionCoefficients['B'],referenceCorrectionCoefficients['C']]
    if measuredReferenceYorN =='yes':
        if len(referenceFileDesiredTuningAndForm) == 0:#TODO: this isn't very good logic, but it allows automatic population of referenceFileDesiredTuningAndForm. The problem is it is reading from file again instead of using the already made ReferenceData object. ABCDetermination and possibly TuningCorrector should be changed so that it can take *either* a ReferenceData object **or** a ReferenceData filename. The function can check if it is receiving a string, and if it's not receiving a string it can assume it's receiving an object.
            if '.csv' in G.referenceFileNamesList[0]:                
               referenceFileDesiredTuningAndForm = [ "ExportedDesiredTuningReferencePattern.csv","xyyy" ] #Take the first item from G.referenceFileNamesList and from G.referenceFormsList.
            if '.tsv' in G.referenceFileNamesList[0]:                
               referenceFileDesiredTuningAndForm = [ "ExportedDesiredTuningReferencePattern.tsv","xyyy" ] #Take the first item from G.referenceFileNamesList and from G.referenceFormsList.
        abcCoefficients, abcCoefficients_cov = ABCDetermination(referenceFileExistingTuningAndForm,referenceFileDesiredTuningAndForm)
        referenceCorrectionCoefficients[0],referenceCorrectionCoefficients[1],referenceCorrectionCoefficients[2]= abcCoefficients
        G.referenceCorrectionCoefficients = referenceCorrectionCoefficients #TODO: Maybe this logic should be changed, since it will result in an exporting of the last coefficients used, whether a person is doing forward tuning or reverse tuning.
        referenceCorrectionCoefficients_cov = abcCoefficients_cov
        G.referenceCorrectionCoefficients_cov = referenceCorrectionCoefficients_cov
    
    referenceabscissa = referenceDataArrayWithAbscissa[:,0] #gets arrays of just data and abscissa
    referenceDataArray = referenceDataArrayWithAbscissa[:,1:]
    nonZeroValueLocations = referenceDataArray > 0 #this returns booleans of "True" and "False", which are the same as 0 and 1 during multiplication.
    referenceDataArray_tuning_uncertainties = referenceDataArray*0.0 #just initializing.
    if list(referenceCorrectionCoefficients) != [0,0,1]:                                                                                    
        for massfrag_counter in range(len(referenceabscissa)):#array-indexed for loop, only the data is altered, based on the abscissa (mass-dependent correction factors)
            factor = referenceCorrectionCoefficients[0]*(referenceabscissa[massfrag_counter]**2)  + referenceCorrectionCoefficients[1]*referenceabscissa[massfrag_counter]+referenceCorrectionCoefficients[2] #obtains the factor from molecular weight of abscissa
            referenceDataArray[massfrag_counter,:] = referenceDataArray[massfrag_counter,:]*factor
            if type(referenceCorrectionCoefficients_cov) != None:
                if sum(numpy.array(referenceCorrectionCoefficients_cov)).all()!=0:
                    if len(numpy.shape(referenceCorrectionCoefficients_cov)) == 1 and (len(referenceCorrectionCoefficients_cov) > 0): #If it's a 1D array/list that is filled, we'll diagonalize it.
                        referenceCorrectionCoefficients_cov = np.diagflat(referenceCorrectionCoefficients_cov) 
                    #Now we need to use a function from XYYYDataFunctionsSG
                    factor_uncertainty = DataFunctions.returnPolyvalEstimatedUncertainties(referenceabscissa[massfrag_counter], referenceCorrectionCoefficients_cov, referenceCorrectionCoefficients_cov)
                    referenceDataArray_tuning_uncertainties[massfrag_counter,:]=referenceDataArray[massfrag_counter,:]*factor_uncertainty                                                                                                                  
    # referenceDataArrayWithAbscissa[:,0] = referenceabscissa
    # referenceDataArrayWithAbscissa[:,1:] = referenceDataArray #This is actually already occuring above because it's a pointer, but this line is just to make more clear what has happened.
    referenceDataArrayWithAbscissa[:,1:] = nonZeroValueLocations*referenceDataArray #we are setting things to zero if they were originally zero (this multiplies the nonzero locations by 1, and everything else by 0)
    referenceDataArrayWithAbscissa[:,1:] = (referenceDataArrayWithAbscissa[:,1:] > 0)*referenceDataArrayWithAbscissa[:,1:]  #multiplying by 1 where the final value is >0, and by 0 everywhere else.
    referenceDataArrayWithAbscissa_tuning_uncertainties = referenceDataArrayWithAbscissa*1.0 #creating copy with abscissa, then will fill.
    #The referenceDataArray_tuning_uncertainties are absolute uncertainties (note that above the factor_uncertainty is multiplied by the actual value)
    referenceDataArrayWithAbscissa_tuning_uncertainties[:,1:] = referenceDataArray_tuning_uncertainties
    return referenceDataArrayWithAbscissa, referenceDataArrayWithAbscissa_tuning_uncertainties 


'''Makes a mixed reference pattern from two reference patterns, including tuning correction.'''
def createReferencePatternWithTuningCorrection(ReferenceData, verbose=True, returnMixedPattern=False):
    #If tuning corrector is off (as of Nov 2021, variable measuredReferenceYorN), then we return ReferenceData unchanged.
    if G.measuredReferenceYorN =='no': 
        return ReferenceData
    
    
    # standardize the reference data columns such that the maximum intensity value per molecule is 100 and everything else is scaled to that maximum value.
    ReferenceData.ExportCollector('StandardizedReferencePattern', use_provided_reference_patterns=False)
    if ReferenceData.referenceFileNameExtension == 'csv':
        ReferenceData.exportReferencePattern('ExportedReferencePatternOriginalAnalysis.csv')
    if ReferenceData.referenceFileNameExtension == 'tsv':
        ReferenceData.exportReferencePattern('ExportedReferencePatternOriginalAnalysis.tsv')

    if G.calculateUncertaintiesInConcentrations == True: 
        if type(G.referenceFileUncertainties) != type(None):
            ReferenceData.update_relative_standard_uncertainties()

    if verbose:
        print('beginning TuningCorrector')
    if type(G.referenceCorrectionCoefficients) == type({}):#check if it's a dictionary.
        G.referenceCorrectionCoefficients = [G.referenceCorrectionCoefficients['A'],G.referenceCorrectionCoefficients['B'],G.referenceCorrectionCoefficients['C']]
    #only apply the tuning correction if the list is not 0 0 1.
    try:
        type(G.referenceCorrectionCoefficients_cov)#If it doesn't exist, we'll make it a None type.
    except:    
        G.referenceCorrectionCoefficients_cov = None
    #Wanted to do something like "if list(G.referenceCorrectionCoefficients) != [0,0,1]:" but can't do it out here. Can do it inside function.       

    if G.measuredReferenceYorN =='no':
       G.createMixedTuningPattern =  False #override the mixed tuning pattern choice if there is no measured reference.

    resetReferenceFileDesiredTuningAndForm = False   #initializing. At end, will reset G.referenceFileDesiredTuningAndForm if it was originally blank.

    if hasattr(G, 'tuningCorrectPatternInternalVsExternal') == False: #Check if the tuningCorrectPatternInternalVsExternal attribute has been filled. If not, then fill it with the default which is 'External'
        G.tuningCorrectPatternInternalVsExternal = 'External' 
    if ((G.createMixedTuningPattern == False) or (G.tuningCorrectPatternInternalVsExternal == 'Internal')):   
        if G.tuningCorrectPatternInternalVsExternal == 'Internal':
            if ((len(G.referenceFileExistingTuningAndForm) == 0) and (len(G.referenceFileStandardTuningAndForm) == 0)): #if no external pattern provided, we will use the coefficients provided directly.
                G.referenceCorrectionCoefficients = G.referenceCorrectionCoefficients
                G.referenceCorrectionCoefficients_cov = G.referenceCorrectionCoefficients_cov
            elif ((G.referenceFileExistingTuningAndForm !=[]) or (G.referenceFileStandardTuningAndForm!=[])) : #in this case, we are going to overwrite any coefficients provided in order to apply the desired tuning to the external pattern.
                referenceFileDesiredTuningAndForm = G.referenceFileDesiredTuningAndForm
                referenceFileExistingTuningAndForm = G.referenceFileExistingTuningAndForm
                if len(referenceFileExistingTuningAndForm) == 0: #Use the standard tuning file if blank.
                    referenceFileExistingTuningAndForm = G.referenceFileStandardTuningAndForm
                ReferenceDataExistingTuning = createReferenceDataObject ( referenceFileExistingTuningAndForm[0],referenceFileExistingTuningAndForm[1], AllMID_ObjectsDict=G.AllMID_ObjectsDict)
                if ReferenceDataExistingTuning.referenceFileNameExtension =='csv':
                    ReferenceDataExistingTuning.exportReferencePattern('ExportedReferencePatternExistingOriginal.csv')
                if ReferenceDataExistingTuning.referenceFileNameExtension =='tsv':
                    ReferenceDataExistingTuning.exportReferencePattern('ExportedReferencePatternExistingOriginal.tsv')
                if len(referenceFileDesiredTuningAndForm) == 0:#TODO: this isn't very good logic, but it allows automatic population of referenceFileDesiredTuningAndForm. The problem is it is reading from file again instead of using the already made ReferenceData object. ABCDetermination and possibly TuningCorrector should be changed so that it can take *either* a ReferenceData object **or** a ReferenceData filename. The function can check if it is receiving a string, and if it's not receiving a string it can assume it's receiving an object.
                    resetReferenceFileDesiredTuningAndForm = True 
                    if ReferenceDataExistingTuning.referenceFileNameExtension =='csv':
                        referenceFileDesiredTuningAndForm = [ "ExportedReferencePatternOriginalAnalysis.csv","xyyy" ] #Take the first item from G.referenceFileNamesList and from G.referenceFormsList.
                    if ReferenceDataExistingTuning.referenceFileNameExtension =='tsv':
                        referenceFileDesiredTuningAndForm = [ "ExportedReferencePatternOriginalAnalysis.tsv","xyyy" ]
                abcCoefficients, abcCoefficients_cov = ABCDetermination(referenceFileExistingTuningAndForm,referenceFileDesiredTuningAndForm)
                referenceCorrectionCoefficients = numpy.zeros(3)
                referenceCorrectionCoefficients[0],referenceCorrectionCoefficients[1],referenceCorrectionCoefficients[2]= abcCoefficients
                G.referenceCorrectionCoefficients = referenceCorrectionCoefficients #TODO: Maybe this logic should be changed, since it will result in an exporting of the last coefficients used, whether a person is doing forward tuning or reverse tuning.
                referenceCorrectionCoefficients_cov = abcCoefficients_cov
                G.referenceCorrectionCoefficients_cov = referenceCorrectionCoefficients_cov

                            
            #if we are not returning a mixed pattern, we are applying the tuning correction directly to the original ReferenceData object.
            ReferenceData.standardized_reference_patterns_tuning_corrected, ReferenceData.standardized_reference_patterns_tuning_uncertainties = TuningCorrector(ReferenceData.standardized_reference_patterns,
                                                                   G.referenceCorrectionCoefficients,G.referenceCorrectionCoefficients_cov)


            #TuningCorrector un-standardizes the patterns, so the patterns have to be standardized again.
            ReferenceData.standardized_reference_patterns=StandardizeReferencePattern(ReferenceData.standardized_reference_patterns_tuning_corrected,len(ReferenceData.molecules))
            #add "TuningCorrected" to SourceOfFragmentationPatterns
            ReferenceData.SourceOfFragmentationPatterns = list(ReferenceData.SourceOfFragmentationPatterns)
            for sourceIndex in range(len(ReferenceData.SourceOfFragmentationPatterns)):
                ReferenceData.SourceOfFragmentationPatterns[sourceIndex] = ReferenceData.SourceOfFragmentationPatterns[sourceIndex] + '-TuningCorrected'
            ReferenceData.ExportCollector('TuningCorrector') #export the pattern.
                
            #Now check if uncertainties already exist, and if they do then the two uncertainties need to be combined. Else, made equal.
            try:
                ReferenceData.absolute_standard_uncertainties[:,1:] = (ReferenceData.absolute_standard_uncertainties[:,1:]**2 + ReferenceData.standardized_reference_patterns_tuning_uncertainties[:,1:]**2)**0.5
            except:
                ReferenceData.absolute_standard_uncertainties = ReferenceData.standardized_reference_patterns_tuning_uncertainties

            ReferenceData.ExportCollector('TuningCorrector_absolute_standard_uncertainties', export_tuning_uncertainties= True) #These are not yet actually standardized. Happens below.
            ReferenceData.ExportCollector('StandardizeReferencePattern_absolute_standard_uncertainties', export_standard_uncertainties= True)
            if G.calculateUncertaintiesInConcentrations == True: 
                if type(G.referenceFileUncertainties) != type(None):                                                      
                    ReferenceData.update_relative_standard_uncertainties()
                    ReferenceData.ExportCollector('StandardizeReferencePattern_relative_standard_uncertainties', export_relative_uncertainties= True)
        elif G.tuningCorrectPatternInternalVsExternal =='External':
            if ((len(G.referenceFileExistingTuningAndForm) == 0) and (len(G.referenceFileStandardTuningAndForm) == 0)): #if no external pattern provided, we will use the coefficients provided directly.
                G.referenceCorrectionCoefficients = G.referenceCorrectionCoefficients
                G.referenceCorrectionCoefficients_cov = G.referenceCorrectionCoefficients_cov
            elif ((G.referenceFileExistingTuningAndForm !=[]) or (G.referenceFileStandardTuningAndForm!=[])) : #in this case, we are going to overwrite any coefficients provided in order to apply the desired tuning to the external pattern.
                referenceFileDesiredTuningAndForm = G.referenceFileDesiredTuningAndForm
                referenceFileExistingTuningAndForm = G.referenceFileExistingTuningAndForm
                if len(referenceFileExistingTuningAndForm) == 0: #Use the standard tuning file if blank.
                        referenceFileExistingTuningAndForm = G.referenceFileStandardTuningAndForm
                if len(referenceFileDesiredTuningAndForm) == 0: #Use the original reference pattern if blank.
                        resetReferenceFileDesiredTuningAndForm = True
                        referenceFileDesiredTuningAndForm = [G.referenceFileNamesList[0],  G.referenceFormsList[0]] #Take the first item from G.referenceFileNamesList and from G.referenceFormsList.
                #We don't use the function GenerateReferenceDataList because that function does more than just making a reference object.
                ReferenceDataExistingTuning = createReferenceDataObject ( referenceFileExistingTuningAndForm[0],referenceFileExistingTuningAndForm[1], AllMID_ObjectsDict=G.AllMID_ObjectsDict)   
                if ReferenceDataExistingTuning.referenceFileNameExtension == 'csv':
                    ReferenceDataExistingTuning.exportReferencePattern('ExportedReferencePatternExistingOriginal.csv')
                if ReferenceDataExistingTuning.referenceFileNameExtension == 'tsv':
                    ReferenceDataExistingTuning.exportReferencePattern('ExportedReferencePatternExistingOriginal.tsv')
                ReferenceDataDesiredTuning = createReferenceDataObject ( referenceFileDesiredTuningAndForm[0],referenceFileDesiredTuningAndForm[1], AllMID_ObjectsDict=G.AllMID_ObjectsDict)   
                if ReferenceDataDesiredTuning.referenceFileNameExtension == 'csv':
                    ReferenceDataDesiredTuning.exportReferencePattern('ExportedReferencePatternDesiredOriginal.csv')
                if ReferenceDataDesiredTuning.referenceFileNameExtension == 'tsv':
                    ReferenceDataDesiredTuning.exportReferencePattern('ExportedReferencePatternDesiredOriginal.tsv')
                abcCoefficients, abcCoefficients_cov = ABCDetermination(referenceFileExistingTuningAndForm,referenceFileDesiredTuningAndForm)
                referenceCorrectionCoefficients = numpy.zeros(3)
                referenceCorrectionCoefficients[0],referenceCorrectionCoefficients[1],referenceCorrectionCoefficients[2]= abcCoefficients
                G.referenceCorrectionCoefficients = referenceCorrectionCoefficients #TODO: Maybe this logic should be changed, since it will result in an exporting of the last coefficients used, whether a person is doing forward tuning or reverse tuning.
                referenceCorrectionCoefficients_cov = abcCoefficients_cov
                G.referenceCorrectionCoefficients_cov = referenceCorrectionCoefficients_cov

                            
            #Now that the coefficients are obtained, we need to to a tuning correction on the external reference pattern.
            ReferenceDataExternalTuningCorrected = copy.deepcopy(ReferenceDataExistingTuning)
            ReferenceDataExternalTuningCorrected.standardized_reference_patterns_tuning_corrected, ReferenceDataExternalTuningCorrected.standardized_reference_patterns_tuning_uncertainties = TuningCorrector(  
                ReferenceDataExternalTuningCorrected.standardized_reference_patterns, G.referenceCorrectionCoefficients,G.referenceCorrectionCoefficients_cov)
            #add "TuningCorrected" to SourceOfFragmentationPatterns
            ReferenceDataExternalTuningCorrected.SourceOfFragmentationPatterns = list(ReferenceDataExternalTuningCorrected.SourceOfFragmentationPatterns)
            for sourceIndex in range(len(ReferenceDataExternalTuningCorrected.SourceOfFragmentationPatterns)):
                ReferenceDataExternalTuningCorrected.SourceOfFragmentationPatterns[sourceIndex] = ReferenceDataExternalTuningCorrected.SourceOfFragmentationPatterns[sourceIndex] + '-TuningCorrected'
      
            #TuningCorrector un-standardizes the patterns, so the patterns have to be standardized again.
            ReferenceDataExternalTuningCorrected.standardized_reference_patterns=StandardizeReferencePattern(ReferenceDataExternalTuningCorrected.standardized_reference_patterns_tuning_corrected)
            
            if ReferenceDataExternalTuningCorrected.referenceFileNameExtension == 'csv':
                ReferenceDataExternalTuningCorrected.exportReferencePattern('ExportedReferencePatternExternalTuningCorrected.csv')
            if ReferenceDataExternalTuningCorrected.referenceFileNameExtension == 'tsv':
                ReferenceDataExternalTuningCorrected.exportReferencePattern('ExportedReferencePatternExternalTuningCorrected.tsv')
                
            #Now check if uncertainties already exist, and if they do then the two uncertainties need to be combined. Else, made equal.
            #We don't export these uncertainties since this is just an intermediate step.
            try:
                ReferenceDataExternalTuningCorrected.absolute_standard_uncertainties = (ReferenceDataExternalTuningCorrected.absolute_standard_uncertainties**2 + ReferenceDataExternalTuningCorrected.standardized_reference_patterns_tuning_uncertainties**2)**0.5
            except:
                ReferenceDataExternalTuningCorrected.absolute_standard_uncertainties = ReferenceDataExternalTuningCorrected.standardized_reference_patterns_tuning_uncertainties
            if G.calculateUncertaintiesInConcentrations == True: 
                if type(G.referenceFileUncertainties) != type(None):                                                      
                    ReferenceDataExternalTuningCorrected.update_relative_standard_uncertainties()
                    ReferenceData.ExportCollector('StandardizeReferencePattern_relative_standard_uncertainties', export_relative_uncertainties= True)
                                                               
    elif ((G.createMixedTuningPattern== True) and (G.tuningCorrectPatternInternalVsExternal =='External')):  #in this case, we are going to apply the current tuning to the external pattern, and also create a mixed pattern. So the ReferenceData pointer will point to a mixed pattern by the end of this if statement.
        if G.measuredReferenceYorN =='yes':
            #First read in the existing tuning patterns.
            referenceFileDesiredTuningAndForm = G.referenceFileDesiredTuningAndForm
            referenceFileExistingTuningAndForm = G.referenceFileExistingTuningAndForm
            if len(referenceFileExistingTuningAndForm) == 0: #Use the standard tuning file if blank.
                    referenceFileExistingTuningAndForm = G.referenceFileStandardTuningAndForm
            if len(referenceFileDesiredTuningAndForm) == 0: #Use the original reference pattern if blank.
                    resetReferenceFileDesiredTuningAndForm = True
                    referenceFileDesiredTuningAndForm = [G.referenceFileNamesList[0],  G.referenceFormsList[0]] #Take the first item from G.referenceFileNamesList and from G.referenceFormsList.
            #We don't use the function GenerateReferenceDataList because that function does more than just making a reference object.
            ReferenceDataExistingTuning = createReferenceDataObject ( referenceFileExistingTuningAndForm[0],referenceFileExistingTuningAndForm[1], AllMID_ObjectsDict=G.AllMID_ObjectsDict)   
            if ReferenceDataExistingTuning.referenceFileNameExtension == 'csv':
                ReferenceDataExistingTuning.exportReferencePattern('ExportedReferencePatternExistingOriginal.csv')
            if ReferenceDataExistingTuning.referenceFileNameExtension == 'tsv':
                ReferenceDataExistingTuning.exportReferencePattern('ExportedReferencePatternExistingOriginal.tsv')            
            ReferenceDataDesiredTuning = createReferenceDataObject ( referenceFileDesiredTuningAndForm[0],referenceFileDesiredTuningAndForm[1], AllMID_ObjectsDict=G.AllMID_ObjectsDict)   
            if ReferenceDataDesiredTuning.referenceFileNameExtension == 'csv':
                ReferenceDataDesiredTuning.exportReferencePattern('ExportedReferencePatternDesiredOriginal.csv')
            if ReferenceDataDesiredTuning.referenceFileNameExtension == 'tsv':
                ReferenceDataDesiredTuning.exportReferencePattern('ExportedReferencePatternDesiredOriginal.tsv')
            abcCoefficients, abcCoefficients_cov = ABCDetermination(referenceFileExistingTuningAndForm,referenceFileDesiredTuningAndForm)
            referenceCorrectionCoefficients = numpy.zeros(3)
            referenceCorrectionCoefficients[0],referenceCorrectionCoefficients[1],referenceCorrectionCoefficients[2]= abcCoefficients
            G.referenceCorrectionCoefficients = referenceCorrectionCoefficients #TODO: Maybe this logic should be changed, since it will result in an exporting of the last coefficients used, whether a person is doing forward tuning or reverse tuning.
            referenceCorrectionCoefficients_cov = abcCoefficients_cov
            G.referenceCorrectionCoefficients_cov = referenceCorrectionCoefficients_cov
            
            #Now that the coefficients are obtained, we need to to a tuning correction on the external reference pattern.
            ReferenceDataExternalTuningCorrected = copy.deepcopy(ReferenceDataExistingTuning)
            ReferenceDataExternalTuningCorrected.standardized_reference_patterns_tuning_corrected, ReferenceDataExternalTuningCorrected.standardized_reference_patterns_tuning_uncertainties = TuningCorrector(  
                ReferenceDataExternalTuningCorrected.standardized_reference_patterns, G.referenceCorrectionCoefficients,G.referenceCorrectionCoefficients_cov)
            #add "TuningCorrected" to SourceOfFragmentationPatterns
            ReferenceDataExternalTuningCorrected.SourceOfFragmentationPatterns = list(ReferenceDataExternalTuningCorrected.SourceOfFragmentationPatterns)
            for sourceIndex in range(len(ReferenceDataExternalTuningCorrected.SourceOfFragmentationPatterns)):
                ReferenceDataExternalTuningCorrected.SourceOfFragmentationPatterns[sourceIndex] = ReferenceDataExternalTuningCorrected.SourceOfFragmentationPatterns[sourceIndex] + '-TuningCorrected'
      
            #TuningCorrector un-standardizes the patterns, so the patterns have to be standardized again.
            ReferenceDataExternalTuningCorrected.standardized_reference_patterns=StandardizeReferencePattern(ReferenceDataExternalTuningCorrected.standardized_reference_patterns_tuning_corrected)
            
            if ReferenceDataExternalTuningCorrected.referenceFileNameExtension == 'csv':
                ReferenceDataExternalTuningCorrected.exportReferencePattern('ExportedReferencePatternExternalTuningCorrected.csv')
            if ReferenceDataExternalTuningCorrected.referenceFileNameExtension == 'tsv':
                ReferenceDataExternalTuningCorrected.exportReferencePattern('ExportedReferencePatternExternalTuningCorrected.tsv')
                
            #Now check if uncertainties already exist, and if they do then the two uncertainties need to be combined. Else, made equal.
            #We don't export these uncertainties since this is just an intermediate step.
            try:
                ReferenceDataExternalTuningCorrected.absolute_standard_uncertainties = (ReferenceDataExternalTuningCorrected.absolute_standard_uncertainties**2 + ReferenceDataExternalTuningCorrected.standardized_reference_patterns_tuning_uncertainties**2)**0.5
            except:
                ReferenceDataExternalTuningCorrected.absolute_standard_uncertainties = ReferenceDataExternalTuningCorrected.standardized_reference_patterns_tuning_uncertainties
            if G.calculateUncertaintiesInConcentrations == True: 
                if type(G.referenceFileUncertainties) != type(None):                                                      
                    ReferenceDataExternalTuningCorrected.update_relative_standard_uncertainties()
                    ReferenceData.ExportCollector('StandardizeReferencePattern_relative_standard_uncertainties', export_relative_uncertainties= True)

            #Now to make the mixed reference pattern using a function that extends one reference pattern by another.
            #The extendReferencePattern function uses existing add molecule and remove molecules functions.
            #The function *automatically* does not extend by molecules that already exist.
            
            ReferenceData, addedReferenceSlice = extendReferencePattern(OriginalReferenceData=ReferenceData,ReferenceDataToExtendBy=ReferenceDataExternalTuningCorrected)
            #extendReferencePattern can't currently add the uncertainties sub-objects, so we do that in the lines below.
            ReferenceData.ExportAtEachStep = G.ExportAtEachStep
            if ReferenceData.referenceFileNameExtension == 'csv':
                ReferenceData.exportReferencePattern('ExportedReferencePatternMixed.csv')
            if ReferenceData.referenceFileNameExtension == 'tsv':
                ReferenceData.exportReferencePattern('ExportedReferencePatternMixed.tsv')            
                
            #If specific molecules is on, then we need to potentially remove molecules from the mixed reference pattern, unless PreparingAllMoleculesReferenceDataList.
            if PreparingAllMoleculesReferenceDataList == False: #THIS IS A GLOBAL FLAG. IT IS AN IMPLIED ARGUMENT.
                if G.specificMolecules.lower() == "yes":
                    moleculesToRemove = []
                    for moleculeName in ReferenceData.molecules:
                        if moleculeName not in G.chosenMoleculesNames:
                            moleculesToRemove.append(moleculeName)
                    if len(moleculesToRemove) > 0:
                        ReferenceData = ReferenceData.removeMolecules(moleculesToRemove)
                    if ReferenceData.referenceFileNameExtension == 'csv':
                        ReferenceData.exportReferencePattern("ExportedReferencePatternMixedChosenMolecules.csv")
                    if ReferenceData.referenceFileNameExtension == 'tsv':
                        ReferenceData.exportReferencePattern("ExportedReferencePatternMixedChosenMolecules.tsv")

                    
                # else:
                    # if ReferenceDataList[ReferenceDataIndex].ExportAtEachStep == 'yes':
                        # ReferenceDataList[ReferenceDataIndex].ExportCollector('StandardizedReferencePattern' + str(ReferenceDataIndex), use_provided_reference_patterns=False)
                    # ReferenceDataList[ReferenceDataIndex].exportReferencePattern("TuningCorrectorMixedPattern" + str(ReferenceDataIndex) + ".csv")


            # #Note: it is assumed that the relative_standard_uncertainties correspond to original reference, so before tuning corrector, thus we do not recalculate that factor.
            # ReferenceData.ExportCollector('StandardizeReferencePattern_absolute_standard_uncertainties', export_standard_uncertainties= True)

            # ReferenceData.update_relative_standard_uncertainties()
            # ReferenceData.ExportCollector('StandardizeReferencePattern_absolute_standard_uncertainties', export_standard_uncertainties= True)
            # ReferenceData.ExportCollector('StandardizeReferencePattern_relative_standard_uncertainties', export_relative_uncertainties= True)
    else:
        print("Error: Incompatible Tuning Choices Detected. Mixed Tuning Pattern is currently only compatible with changing the external pattern.")
    if resetReferenceFileDesiredTuningAndForm == True: #Reset this variable if needed.
        G.referenceFileDesiredTuningAndForm = []               
    return ReferenceData
    

    
#tuningCorrectorGasMixture will correct all items in the ReferenceDataList, but uses a single tuning correction that it applies to each of the ReferenceData objects in the list.
#This function takes a measured gas mixture set of signals, then creates simulated signals (from NIST patterns, for example) to get a tuning correction for the spectrometer.
def tuningCorrectorGasMixture(ReferenceDataList, G, ExperimentData=None): #making it clear that there are UserInput choices that are needed for this feature.
        if ReferenceDataList[0].referenceFileNameExtension == 'csv':
            ReferenceDataList[0].exportReferencePattern('ExportedReferencePatternOriginalAnalysis.csv')
        if ReferenceDataList[0].referenceFileNameExtension == 'tsv':
            ReferenceDataList[0].exportReferencePattern('ExportedReferencePatternOriginalAnalysis.tsv')
        
        #First read in the existing tuning patterns.
        #We don't use the function GenerateReferenceDataList because that function does more than just making a reference object.
        ReferenceDataExistingTuning = createReferenceDataObject ( G.referenceFileExistingTuningAndForm[0],G.referenceFileExistingTuningAndForm[1], AllMID_ObjectsDict=G.AllMID_ObjectsDict)   
        ReferenceDataExistingTuning.ExportAtEachStep = G.ExportAtEachStep
    
        #Copying what is in GenerateReferenceDataList and ReferenceInputPreProcessing, we need to add the following block of code in if we want to allow uncertainties.   
        if G.calculateUncertaintiesInConcentrations == True:
            if type(G.referenceFileUncertainties) != type(None):
                if type(G.referenceFileUncertainties) == type(float(5)) or  type(G.referenceFileUncertainties) == type(int(5)) :
                    #TODO: Low priority. The below results in "nan" values. It could be better to change it to make zeros using a numpy "where" statement.
                    G.referenceFileUncertainties = float(G.referenceFileUncertainties) #Make sure we have a float.
                    #Get what we need.
                    provided_reference_patterns = ReferenceDataExistingTuning.provided_reference_patterns
                    provided_reference_patterns_without_masses = ReferenceDataExistingTuning.provided_reference_patterns[:,1:] #[:,0] is mass fragments, so we slice to remove those. 
                    #Make our variables ready.
                    absolute_standard_uncertainties = provided_reference_patterns*1.0 #First we make a copy.
                    absolute_standard_uncertainties_without_masses = (provided_reference_patterns_without_masses/provided_reference_patterns_without_masses)*G.referenceFileUncertainties #We do the division to get an array of ones. Then we multiply to get an array of the same size.
                    #now we populate our copy's non-mass area.
                    absolute_standard_uncertainties[:,1:] = absolute_standard_uncertainties_without_masses                                        
                    ReferenceDataExistingTuning.absolute_standard_uncertainties = absolute_standard_uncertainties
                    #We can't convert to relative uncertainties yet because the file may not be standardized yet.
                if type(G.referenceFileUncertainties) == type('string'):
                    if '.csv' in G.referenceFileExistingTuningAndForm[0]:
                        provided_reference_patterns_absolute_uncertainties, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form = readReferenceFile(G.referenceFileExistingTuningAndForm[0][:-4]+"_absolute_uncertainties.csv",G.referenceFileExistingTuningAndForm[1])
                    if '.tsv' in G.referenceFileExistingTuningAndForm[0]:
                        provided_reference_patterns_absolute_uncertainties, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form = readReferenceFile(G.referenceFileExistingTuningAndForm[0][:-4]+"_absolute_uncertainties.tsv",G.referenceFileExistingTuningAndForm[1])
                    ReferenceDataExistingTuning.absolute_standard_uncertainties = provided_reference_patterns_absolute_uncertainties #Just initializing the variable before filling it properly.
                    maximum_absolute_intensities = numpy.amax(ReferenceDataExistingTuning.provided_reference_patterns[:,1:], axis = 0) #Find the maximum intensity for each molecule.
                    ReferenceDataExistingTuning.absolute_standard_uncertainties[:,1:] = 100*ReferenceDataExistingTuning.absolute_standard_uncertainties[:,1:]/maximum_absolute_intensities
        if G.calculateUncertaintiesInConcentrations == True: 
            if type(G.referenceFileUncertainties) != type(None):
                ReferenceDataExistingTuning.update_relative_standard_uncertainties()
        
        #TODO: For above function call, Still need to put the last argument in later which is the ionization information: AllMID_ObjectsDict={})
        #Currently, in a "regular" MSRESOLVE run, the calculation of the correction values normally occurs inside  the ReferenceInputPreProcessing. That is a problem because 
        #the funciton ReferenceInputPreProcessing applies a tuning correction. So we need to get the CorrectionValues directly.
        ReferenceDataExistingTuning.correction_values, ReferenceDataExistingTuning.correction_values_relative_uncertainties = CorrectionValuesObtain(ReferenceDataExistingTuning)
        #We don't need a desired tuning data object right now, but we will make one since we'll need it for creating the mixed reference pattern later.
        #By default, the regular reference data object will be considered the desired tuning, if not specified.
        if len(G.referenceFileDesiredTuningAndForm) == 0:
            #Take the first item from ReferenceDataList (which may have been extracted from the experiment file)
            ReferenceDataDesiredTuning = ReferenceDataList[0]
            if ReferenceDataDesiredTuning.referenceFileNameExtension == 'csv':
                ReferenceDataDesiredTuning.exportReferencePattern("ExportedReferencePatternOriginal.csv") #Exported for during TuningCorrection.     
                referenceFileDesiredTuningAndForm = ["ExportedReferencePatternOriginal.csv", "xxxy"]
            if ReferenceDataDesiredTuning.referenceFileNameExtension == 'tsv':
                ReferenceDataDesiredTuning.exportReferencePattern("ExportedReferencePatternOriginal.tsv") #Exported for during TuningCorrection.     
                referenceFileDesiredTuningAndForm = ["ExportedReferencePatternOriginal.tsv", "xxxy"]
        else:
            ReferenceDataDesiredTuning = createReferenceDataObject ( G.referenceFileDesiredTuningAndForm[0],G.referenceFileDesiredTuningAndForm[1], AllMID_ObjectsDict=G.AllMID_ObjectsDict)   
            ReferenceDataDesiredTuning.ExportAtEachStep = G.ExportAtEachStep

        #We are going to make the simulated gas mixture first, because that turns out to be easier (for data structure reasons) when we want to extract the gas mixture signals from the data.
        
        #We need to make prepare an appropriately made list for the concentrations to use in simulating data using the ExistingTuningReferencePattern 
        #The concentrations must match the molecular order of the literature file.
        #  Firstly, the literature file can have molecules not part of the concentrations, so we should add concentrations of zero for those. (It's easier than removing those molecules from the literature object).
        #  Secondly, we need to re-order the concentrations to match the those of the existing reference pattern.
        knownConcentrationsArray = list(numpy.zeros(len(ReferenceDataExistingTuning.molecules))) #initializing to number of molecules. Making it a list so we can insert later.
        #The desired molecules list is the one in ReferenceDataExistingTuning.molecules, so we will loop across those.
        for moleculeIndex,moleculeName in enumerate(ReferenceDataExistingTuning.molecules):
            desiredConcentrationIndex = moleculeIndex  
            #Take the concentration if it's named:
            if moleculeName in G.UserChoices['measuredReferenceYorN']['tuningCorrectorGasMixtureMoleculeNames']:
                currentConcentrationIndex = list(G.UserChoices['measuredReferenceYorN']['tuningCorrectorGasMixtureMoleculeNames']).index(moleculeName)
                knownConcentrationsArray[desiredConcentrationIndex] = G.UserChoices['measuredReferenceYorN']['tuningCorrectorGasMixtureConcentrations'][currentConcentrationIndex]
            else: #else set it to zero.
                knownConcentrationsArray[desiredConcentrationIndex] = 0
        #The first item in the concentrations array is supposed to be the time. We will simply put the integer "1" there, since we will have one point. 
        knownConcentrationsArray.insert(0, 1)
        #knownConcentrationsArray needs to be nested in a numpy array for expected dimensionality when using RawSignalsSimulation
        knownConcentrationsArray = numpy.array([knownConcentrationsArray])
        #Before simulation, we also need the reciprocal_matching_correction_values array. In order to make the reciprocal_matching_correction_values array, we need to know which masses we need. We can actually just simulate all of the masses from the existingReferencePattern, and then let tuning corrector use whichever masses are useful
        #Below we directly call Populate_reciprocal_matching_correction_values because PrepareReferenceObjectsAndCorrectionValues could potentially apply an undesired "before comparisons" tuning factor correction. We will need this done before we simulate any signals.
        ReferenceDataExistingTuningMassFragments = ReferenceDataExistingTuning.standardized_reference_patterns[:,0]
        ReferenceDataExistingTuning = Populate_reciprocal_matching_correction_values(ReferenceDataExistingTuningMassFragments, ReferenceDataExistingTuning)
        #Now need to make the inputs for simulating raw signals of the gas mixture. A properly ordered and formatted concentration array, as well as properly formatted reciprocal_matching_correction_values.
        #reciprocal_matching_correction_values needs to be nested in a numpy array for expected dimensionality when using RawSignalsSimulation
        reciprocal_matching_correction_values_array = numpy.array([ReferenceDataExistingTuning.reciprocal_matching_correction_values]) 
        #now need to do the actual simulation of the gas mixture signals.
        simulateddata = RawSignalsSimulation(knownConcentrationsArray, reciprocal_matching_correction_values_array)
        ExportXYYYData("TuningCorrectorGasMixtureHypotheticalSimulatedSignals.csv", simulateddata, ReferenceDataExistingTuningMassFragments, abscissaHeader = "Time", fileSuffix = G.iterationSuffix, dataType = 'simulated')
        #now need to make a fake reference file for that. We'll make an MSReferenceObject and then use the internally built exportReferencePattern class function.
        simulateddata_intensities = simulateddata[0][1:] #The "0" is because it's a nested object. THe slicing is becuase the abscissa (time) value needs to be removed.
        #Now to transpose the simulateddata to make the simulated_reference_pattern
        simulated_reference_pattern = numpy.vstack((ReferenceDataExistingTuningMassFragments, simulateddata_intensities)).transpose()
        ReferenceDataTuningCorrectorGasMixtureSimulatedHypothetical = MSReference(simulated_reference_pattern, electronnumbers=[1], molecules=["GasMixture"], molecularWeights=[1], SourceOfFragmentationPatterns=["simulated"], sourceOfIonizationData=["referenceFileExistingTuningAndForm"], relativeIonizationEfficiencies=['GasMixture'], moleculeIonizationType=["GasMixture"]) #We fill some variables with the word "GasMixture" because there is no single ionization factor for the gas mixture.
        ReferenceDataTuningCorrectorGasMixtureSimulatedHypothetical.referenceFileNameExtension = ReferenceDataList[0].referenceFileNameExtension #This needs to be populated since there was no filename initially.
        if ReferenceDataTuningCorrectorGasMixtureSimulatedHypothetical.referenceFileNameExtension == 'csv':
            ReferenceDataTuningCorrectorGasMixtureSimulatedHypothetical.exportReferencePattern("ExportedReferencePatternGasMixtureSimulatedHypothetical.csv")
        if ReferenceDataTuningCorrectorGasMixtureSimulatedHypothetical.referenceFileNameExtension == 'tsv':
            ReferenceDataTuningCorrectorGasMixtureSimulatedHypothetical.exportReferencePattern("ExportedReferencePatternGasMixtureSimulatedHypothetical.tsv")
        
        if len(G.UserChoices['measuredReferenceYorN']['tuningCorrectorGasMixtureSignals']) == 0: #if the length is zero, we will populate this with the max and min of the times.
            G.UserChoices['measuredReferenceYorN']['tuningCorrectorGasMixtureSignals'] = [ min(ExperimentData.times), max(ExperimentData.times)]
        if len(G.UserChoices['measuredReferenceYorN']['tuningCorrectorGasMixtureSignals']) == 2: #We assume in this case that we have a pair of times.
            ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical = copy.deepcopy(ReferenceDataTuningCorrectorGasMixtureSimulatedHypothetical) #just initializing.
            ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.SourceOfFragmentationPatterns = ["measured"]
            #Get the mass fragments we need.
            expectedMassFragments = ReferenceDataTuningCorrectorGasMixtureSimulatedHypothetical.standardized_reference_patterns[:,0]
            
            #need to get the masses in common between the experiment data and ReferenceDataExistingTuningMassFragments            
            commonMF = []
            for massValue in numpy.array(ExperimentData.mass_fragment_numbers, dtype=float):
                if massValue in numpy.array(expectedMassFragments, dtype=float): #the 0 is because this object is nested.
                    commonMF.append(massValue)

            #Fill all of the fragments intensities with ones, this will make the extraction easier.
            ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.provided_reference_patterns[:,1:] = numpy.ones(numpy.shape(ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.provided_reference_patterns[:,1:]))
            ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.provided_mass_fragments = ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.provided_reference_patterns[:,0]
            #note that the last two arguments need to be further nested.
            ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.standardized_reference_patterns,ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical_uncertainties  = ExtractReferencePatternFromData(ExperimentData=ExperimentData, referenceData=ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical,rpcChosenMolecules=["GasMixture"],rpcChosenMoleculesMF=[commonMF],rpcTimeRanges=[G.UserChoices['measuredReferenceYorN']['tuningCorrectorGasMixtureSignals']])
            
            if ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.referenceFileNameExtension == 'csv':
                ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.exportReferencePattern("ExportedReferencePatternGasMixtureMeasuredHypothetical.csv")
            if ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.referenceFileNameExtension == 'tsv':
                ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.exportReferencePattern("ExportedReferencePatternGasMixtureMeasuredHypothetical.tsv")
        
        if G.measuredReferenceYorN =='yes':
            #we are using specific files for obtaining the tuning correction coefficients in the case of GasMixtureTuningCorrector
            if ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.referenceFileNameExtension == 'csv':
                referenceFileExistingTuningAndForm=["ExportedReferencePatternGasMixtureSimulatedHypothetical.csv","XYYY"]
                referenceFileDesiredTuningAndForm=["ExportedReferencePatternGasMixtureMeasuredHypothetical.csv", "XYYY"]
            if ReferenceDataTuningCorrectorGasMixtureMeasuredHypothetical.referenceFileNameExtension == 'tsv':
                referenceFileExistingTuningAndForm=["ExportedReferencePatternGasMixtureSimulatedHypothetical.tsv","XYYY"]
                referenceFileDesiredTuningAndForm=["ExportedReferencePatternGasMixtureMeasuredHypothetical.tsv", "XYYY"]
            abcCoefficients, abcCoefficients_cov = ABCDetermination(referenceFileExistingTuningAndForm,referenceFileDesiredTuningAndForm)
            referenceCorrectionCoefficients = numpy.zeros(3)
            referenceCorrectionCoefficients[0],referenceCorrectionCoefficients[1],referenceCorrectionCoefficients[2]= abcCoefficients
            G.referenceCorrectionCoefficients = referenceCorrectionCoefficients #TODO: Maybe this logic should be changed, since it will result in an exporting of the last coefficients used, whether a person is doing forward tuning or reverse tuning.
            referenceCorrectionCoefficients_cov = abcCoefficients_cov
            G.referenceCorrectionCoefficients_cov = referenceCorrectionCoefficients_cov
        ReferenceDataExistingTuningAfterCorrection = copy.deepcopy(ReferenceDataExistingTuning) #just initializing here, then actual tuning correction occurs in next lines.                
        #Tuning corrector does not operate on a ReferenceData object, it operates on standardized reference patterns.
        ReferenceDataExistingTuningAfterCorrection.standardized_reference_patterns, ReferenceDataExistingTuningAfterCorrection.standardized_reference_patterns_tuning_uncertainties = TuningCorrector(ReferenceDataExistingTuning.standardized_reference_patterns, G.referenceCorrectionCoefficients, G.referenceCorrectionCoefficients_cov)
        # TuningCorrectorGasMixtureCorrectedReferenceDataObject = MSReference(ReferenceDataExistingTuning.standardized_reference_patterns, electronnumbers=ReferenceDataExistingTuning.electronnumbers, molecules=ReferenceDataExistingTuning.molecules, molecularWeights=ReferenceDataExistingTuning.molecularWeights, SourceOfFragmentationPatterns=ReferenceDataExistingTuning.SourceOfFragmentationPatterns, sourceOfIonizationData=ReferenceDataExistingTuning.sourceOfIonizationData, relativeIonizationEfficiencies=ReferenceDataExistingTuning.relativeIonizationEfficiencies, moleculeIonizationType=ReferenceDataExistingTuning.moleculeIonizationType)
        ReferenceDataExistingTuningAfterCorrection.addSuffixToSourceOfFragmentationPatterns("_TuningCorrected")

        
        # #Now check if uncertainties already exist, and if they do then the two uncertainties need to be combined. Else, made equal.
        if G.calculateUncertaintiesInConcentrations == True: 
            if type(G.referenceFileUncertainties) != type(None):
                try:
                    ReferenceDataExistingTuningAfterCorrection.absolute_standard_uncertainties = (ReferenceDataExistingTuningAfterCorrection.absolute_standard_uncertainties**2 + ReferenceDataExistingTuningAfterCorrection.standardized_reference_patterns_tuning_uncertainties**2)**0.5
                except:
                    ReferenceDataExistingTuningAfterCorrection.absolute_standard_uncertainties = ReferenceDataExistingTuningAfterCorrection.standardized_reference_patterns_tuning_uncertainties
                ReferenceDataExistingTuningAfterCorrection.update_relative_standard_uncertainties()
                
        #TuningCorrector un-standardizes the patterns, so the patterns have to be standardized again.
        ReferenceDataExistingTuningAfterCorrection.standardized_reference_patterns=StandardizeReferencePattern(ReferenceDataExistingTuningAfterCorrection.standardized_reference_patterns,len(ReferenceDataExistingTuningAfterCorrection.molecules))
        
        if ReferenceDataExistingTuningAfterCorrection.referenceFileNameExtension == 'csv':
            ReferenceDataExistingTuningAfterCorrection.exportReferencePattern('ExportedReferencePatternExternalTuningCorrected.csv')
        if ReferenceDataExistingTuningAfterCorrection.referenceFileNameExtension == 'tsv':
            ReferenceDataExistingTuningAfterCorrection.exportReferencePattern('ExportedReferencePatternExternalTuningCorrected.tsv')
        
        ###The below code block is probably deprecated###
        # #Copying what is in GenerateReferenceDataList and ReferenceInputPreProcessing, we need to add the following block of code in if we want to allow uncertainties.   
        # if G.calculateUncertaintiesInConcentrations == True:
            # if type(G.referenceFileUncertainties) != type(None):
                # if type(G.referenceFileUncertainties) == type(float(5)) or  type(G.referenceFileUncertainties) == type(int(5)) :
                    # #TODO: Low priority. The below results in "nan" values. It could be better to change it to make zeros using a numpy "where" statement.
                    # G.referenceFileUncertainties = float(G.referenceFileUncertainties) #Maks sure we have a float.
                    # #Get what we need.
                    # provided_reference_patterns = ReferenceDataExistingTuningAfterCorrection.provided_reference_patterns
                    # provided_reference_patterns_without_masses = ReferenceDataExistingTuningAfterCorrection.provided_reference_patterns[:,1:] #[:,0] is mass fragments, so we slice to remove those. 
                    # #Make our variables ready.
                    # absolute_standard_uncertainties = provided_reference_patterns*1.0 #First we make a copy.
                    # absolute_standard_uncertainties_without_masses = (provided_reference_patterns_without_masses/provided_reference_patterns_without_masses)*G.referenceFileUncertainties #We do the division to get an array of ones. Then we multiply to get an array of the same size.
                    # #now we populate our copy's non-mass area.
                    # absolute_standard_uncertainties[:,1:] = absolute_standard_uncertainties_without_masses                                        
                    # ReferenceDataExistingTuningAfterCorrection.absolute_standard_uncertainties = absolute_standard_uncertainties
                    # #We can't convert to relative uncertainties yet because the file may not be standardized yet.
                # if type(G.referenceFileUncertainties) == type('string'):
                    # provided_reference_patterns_absolute_uncertainties, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form = readReferenceFile("ReferenceDataExistingTuningAfterCorrection.csv"[:-4]+"_absolute_uncertainties.csv","xyyy")
                    # TuningCorrectorGasMixtureCorrectedReferenceDataObject.absolute_standard_uncertainties = provided_reference_patterns_absolute_uncertainties #Just initializing the variable before filling it properly.
                    # maximum_absolute_intensities = numpy.amax(TuningCorrectorGasMixtureCorrectedReferenceDataObject.provided_reference_patterns[:,1:], axis = 0) #Find the maximum intensity for each molecule.
                    # TuningCorrectorGasMixtureCorrectedReferenceDataObject.absolute_standard_uncertainties[:,1:] = 100*TuningCorrectorGasMixtureCorrectedReferenceDataObject.absolute_standard_uncertainties[:,1:]/maximum_absolute_intensities
        ###The above code block is probably deprecated###
        
        #Now to make the mixed reference pattern using a function that extends one reference pattern by another.
        #In this loop the extension occurs for each item in ReferenceDataList.
        #The extendReferencePattern function uses existing add molecule and remove molecules functions.
        for ReferenceDataIndex in range(len(ReferenceDataList)):
            ReferenceDataList[ReferenceDataIndex], addedReferenceSlice = extendReferencePattern(ReferenceDataList[ReferenceDataIndex],ReferenceDataExistingTuningAfterCorrection)
            #extendReferencePattern also does add the uncertainties to the subobjects.
            
            #If specific molecules is on, we need to remove any molecules that were added that should not be added.
            if G.specificMolecules.lower() == "yes":
                moleculesToRemove = []
                for moleculeName in ReferenceDataList[ReferenceDataIndex].molecules:
                    if moleculeName not in G.chosenMoleculesNames:
                        moleculesToRemove.append(moleculeName)
                if len(moleculesToRemove) > 0:
                    ReferenceDataList[ReferenceDataIndex] = ReferenceDataList[ReferenceDataIndex].removeMolecules(moleculesToRemove)
            
            ReferenceDataList[ReferenceDataIndex].ExportAtEachStep = G.ExportAtEachStep
            if ReferenceDataList[ReferenceDataIndex].ExportAtEachStep == 'yes':
                ReferenceDataList[ReferenceDataIndex].ExportCollector('StandardizedReferencePattern', use_provided_reference_patterns=False)
            ReferenceDataList[ReferenceDataIndex].ClearZeroRowsFromStandardizedReferenceIntensities()
            if len(ReferenceDataList) == 1:
                if ReferenceDataList[ReferenceDataIndex].ExportAtEachStep == 'yes':
                    ReferenceDataList[ReferenceDataIndex].ExportCollector('StandardizedReferencePattern', use_provided_reference_patterns=False)
                if ReferenceDataList[ReferenceDataIndex].referenceFileNameExtension == 'csv':
                    ReferenceDataList[ReferenceDataIndex].exportReferencePattern('ExportedReferencePatternMixed.csv')
                if ReferenceDataList[ReferenceDataIndex].referenceFileNameExtension == 'tsv':
                    ReferenceDataList[ReferenceDataIndex].exportReferencePattern('ExportedReferencePatternMixed.tsv')
            else:
                if ReferenceDataList[ReferenceDataIndex].ExportAtEachStep == 'yes':
                    ReferenceDataList[ReferenceDataIndex].ExportCollector('StandardizedReferencePattern' + str(ReferenceDataIndex), use_provided_reference_patterns=False)
                if ReferenceDataList[ReferenceDataIndex].referenceFileNameExtension == 'csv':
                    ReferenceDataList[ReferenceDataIndex].exportReferencePattern("ExportedReferencePatternMixed" + str(ReferenceDataIndex) + ".csv")
                if ReferenceDataList[ReferenceDataIndex].referenceFileNameExtension == 'tsv':
                    ReferenceDataList[ReferenceDataIndex].exportReferencePattern("ExportedReferencePatternMixed" + str(ReferenceDataIndex) + ".tsv")
                    
        return ReferenceDataList

#this function eliminates (neglects) reference intensities that are below a certain threshold. Useful for solving 
#data that is giving negatives or over emphasizing small mass fragments,by assuming no contribution from the molecule at that mass fragment.
def ReferenceThresholdFilter(referenceDataArrayWithAbscissa,referenceValueThreshold):
    numMolecules = len(referenceDataArrayWithAbscissa[0]-1)
    if numMolecules > len(list(referenceValueThreshold)):
        referenceValueThreshold = list(referenceValueThreshold)* numMolecules
    referenceDataArray = referenceDataArrayWithAbscissa[:,1:] #all the data except the line of abscissa- mass fragment numbers
    for columncounter in range(len(referenceDataArray[0,:])):#goes through all columns in all rows in reference (this loop is one molecule at a time)
        for rowcounter in range(len(referenceDataArray[:,0])):#goes through all rows in references (one mass fragment at a time)        
            if len(referenceValueThreshold) == 1: #this is for if a single value was provided for referenceValueThreshold
                if referenceDataArray[rowcounter,columncounter] < referenceValueThreshold[0]:
                    referenceDataArray[rowcounter,columncounter] = 0 #made to be equal to zero
                # (len(referenceDataArray[:,0])) #this is masses.
                # (len(referenceDataArray[0,:])) #this is molecules
            else:  #this is for if values of referenceValueThreshold were provided for each molecule.
                if referenceDataArray[rowcounter,columncounter] < referenceValueThreshold[columncounter]: 
                    referenceDataArray[rowcounter,columncounter] = 0 #made to be equal to zero
    referenceDataArrayWithAbscissa[:,1:] = referenceDataArray #this puts changed referenceData back with mass fragment numbers
    return referenceDataArrayWithAbscissa
    
    
#The function just uses the input sheet to get values to change the reference sheet
#this is done by getting a list of all the values needed from the collected sheet at
#the right time and then gets each number following the first, and finds its ratio
#with the first, and multiplies that number by the number in the reference sheet in 
#order to change the second mass fragments number in the table
#TODO: add a check to pop out or warn people if someone chooses a mass that is not in experimental data.
def ExtractReferencePatternFromData(ExperimentData, referenceData, rpcChosenMolecules,rpcChosenMoleculesMF,rpcTimeRanges):
    copyOfreferenceData = copy.deepcopy(referenceData)    
    ExperimentData.mass_fragment_numbers = numpy.array(ExperimentData.mass_fragment_numbers, dtype=float)
    for moleculeIndex in range(len(rpcChosenMoleculesMF)):
        rpcChosenMoleculesMF[moleculeIndex] = numpy.array(rpcChosenMoleculesMF[moleculeIndex], dtype=float)
    for moleculeIndex in range(len(rpcChosenMoleculesMF)):
        currentMassesList = rpcChosenMoleculesMF[moleculeIndex]
        for massChoice in currentMassesList:
            if massChoice not in ExperimentData.mass_fragment_numbers:
                print("ERROR: You have selected to extract the reference pattern intensity of mass " + str(massChoice) + " from experimental data for one of your molecules, but this mass was not monitored. Please adjust your choices and run MSRESOLVE again."); sys.exit()
            if massChoice not in copyOfreferenceData.standardized_reference_patterns[:,0]:
                print("ERROR: You have selected to extract the reference pattern intensity of mass " + str(massChoice) + " from experimental data for one of your molecules, but this mass did not have a value greater than 0 in in your original reference pattern. Please adjust your inputs and run MSRESOLVE again."); sys.exit()
                
                
    try:
        type(copyOfreferenceData.provided_reference_patterns_absolute_uncertainties) #this checks if the variable exists.
    except: #create the variable if it does not exist.
        copyOfreferenceData.provided_reference_patterns_absolute_uncertainties = copyOfreferenceData.provided_reference_patterns*0.0
    for chosenmoleculescounter in range(len(rpcChosenMolecules)):#array-indexed for loop
        extractedIntensities = []
        allExtractedIntensities = []
        massfragindexerRef = [] #This is for the indices of the desired mass fragments in the reference pattern.
        for moleculecounter in range(len(copyOfreferenceData.molecules)):#array-indexed for loop
            if copyOfreferenceData.molecules[moleculecounter] == rpcChosenMolecules[chosenmoleculescounter]:#finds index of molecule
                if len(rpcChosenMoleculesMF[chosenmoleculescounter]) == 1:#if only one number is given then we alter the list so that the function will change (extract) all the other values scaled according to this one.
                    for x in range(len(ExperimentData.mass_fragment_numbers)):#array-indexed for loop
                        if ExperimentData.mass_fragment_numbers[x] != 0:#finds the other mass fragments and appends them. The !=0 is probably not necessary.
                            if ExperimentData.mass_fragment_numbers[x] != rpcChosenMoleculesMF[chosenmoleculescounter][0]:#We only append for masses different than the one we're scaling by.
                                rpcChosenMoleculesMF[chosenmoleculescounter].append(ExperimentData.mass_fragment_numbers[x])
                                rpcTimeRanges[chosenmoleculescounter].append(rpcTimeRanges[chosenmoleculescounter][0])
                #For a given molecule, and a given mass fragment to be changed, the below for loop finds the index of the fragment to be changed, with respect to the reference pattern's mass fragments.
                for eachChosenMoleculesMFIndex in range(len(rpcChosenMoleculesMF[chosenmoleculescounter])):#array-indexed for loop
                    for refMFCounter in range(len(copyOfreferenceData.standardized_reference_patterns[:,0])): #checks whole input list (or list that was made by previous loop)
                        if rpcChosenMoleculesMF[chosenmoleculescounter][eachChosenMoleculesMFIndex] == copyOfreferenceData.standardized_reference_patterns[:,0][refMFCounter]:#gets index of mass fragment number
                            massfragindexerRef.append(refMFCounter)
                for eachChosenMoleculesMFIndex in range(len(rpcChosenMoleculesMF[chosenmoleculescounter])):
                    for expMFCounter in range(len(ExperimentData.mass_fragment_numbers)):
                        if rpcChosenMoleculesMF[chosenmoleculescounter][eachChosenMoleculesMFIndex] == ExperimentData.mass_fragment_numbers[expMFCounter]: #This if statment means that the masss fragment number that is desired for this molecule is at the current expMFCounter index.
                            for timecounter in range(len(ExperimentData.times)):#array-indexed for loop. This for loop extracts the data.
                                if (rpcTimeRanges[chosenmoleculescounter][0] <= ExperimentData.times[timecounter]) and (rpcTimeRanges[chosenmoleculescounter][1] >= ExperimentData.times[timecounter]):#gets index of time
                                    extractedIntensities.append(ExperimentData.workingData[timecounter,expMFCounter])
                            #Place current extractedIntensities in a larger list using copy so it can be cleared without affecting allExtractedIntensities
                            allExtractedIntensities.append(copy.copy(extractedIntensities))
                            #clear extracted intensities for the next mass fragment
                            extractedIntensities.clear()
                #Convert list to an array
                allExtractedIntensitiesArray = numpy.array(allExtractedIntensities)
                #Initialize empty list to store average values
                allExtractedIntensitiesAveragesList = []
                allExtractedIntensitiesStandardDeviationsList = []
                #For loop to find the average intensity values
                #Then store the average value in the allExtractedIntensityAverage list
                for eachChosenMoleculesMFIndex in range(len(allExtractedIntensitiesArray)):
                    intensitiesAverage = numpy.average(allExtractedIntensitiesArray[eachChosenMoleculesMFIndex])
                    intensitiesStandardDeviation = numpy.std(allExtractedIntensitiesArray[eachChosenMoleculesMFIndex])
                    allExtractedIntensitiesAveragesList.append(intensitiesAverage)
                    allExtractedIntensitiesStandardDeviationsList.append(intensitiesStandardDeviation)
                #For loop to overwrite a chosen mass fragment's signal in the reference file with the product of the extracted ratios and the reference signal of the base mass fragment (that is, to make a reference pattern with a ratio matching the extracted ratio)
                normalizationFactor = copyOfreferenceData.provided_reference_patterns[massfragindexerRef[0],moleculecounter+1] #The +1 is to skip the column that has mass fragment values.
                if normalizationFactor == 0: #This is just for the case that there's an 'exception' with no reasonable value provided by the line above.
                    normalizationFactor = 1
                    print("WARNING: There may be an error occurring when extracting the reference pattern. Search for normalizationFactor == 0 in the code for debugging.")
                for eachChosenMoleculesMFIndex in range(len(rpcChosenMoleculesMF[chosenmoleculescounter])): #The +1 below is b/c first column is mass frag
                    copyOfreferenceData.provided_reference_patterns[massfragindexerRef[eachChosenMoleculesMFIndex],moleculecounter+1] = (allExtractedIntensitiesAveragesList[eachChosenMoleculesMFIndex]/allExtractedIntensitiesAveragesList[0])*normalizationFactor
                    copyOfreferenceData.provided_reference_patterns_absolute_uncertainties[massfragindexerRef[eachChosenMoleculesMFIndex],moleculecounter+1] = (allExtractedIntensitiesStandardDeviationsList[eachChosenMoleculesMFIndex]/allExtractedIntensitiesAveragesList[0])*normalizationFactor #Need to divide by the same thing and multiply by the same thing as previous line to get right scaling.
    return copyOfreferenceData.provided_reference_patterns, copyOfreferenceData.provided_reference_patterns_absolute_uncertainties

'''
RemoveUnreferencedMasses() is used to prune ExperimentData.workingData and ExperimentData.mass_fragment_numbers 
in accordance with the available reference data. If there is no reference data for a particular mass present
in the ExperimentData.workingData then delete that mass column from ExperimentData.workingData. Also delete that 
mass from the list ExperimentData.mass_fragment_numbers
'''
def RemoveUnreferencedMasses(ExperimentData, reference):  ## DEPRECATED Replaced by KeepOnlySelectedYYYYColumns() Dec 2017
    # masses available in the reference data
    referenceDataMassAbscissa = reference[:,0]

    # convenience variables
    workingData = ExperimentData.workingData
    massFragmentNumbers = ExperimentData.mass_fragment_numbers

    # Storage for deletion (will cause trouble if we delete as we iterate)
    deletion_indices = []

    # loop through the massFragmentNumbers, if they aren't in the
    # reference data get rid of them from both massFragmentNumbers
    # and get rid of that mass column in workingData
    for massFragIndex, massFragment in enumerate(massFragmentNumbers):
        if massFragment not in referenceDataMassAbscissa:
            deletion_indices.append(massFragIndex)

    # Now delete the unreferenced masses from massFragmentNumbers and workingData
    massFragmentNumbers = numpy.delete(massFragmentNumbers.astype(int), deletion_indices)
    workingData = numpy.delete(workingData, deletion_indices, 1)

    return massFragmentNumbers, workingData


'''
MassFragChooser() compares chosenMassFragments and ExperimentData.mass_fragment_numbers. If there is a 
mass in mass_fragment_numbers that isn't in chosenMassFragments it removes that mass from
mass_fragment_numbers. It also removes the column of data that corresponds to that mass
from ExperimentData.workingData.
NOTE: This code was changed algorithmically by Clint on 171204. We believe that it is functioning
properly due to the similar outputs prior to and after the change. It has not, however, been
rigorously tested and the old code is still below in the function body but commented out.
'''
def MassFragChooser (ExperimentData, chosenMassFragments):    ## DEPRECATED Replaced by KeepOnlySelectedYYYYColumns() Dec 2017

    # convenience variables
    workingData = ExperimentData.workingData
    massFragmentNumbers = ExperimentData.mass_fragment_numbers

    # Storage for deletion (will cause trouble if we delete as we iterate)
    deletion_indices = []

    # loop through massFragmentNumbers, if there is a mass in the list
    # that isnt in choosenMassFragments then delete it from massFragmentNumbers.
    # Further delete the column that corresponds to that mass from workingData
    for massFragIndex, massFragment in enumerate(massFragmentNumbers):
        if massFragment not in chosenMassFragments:
            deletion_indices.append(massFragIndex)

    # Now remove the unreferened data from massFragmentNumbers and workingData
    massFragmentNumbers = numpy.delete(massFragmentNumbers.astype(int), deletion_indices)
    workingData = numpy.delete(workingData, deletion_indices, 1)

    return massFragmentNumbers,workingData


    
    ## OLD VERSION: Replaced by code directly above (in MassFragChoosed)
    # mass_fragment_length = len(ExperimentData.mass_fragment_numbers) 
    # mass_fragment_numbers_holder = []
    # if specificMassFragments == 'yes': #from data edit file
    #     for masscounter in range(mass_fragment_length):#array-indexed for loop
    #         for mfcounter in range(len(chosenMassFragments)):#array-indexed for loop
    #             if chosenMassFragments[mfcounter] == ExperimentData.mass_fragment_numbers[masscounter]: #gets an index
    #                 mass_fragment_numbers_holder.append(chosenMassFragments[mfcounter])
    #     mass_fragment_numbers2 = mass_fragment_numbers_holder 
    #     place_holder = 0#keeps array index with deletion
    #     for massFragmentIndex in range(len(ExperimentData.mass_fragment_numbers)):#this loop deletes any collected columns that belong to the mass fragment numbers that got deleted
    #         summer = 0#resets sum each loop
    #         for mFIndex in range(len(mass_fragment_numbers2)):
    #             if ExperimentData.mass_fragment_numbers[massFragmentIndex] == mass_fragment_numbers2[mFIndex]:#if there is a matching number in the mass fragments chosen then it is remembered
    #                 summer = summer+1
    #             if mFIndex == len(mass_fragment_numbers2)-1:#at the end of the loop
    #                 if summer == 0:#if none of the numbers were equal
    #                     ExperimentData.workingData = numpy.delete(ExperimentData.workingData,(massFragmentIndex-place_holder),axis = 1)
    #                     place_holder = place_holder + 1

                        
    ## OLD BELOW HERE: this functionality is now in a seperate function above 'RemoveUnreferencedMasses'
    
    # else:#the checks later need to check against the reference data to delete any not present masses
    #     mass_fragment_numbers2 = ExperimentData.mass_fragment_numbers
    # referenceabscissa = reference[:,0]#getting the abscissa
    # place_holder = 0#saving the index
    # for mFIndex2 in range(len(mass_fragment_numbers2)):#array indexed loop
    #     summer = 0#each loop the sum resets
    #     for refAbscissaIndex in range(len(referenceabscissa)):#checks all of abscissa array for the current mass fragment value
    #         if mass_fragment_numbers2[mFIndex2-place_holder] == referenceabscissa[refAbscissaIndex]:#If there is one equal, the summer becomes one
    #             summer = summer+1
    #         if refAbscissaIndex == len(referenceabscissa)-1:#At the end of the loop
    #             if summer == 0:#if this value is not present, it is deleted from the mass fragment numbers list and 
    #                 mass_fragment_numbers2 = numpy.delete(mass_fragment_numbers2,(mFIndex2-place_holder))
    #                 ExperimentData.workingData = numpy.delete(ExperimentData.workingData,(mFIndex2-place_holder),axis = 1)
    #                 place_holder = place_holder + 1
                    
    #return [mass_fragment_numbers2,ExperimentData.workingData]
    
#This function "combines" two reference patterns by adding to the first pattern any molecules that are only the 2nd pattern.
#This is *not* a merging function. If the molecule is only in the first pattern, then it keeps the version in the first pattern.
def extendReferencePattern(OriginalReferenceData, ReferenceDataToExtendBy):
    #First figure out which molecules are needed and which are duplicates.
    moleculesToExtendBy = [] #first initialize, then we'll populate below
    duplicateMolecules = [] #first initialize, then populate below
    for moleculeName in ReferenceDataToExtendBy.molecules:
        if moleculeName in OriginalReferenceData.molecules:
            duplicateMolecules.append(moleculeName)
        else: #implies not in OriginalReferenceData.molecules)
            moleculesToExtendBy.append(moleculeName)
    if len(moleculesToExtendBy) == 0:
        return OriginalReferenceData, None
    ReferenceDataToExtendBy = ReferenceDataToExtendBy.removeMolecules(duplicateMolecules) #this removes the duplicates and makes a new object. It does not affect the original reference data object.
    #we will assume that the pattern that we're going to use is the standardized_reference_patterns rather than the provided_reference_patterns
    if hasattr(ReferenceDataToExtendBy, 'standardized_reference_patterns'):
        patternToExtendBy = ReferenceDataToExtendBy.standardized_reference_patterns
    else:
        patternToExtendBy = ReferenceDataToExtendBy.provided_reference_patterns
    extendedReferenceData = copy.deepcopy(OriginalReferenceData)  #We need to make a copy so we don't change anything in the original reference data object, since that may not be desired.
    #The below command will modify the  extendedReferenceData, so we don't need to return anything from the function though we will do so.
    extendedReferenceData = extendedReferenceData.addMolecules(provided_reference_patterns=ReferenceDataToExtendBy.standardized_reference_patterns, electronnumbers=ReferenceDataToExtendBy.electronnumbers, molecules=ReferenceDataToExtendBy.molecules, molecularWeights=ReferenceDataToExtendBy.molecularWeights, SourceOfFragmentationPatterns=ReferenceDataToExtendBy.SourceOfFragmentationPatterns, sourceOfIonizationData=ReferenceDataToExtendBy.sourceOfIonizationData, relativeIonizationEfficiencies=ReferenceDataToExtendBy.relativeIonizationEfficiencies, moleculeIonizationType=ReferenceDataToExtendBy.moleculeIonizationType)
    ##TODO: Currently, the addMolecules feature does not extend the absolute_standard_uncertainties. This should probably be changed with an additoinal argument. For now, we extend the absolute_standard_uncertainties in a separate line.
    #print("line 1114", extendedReferenceData.absolute_standard_uncertainties, ReferenceDataToExtendBy.absolute_standard_uncertainties)
    if hasattr(extendedReferenceData, "absolute_standard_uncertainties"):
        extendedReferenceData.absolute_standard_uncertainties = DataFunctions.addXYYYtoXYYY(extendedReferenceData.absolute_standard_uncertainties, ReferenceDataToExtendBy.absolute_standard_uncertainties)
        extendedReferenceData.update_relative_standard_uncertainties()
    #print("line 1115", extendedReferenceData.absolute_standard_uncertainties); sys.exit()
    addedReferenceSlice = ReferenceDataToExtendBy #the addition has finished.
    extendedReferenceData.populateIonizationEfficiencies(extendedReferenceData.AllMID_ObjectsDict)
    extendedReferenceData.exportIonizationInfo()
    return extendedReferenceData, addedReferenceSlice


#desiredMoleculesOrder should be a list of strings of molecule names.
def rearrangeReferenceData(ReferenceData, desiredMoleculesOrder):
    #We are going to make separte ReferenceData objects for each molecule, and then 'concatenate' those objects to make the rearranged molecular object.
    #This function can thus also be used for truncating a ReferenceData object.
    #We will first make a list where we will place ReferenceData object copy for each of the desiredMolecules.
    ReferenceData = copy.deepcopy(ReferenceData)
    individualMoleculeReferenceDataList = []
    for moleculeIndex,moleculeName in enumerate(desiredMoleculesOrder):
        #we are going to make ReferenceData object with only moleculeName.
        #make a list of the other molecules so we can remove them.
        otherMolecules = copy.deepcopy(list(desiredMoleculesOrder)) #initializing.
        otherMolecules.remove(moleculeName) #remove modifies the original ilst.
        #Below creates the individualMoleculeReferenceData and then appends it with the same moleculeIndex to the individualMoleculeReferenceDataList
        individualMoleculeReferenceData = ReferenceData.removeMolecules(otherMolecules)
        individualMoleculeReferenceDataList.append(individualMoleculeReferenceData)
    #Now having a list of all the desired molecules individual ReferenceData objects, we just need to put them together.
    combinedReferenceData  = individualMoleculeReferenceDataList[0] #initialize with the first molecule.
    for individualMoleculeReferenceData in individualMoleculeReferenceDataList[1:]:#note that we are skipping index 0.
        combinedReferenceData, addedReferenceSlice = extendReferencePattern(OriginalReferenceData= combinedReferenceData, ReferenceDataToExtendBy= individualMoleculeReferenceData)
    combinedReferenceData.populateIonizationEfficiencies(combinedReferenceData.AllMID_ObjectsDict)
    return combinedReferenceData
    
#This function operates in a parallel way to trimDataMasses, but it operates on the reference data and all of it's constituent variables  
def trimDataMoleculesToMatchChosenMolecules(ReferenceData, chosenMolecules):
    
    print("MoleculeChooser")
    #getting a list of all molecules (a header) to compare to during trimming.	    
    allMoleculesList = ReferenceData.molecules
    
    #initializing object that will become the trimmed copy of ReferenceData
    trimmedRefererenceData = copy.deepcopy(ReferenceData)
    #trim the reference fragmentation patterns to only the selected molecules 
    
    #unused trimmed copy molecules is just a place holder to dispose of a function return that is not needed
    trimmedReferenceIntensities, trimmedMoleculesList = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.provided_reference_patterns[:,1:],
                                                                                                                    allMoleculesList, chosenMolecules, header_dtype_casting=str)  
    #add a second dimension to the reference data
    trimmedReferenceMF = numpy.reshape(trimmedRefererenceData.standardized_reference_patterns[:,0],(-1,1))
    #Add the abscissa back into the reference values
    trimmedRefererenceData.provided_reference_patterns = numpy.hstack((trimmedReferenceMF,trimmedReferenceIntensities))
    
    #repeat the process for the standardized_reference_patterns
    if hasattr(trimmedRefererenceData, 'standardized_reference_patterns'):
        standardized_reference_patterns_mass_fragments = trimmedRefererenceData.standardized_reference_patterns[:,0] 
        trimmedReferenceStandardizedIntensities, trimmedMoleculesList = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.standardized_reference_patterns[:,1:],
                                                                                                                        allMoleculesList, chosenMolecules, header_dtype_casting=str) 
        trimmedReferenceMF = numpy.reshape(standardized_reference_patterns_mass_fragments,(-1,1))
        trimmedRefererenceData.standardized_reference_patterns =numpy.hstack((trimmedReferenceMF,trimmedReferenceStandardizedIntensities))
    
    #TODO: The below line works with provided_reference_patterns. This is because trimDataMoleculesToMatchChosenMolecules
    #TODO continued: is currently working prior to standardized Reference patterns existing, and also because it is occurring
    #TODO continued: Before we have the monitored mass fragments (which also occurs later data analysis).
    #TODO continued: The best solution is probably to do the standardization earlier and then do this trimming after that.

    #Shorten the electronnumbers to the correct values, using the full copy of molecules. Do the same for molecularWeights and sourceInfo
    trimmedRefererenceData.electronnumbers, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.electronnumbers, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    trimmedRefererenceData.molecularWeights, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.molecularWeights, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    trimmedRefererenceData.moleculeIonizationType, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.moleculeIonizationType, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    trimmedRefererenceData.relativeIonizationEfficiencies, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.relativeIonizationEfficiencies, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    trimmedRefererenceData.SourceOfFragmentationPatterns, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.SourceOfFragmentationPatterns, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    trimmedRefererenceData.sourceOfIonizationData, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.sourceOfIonizationData, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    if G.calculateUncertaintiesInConcentrations == True: 
        if type(G.referenceFileUncertainties) != type(None): #The [:,1:] is related to the mass fragements.  Just copying the two lines syntax from above.
            trimmedAbsoluteUncertainties, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.absolute_standard_uncertainties[:,1:], 
                                                                                                                    allMoleculesList, chosenMolecules, header_dtype_casting=str)  
        #add a second dimension as above.
        trimmedRefererenceData.absolute_standard_uncertainties = numpy.hstack((trimmedReferenceMF,trimmedAbsoluteUncertainties))

    trimmedRefererenceData.molecules = trimmedMoleculesList
    trimmedRefererenceData.ExportCollector("MoleculeChooser", use_provided_reference_patterns=True)    
    
    #remove any zero rows that may have been created #this should not be here. Zero rows can be cleared outside of this function. There are applications where this removing of zero rows causes problems.
    #trimmedRefererenceData.ClearZeroRowsFromProvidedReferenceIntensities()
    if hasattr(trimmedRefererenceData, 'standardized_reference_patterns'):
        trimmedRefererenceData.standardized_reference_patterns=StandardizeReferencePattern(trimmedRefererenceData.standardized_reference_patterns,len(trimmedRefererenceData.molecules))
        trimmedRefererenceData.ExportCollector("MoleculeChooser", use_provided_reference_patterns=False)    
    trimmedRefererenceData.populateIonizationEfficiencies(trimmedRefererenceData.AllMID_ObjectsDict)
    return trimmedRefererenceData
    
'''
trimDataMassesToMatchChosenMassFragments() and trimDataMassesToMatchReference() are just a wrapper functions for calls to DataFunctions.KeepOnlySelectedYYYYColumns(). 
Both of the functions trim ExperimentData.workingData and ExperimentData.mass_fragment_numbers. 
The first function trims the data according to the mass fragment selections in G.chosenMassFragments.
The second function trims the data to remove any mass fragments for which there is no ReferenceData. 

Parameters:
ExperimentData - of type MSData, the one instantiated in main() named ExperimentData is a good example of one
    that will work here
ReferenceData - of type MSReference, ReferenceData from main() is a good example
chosenMassFragments  - list of integers, like the one created in UserInput  
'''
def trimDataMassesToMatchChosenMassFragments(ExperimentData, chosenMassFragments):
    # If we are only interested in a subset of the MS data
    # and that subset is a subset of the loaded data
    # remove the irrelevant mass data series from ExperimentData.mass_fragment_numbers
    # and the corresponding colums from ExperimentData.workingData
    trimmedExperimentData = copy.deepcopy(ExperimentData)
    trimmedExperimentData2 = copy.deepcopy(ExperimentData) 
    
    print("MassFragChooser")
    (trimmedExperimentData.workingData, trimmedExperimentData.mass_fragment_numbers) = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedExperimentData.workingData,
                                                                                                            trimmedExperimentData.mass_fragment_numbers,
                                                                                                            chosenMassFragments, header_dtype_casting=float)
    if G.calculateUncertaintiesInConcentrations:
        if type(G.collectedFileUncertainties) != type(None): #TODO: As of Feb 5th 2020, this if statement has been added but not tested. It simply mimics the above code.   The below function   trimDataMassesToMatchReference is similar and has been tested.
                (trimmedExperimentData2.rawsignals_absolute_uncertainties, trimmedExperimentData2.mass_fragment_numbers) = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedExperimentData2.rawsignals_absolute_uncertainties,
                                                                                                                trimmedExperimentData2.mass_fragment_numbers,
                                                                                                                chosenMassFragments, header_dtype_casting=float)
                trimmedExperimentData.rawsignals_absolute_uncertainties = trimmedExperimentData2.rawsignals_absolute_uncertainties
    trimmedExperimentData.ExportCollector("MassFragChooser")
    
    return trimmedExperimentData

def trimDataMassesToMatchReference(ExperimentData, ReferenceData):
    
    trimmedExperimentData = copy.deepcopy(ExperimentData)
    trimmedExperimentData2 = copy.deepcopy(ExperimentData) #needed to make a separate temporary data structure because it gets modified during DataFunctions.
    
    # Remove elements of ExperimentData.mass_fragment_numbers for which there is no matching mass in the reference data.
    # Also remove the corresponding mass data column from Experiment.workingData.
    (trimmedExperimentData.workingData, trimmedExperimentData.mass_fragment_numbers) = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedExperimentData.workingData,
                                                                                                        trimmedExperimentData.mass_fragment_numbers,
                                                                                                        ReferenceData.standardized_reference_patterns[:,0], header_dtype_casting=float)
    if G.calculateUncertaintiesInConcentrations:
        if (type(G.collectedFileUncertainties) != type(None)) and (str(G.collectedFileUncertainties).lower() != str('None').lower()): #Note: As of Feb 5th 2020, this code has been added and tested.
            if hasattr(trimmedExperimentData2, 'rawsignals_absolute_uncertainties'):
                (trimmedExperimentData2.rawsignals_absolute_uncertainties, trimmedExperimentData2.mass_fragment_numbers) = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedExperimentData2.rawsignals_absolute_uncertainties,
                                                                                                                    trimmedExperimentData2.mass_fragment_numbers,
                                                                                                                    ReferenceData.standardized_reference_patterns[:,0], header_dtype_casting=float)
                trimmedExperimentData.rawsignals_absolute_uncertainties = trimmedExperimentData2.rawsignals_absolute_uncertainties
                                                                                                            
    return trimmedExperimentData
    
#The correction values are calculated based on each molecule's molecular weight, number of electrons and mass fragmentation pattern
#This obtains the ionization energy, transmission and electron multiplier gains from the mw, e- and the application of simple equations
#Then it gets the sum of each molecules frag pattern and uses all these values to get each moleule's correction value via the
#system created by Madix and Ko and puts this in an answers array
def CorrectionValuesObtain(ReferenceData):
    ReferenceData.ExportCollector("ReferencePatternOriginalForCorrectionValues",use_provided_reference_patterns=False)
    if len(G.referenceFileStandardTuningAndForm) == 0: #This is the simple case, if there is no tuningCorrectionIntensity feature being used.
        ReferenceDataForCorrectionValues = copy.deepcopy(ReferenceData)
    if len(G.referenceFileStandardTuningAndForm) > 0: #if referenceFileStandardTuningAndForm is provided, that means that the tuningCorrectorIntensity feature will be used.
        referenceFileStandardTuningAndForm = G.referenceFileStandardTuningAndForm
        if ReferenceData.referenceFileNameExtension == 'csv':
            ReferenceData.exportReferencePattern("ExportedReferencePatternOriginalForCorrectionValues.csv") #Separate file, outside of export collector.
        if ReferenceData.referenceFileNameExtension == 'tsv':
            ReferenceData.exportReferencePattern("ExportedReferencePatternOriginalForCorrectionValues.tsv") #Separate file, outside of export collector.
        #Below, we are going to get the tuningCorrector coefficients for Desired towards Standard tuning. We need to **intentionally** reverse the arguments below, because are going to be making a mixed reference pattern which has some molecules from the Standard Tuning and also some where where the Original Tuning is changed to Standard Tuning.
        
        #NEED AN IF STATEMENT TO DETERMINE THESE COEFFICENTS FROM GAS MIXTURE IF USING THAT.
        if str(G.referenceFileStandardTuningAndForm[0]).lower() != 'gaxmixture':
            if len(G.referenceFileDesiredTuningAndForm) == 0: #Use the original reference pattern if blank.
                referenceFileDesiredTuningAndForm = [G.referenceFileNamesList[0],  G.referenceFormsList[0]] #Take the first item from G.referenceFileNamesList and from G.referenceFormsList.
            else:
                referenceFileDesiredTuningAndForm = G.referenceFileDesiredTuningAndForm
            abcCoefficients, abcCoefficients_covmat = ABCDetermination(ReferencePatternExistingTuning_FileNameAndForm=referenceFileDesiredTuningAndForm, ReferencePatternDesiredTuning_FileNameAndForm = G.referenceFileStandardTuningAndForm, exportCoefficients=False) #We will separately export the coefficents for this usage.
        #The above abcCoefficients and abcCoefficients_covmat are for the direction that if we multiply the existing tuning we will get the standard tuning. This factor is actually the same as what we want for our tuning factor intensity correction: if the "existing" pattern is low, that means that the actual amount of the ion is "higher".  The second question is whether the correctionValues should be multiplied by or divided by this factor.  In general, the equation is that "signals*correctionValue = concentration".  So to bring a small signal up, this is already in the right direction.
        #We will also use this to create a mixed pattern for getting the standard tuning correction values.
        #Export the TuningCorrectorCoefficients as a list.
        with open('ExportedTuningCorrectorIntensityCoefficients.txt', 'w') as the_file:
            the_file.write(str(list(abcCoefficients))) 
        #Export the TuningCorrectorCoefficients_cov to a csv.
        numpy.savetxt('ExportedTuningCorrectorIntensityCoefficients_covmat.csv', abcCoefficients_covmat, delimiter=",")
        
        #We are going to make a standard tuning type reference pattern in ReferenceDataOriginalStandardTuning.
        ##Start block of making mixed tuning##
        listOfNeededMolecules = ReferenceData.molecules #intentionally grab molecule list from the original data object fed in.
        ReferenceDataOriginalStandardTuning = copy.deepcopy(ReferenceData)  #This is mainly needed for G.referenceFileStandardTuningAndForm
        ReferenceDataStandardTuning = createReferenceDataObject ( G.referenceFileStandardTuningAndForm[0], G.referenceFileStandardTuningAndForm[1], AllMID_ObjectsDict=G.AllMID_ObjectsDict)
        #TODO: Make the below block of code a function. This block is actually the copy of a block from tuning correction.
        if (G.calculateUncertaintiesInConcentrations == True) and (type(G.referenceFileUncertainties) != type(None)):
            if type(G.referenceFileUncertainties) == type(float(5)) or  type(G.referenceFileUncertainties) == type(int(5)) :
                #TODO: Low priority. The below results in "nan" values. It could be better to change it to make zeros using a numpy "where" statement.
                G.referenceFileUncertainties = float(G.referenceFileUncertainties) #Make sure we have a float.
                #Get what we need.
                provided_reference_patterns = ReferenceDataStandardTuning.provided_reference_patterns
                provided_reference_patterns_without_masses = ReferenceDataStandardTuning.provided_reference_patterns[:,1:] #[:,0] is mass fragments, so we slice to remove those. 
                #Make our variables ready.
                absolute_standard_uncertainties = provided_reference_patterns*1.0 #First we make a copy.
                absolute_standard_uncertainties_without_masses = (provided_reference_patterns_without_masses/provided_reference_patterns_without_masses)*G.referenceFileUncertainties #We do the division to get an array of ones. Then we multiply to get an array of the same size.
                #now we populate our copy's non-mass area.
                absolute_standard_uncertainties[:,1:] = absolute_standard_uncertainties_without_masses                                        
                ReferenceDataStandardTuning.absolute_standard_uncertainties = absolute_standard_uncertainties
                #We can't convert to relative uncertainties yet because the file may not be standardized yet.
            else:
                print("WARNING: Line 1325 of MSRESOLVE.py else statement has not been programmed. ReferenceDataStandardTuning currently only receives uncertainties if there is an integer in the referenceFileUncertainties.")
                ReferenceDataStandardTuning.absolute_standard_uncertainties = ReferenceDataStandardTuning.provided_reference_patterns*1.0 #First we make a copy.
                ReferenceDataStandardTuning.absolute_standard_uncertainties[:,1:] = 0
        listOfStandardTuningMoleculePatternsAvailable = copy.deepcopy(ReferenceDataStandardTuning.molecules)
        #get a list of molecules to remove from the ReferenceDataStandardTuning file.
        listOfMoleculesToRemove = []
        for moleculeName in listOfStandardTuningMoleculePatternsAvailable:
            if moleculeName not in listOfNeededMolecules:
                listOfMoleculesToRemove.append(moleculeName)
        ReferenceDataStandardTuning = ReferenceDataStandardTuning.removeMolecules(listOfMoleculesToRemove)
        #Now have to convert the original reference data for anlaysis to standard tuning to make the mixed pattern.
        referenceDataArrayWithAbscissa, referenceDataArrayWithAbscissa_tuning_uncertainties = TuningCorrector(ReferenceDataOriginalStandardTuning.standardized_reference_patterns,abcCoefficients, abcCoefficients_covmat, referenceFileExistingTuningAndForm=None,referenceFileDesiredTuningAndForm=None,measuredReferenceYorN="no")
        #TuningCorrector un-standardizes the patterns, so the patterns have to be standardized again.
        ReferenceDataOriginalStandardTuning.standardized_reference_patterns = StandardizeReferencePattern(referenceDataArrayWithAbscissa)
        ReferenceDataOriginalStandardTuning.standardized_reference_patterns_tuning_uncertainties = referenceDataArrayWithAbscissa_tuning_uncertainties
        #add "TuningCorrected" to SourceOfFragmentationPatterns
        ReferenceDataOriginalStandardTuning.SourceOfFragmentationPatterns = list(ReferenceDataOriginalStandardTuning.SourceOfFragmentationPatterns)
        for sourceIndex in range(len(ReferenceDataOriginalStandardTuning.SourceOfFragmentationPatterns)):
            ReferenceDataOriginalStandardTuning.SourceOfFragmentationPatterns[sourceIndex] = ReferenceDataOriginalStandardTuning.SourceOfFragmentationPatterns[sourceIndex] + '-TuningCorrected'
        #Now check if uncertainties already exist, and if they do then the two uncertainties need to be combined. Else, made equal.
        #We don't export these uncertainties since this is just an intermediate step.
        try:
            ReferenceDataOriginalStandardTuning.absolute_standard_uncertainties[:,1:] = (ReferenceDataOriginalStandardTuning.absolute_standard_uncertainties[:,1:]**2 + ReferenceDataOriginalStandardTuning.standardized_reference_patterns_tuning_uncertainties[:,1:]**2)**0.5
        except:
            ReferenceDataOriginalStandardTuning.absolute_standard_uncertainties = ReferenceDataOriginalStandardTuning.standardized_reference_patterns_tuning_uncertainties
        #Export (won't be needed).
        if ReferenceDataOriginalStandardTuning.referenceFileNameExtension == ".csv":
            ReferenceDataOriginalStandardTuning.exportReferencePattern("ExportedReferencePatternOriginalForCorrectionValuesStandardTuning.csv")
        if ReferenceDataOriginalStandardTuning.referenceFileNameExtension == ".tsv":
            ReferenceDataOriginalStandardTuning.exportReferencePattern("ExportedReferencePatternOriginalForCorrectionValuesStandardTuning.tsv")
        #Now create the mixed reference pattern.  The pattern for ReferenceDataStandardTuning will be extended by ReferenceDataOriginalStandardTuning.
        ReferenceDataStandardTuning, addedReferenceSlice = extendReferencePattern(ReferenceDataStandardTuning, ReferenceDataOriginalStandardTuning)
        #This is the ReferenceDataObject to carry forward.
        if ReferenceDataStandardTuning.referenceFileNameExtension == ".csv":
            ReferenceDataStandardTuning.exportReferencePattern("ExportedReferencePatternStandardForCorrectionValuesMixedStandardTuning.csv")
        if ReferenceDataStandardTuning.referenceFileNameExtension == ".tsv":
            ReferenceDataStandardTuning.exportReferencePattern("ExportedReferencePatternStandardForCorrectionValuesMixedStandardTuning.tsv")
        if G.createMixedTuningPattern == False:
            listOfMoleculesToRemove = []
            for moleculeName in ReferenceDataStandardTuning.molecules:
                if moleculeName not in ReferenceData.molecules:
                    listOfMoleculesToRemove.append(moleculeName)
            ReferenceDataStandardTuning = ReferenceDataStandardTuning.removeMolecules(listOfMoleculesToRemove)
        #rearrange the molecules to the right order.
        ReferenceDataStandardTuning = rearrangeReferenceData(ReferenceData=ReferenceDataStandardTuning, desiredMoleculesOrder=ReferenceData.molecules)
        if ReferenceDataStandardTuning.referenceFileNameExtension == 'csv':
            ReferenceDataStandardTuning.exportReferencePattern("ExportedReferencePatternStandardForCorrectionValuesMixedStandardTuningRearranged.csv")
        if ReferenceDataStandardTuning.referenceFileNameExtension == 'tsv':
            ReferenceDataStandardTuning.exportReferencePattern("ExportedReferencePatternStandardForCorrectionValuesMixedStandardTuningRearranged.tsv")        
        #move the pointer.
        ReferenceDataForCorrectionValues = ReferenceDataStandardTuning
    #Before the below, ReferenceDataForCorrectionValues and ReferenceData must be made the same size and width. To do that, we will extend each to have any rows missing that the other has.
    ReferenceDataMassFragments = ReferenceData.standardized_reference_patterns[:,0]
    ReferenceDataForCorrectionValuesMassFragments = ReferenceDataForCorrectionValues.standardized_reference_patterns[:,0]
    ReferenceData.extendMassFragments(ReferenceDataForCorrectionValuesMassFragments)
    ReferenceDataForCorrectionValues.extendMassFragments(ReferenceDataMassFragments)   
    reference_width = len(ReferenceDataForCorrectionValues.standardized_reference_patterns[0,:])  #This is number of molecules plus 1 because of the mass fragments column.
    reference_height = len(ReferenceDataForCorrectionValues.standardized_reference_patterns[:,0]) #this is the number of mass fragments.
    correction_values_direct = ReferenceDataForCorrectionValues.standardized_reference_patterns*1.0 #just initializing as same size, but note that this has the mass fragments column.
    correction_values_direct_relative_uncertainty = correction_values_direct*1.0
    if (G.calculateUncertaintiesInConcentrations == True) and (type(G.referenceFileUncertainties) != type(None)):  #This means we'll need uncertainties on each correction value.
        ReferenceDataForCorrectionValues.relative_correction_uncertainties = ReferenceDataForCorrectionValues.relative_standard_uncertainties[:,1:]*0.0 #make an array with the same shape, without the mass fragments.
    else: 
        correction_values_relative_uncertainties = None #If uncertainties are not going to be calculated, we'll use none here and just return that.
    #If we are using the tuningCorrectorIntensity feature, then we would want to calculate uncertainties when they are desired.
    if (G.calculateUncertaintiesInConcentrations == True) and len(G.referenceFileStandardTuningAndForm) > 0: 
        correction_values_relative_uncertainties = True
    
    #the first for loop here gets all of the values for e- and mw and uses them to get the
    #respective values that can find the correction factor for each mass fragment of each molecule
    for column_counter in range(1,reference_width): #array-indexed for loop, skips column one b/c that is the mass fragment numbers, not relative intensities. This loops across each molecule.
        moleculeName = ReferenceDataForCorrectionValues.molecules[column_counter-1]
        ionization_efficiency = ReferenceDataForCorrectionValues.relativeIonizationEfficiencies[column_counter-1]
        current_molecule_correction_factors_list = []
        quotients = numpy.zeros(len(ReferenceDataForCorrectionValues.standardized_reference_patterns[:,0])) #This means the same length as the mass fragments.
        quotients_relative_uncertainties = quotients*1.0 #we make this variable but only populate it if we're going to calculate uncertainties.
        #This first loop goes through the row and gets all of the mass fragment's Tms and Gms
        #and then uses those along with the relative intensity of the mass fragment itself to
        #calculate the sum of each column's Fm/(Gm*Tm)
        for row_counter in range(0,reference_height):  #array-indexed for loop
            fragment_intensity = ReferenceDataForCorrectionValues.standardized_reference_patterns[row_counter,column_counter] #This is actually the fragment's intensity in the standardized reference pattern.
            if fragment_intensity != 0: #only gets the Gm and Tm if relative intensity is not equal to zero
                electron_multiplier_gain = (28/ReferenceDataForCorrectionValues.standardized_reference_patterns[row_counter,0])**0.5 
                if (ReferenceDataForCorrectionValues.standardized_reference_patterns[row_counter,0] < 30): #transmission gain depends on the mass fragment mass
                    transmission_gain = 1
                else:
                    transmission_gain = 10**((30-ReferenceDataForCorrectionValues.standardized_reference_patterns[row_counter,0])/155)
                quotient = fragment_intensity/(transmission_gain*electron_multiplier_gain)
                quotients[row_counter] = quotient

            #calculation of uncertainties is being separated from above if statement to avoid making the code harder to read.
            if (G.calculateUncertaintiesInConcentrations == True) and (type(G.referenceFileUncertainties) != type(None)):    
                if fragment_intensity != 0: #We only get the Gm and Tm if relative intensity is not equal to zero, so we only need uncertainties for those values.
                    fragment_intensity_relative_uncertainty = ReferenceDataForCorrectionValues.relative_standard_uncertainties[row_counter,column_counter] #row_counter and column_counter already skip the mass fragments.
                    quotient_relative_uncertainty = abs(fragment_intensity_relative_uncertainty*1.0) #there is no change since this is relative uncertainty divided by some constants (the two gains are dependent on molecular weight, but not on intensity).
                    quotients_relative_uncertainties[row_counter] = quotient_relative_uncertainty
                    
        a = sum(quotients)
        if (G.calculateUncertaintiesInConcentrations == True) and (type(G.referenceFileUncertainties) != type(None)): #if uncertainties are going to be calculated.
            quotients_absolute_uncertainties = quotients*quotients_relative_uncertainties #We will need the absolute uncertainties because now there is a summation.
            #now to calculate the a_uncertainty, we we will use a_uncertainty = sqrt( uncertainty_1^2 + uncertainty_2^2 + uncertainty_3^2)
            #TODO: there is actually correlation in the errors, so the real number should be smaller than this.
            sum_in_sqrt = 0 #just initailizing this to add to it.
            for uncertainty_term in quotients_absolute_uncertainties:
                sum_in_sqrt =  sum_in_sqrt + uncertainty_term**2
            a_absolute_uncertainty = (sum_in_sqrt)**(0.5) 
            
        # This second part of the loop must be separate only because the sums of each column's 
        #Fm/(Gm*Tm) has to be found first, but then the loop gets each mass fragment's correction
        #value by dividing the sum by that relative intensity and ionization efficiency for the 
        #molecule, then creating an correction_values of correction values
        #Before starting the below loop, we can now trim the ReferenceDataForCorrectionValues (and the associated other arrays) if it has too many masses.
        #ReferenceData.ClearZeroRowsFromStandardizedReferenceIntensities()
        
        #First calculate any tuning correction intensity factors and their uncertainties.
        if len(G.referenceFileStandardTuningAndForm) > 0: #if referenceFileStandardTuningAndForm is provided, that means that the tuningCorrectorIntensity feature will be used. These factors could be calculated one at a time inside the for loop below, but due to how the function works and to make the code somewhat more efficient we are calculating this uncertainty array ahead of time and then will popualte other variables with them in the loop below.
            tuningCorrectionIntensityFactors, tuningCorrectionIntensityFactors_absolute_uncertainties = DataFunctions.returnPolyvalEstimatesAndUncertainties(x_values=ReferenceData.standardized_reference_patterns[:,0].flatten(), abcCoefficients=abcCoefficients, abcCoefficients_covMat= abcCoefficients_covmat)
            tuningCorrectionIntensityFactors_relative_uncertainties = tuningCorrectionIntensityFactors_absolute_uncertainties/tuningCorrectionIntensityFactors
        #The row_counter below is the mass fragment index.        
        for row_counter in range(0,reference_height):#array-indexed for loop, there is a second loop for each column, because the first only gets the sum of (Fm)/(Gm*Tm)
            fragment_intensity = ReferenceData.standardized_reference_patterns[row_counter,column_counter] 
            if fragment_intensity != 0: #once again, if the relative intensity is zero, then the correction value will be zero as well, this is done by the else statement
                correction = a/(ionization_efficiency*fragment_intensity)
                if len(G.referenceFileStandardTuningAndForm) > 0: #if referenceFileStandardTuningAndForm is provided, that means that the tuningCorrectorIntensity feature will be used.
                    correction = correction*tuningCorrectionIntensityFactors[row_counter]
                current_molecule_correction_factors_list.append(correction)
                correction_values_direct[row_counter,column_counter] = correction
            else: 
                correction = 0
                current_molecule_correction_factors_list.append(correction)
                correction_values_direct[row_counter,column_counter] = 0
            ##Here we calculate the correction factor uncertainties if desired.
            #To calculate the uncertainty associated with this operation: correction = a/(ionization_efficiency*fragment_intensity)
            #we use the following formula:  if Q = (a*b)/(c*d), Q_relative_uncertainty = sqrt( a_relative_uncertainty^2 + b_relative_uncertainty^2 + c_relative_uncertainty^2 + d_relative_uncertainty^2)
            #Note that there is no difference between how we propagate for multiplication and division in this formula: they can be done in the same propagation step.
            if (G.calculateUncertaintiesInConcentrations == True) and (type(G.referenceFileUncertainties) != type(None)): #if uncertainties are going to be calculated.
                if fragment_intensity != 0:
                    #For uncertainty associatesd with correction = a/(ionization_efficiency*fragment_intensity) we first gather the various separate uncertainties.
                    #We already have a_absolute_uncertainty but we need a_relative_uncertainty
                    a_relative_uncertainty = a_absolute_uncertainty/a                    
                    ionization_efficiency_relative_uncertainty = 0 #For now we make it 0. #TODO add in effects of ionization_efficiency_relative_uncertainty
                    fragment_intensity_relative_uncertainty = ReferenceData.relative_standard_uncertainties[row_counter,column_counter] #this is just reproduced from above.
                    if len(G.referenceFileStandardTuningAndForm) > 0:
                        tuningCorrectionIntensityFactor_relative_uncertainty = tuningCorrectionIntensityFactors_relative_uncertainties[row_counter] #This is already in relative uncertainties. When the equation is written out, the factor is multplied against the correction factor, such that the relative uncertainty of tuningCorrectionIntensityFactors propagates to a relative uncertainty in the correction factors.
                    else:
                        tuningCorrectionIntensityFactor_relative_uncertainty = 0
                    correction_relative_uncertainty = (a_relative_uncertainty**2 + ionization_efficiency_relative_uncertainty**2 + fragment_intensity_relative_uncertainty**2 + tuningCorrectionIntensityFactor_relative_uncertainty**2)**0.5
                    correction_values_direct_relative_uncertainty[row_counter,column_counter] = correction_relative_uncertainty
                else:
                    correction_relative_uncertainty = 0 #There is no correction factor so we use a relative uncertainty is 0.
                    correction_values_direct_relative_uncertainty[row_counter,column_counter] = correction_relative_uncertainty
    #correction_values_direct[:,0] This is the mass fragments.
    correction_values = numpy.transpose(correction_values_direct[:,1:]) #skipping one column because it has mass fragments.
    correction_values_relative_uncertainties = numpy.transpose(correction_values_direct_relative_uncertainty[:,1:])
    #G.correction_values_relative_uncertainties = correction_values_relative_uncertainties
    return correction_values, correction_values_relative_uncertainties


#Populate_reciprocal_matching_correction_values-> for loops that build the reciprocal_matching_correction_values and raw signals into multiple arrays, inside of lists
#that are ready to be solved in the next section. It does this by using a itertools combination function in order
#to create all the combinations of row and thereby creating all the possible answer arrays for the data given
#need these indexes for later
def Populate_reciprocal_matching_correction_values(mass_fragment_numbers, ReferenceData):
    ReferenceData.referenceabscissa = ReferenceData.standardized_reference_patterns[:,0]
    referenceDataArray = ReferenceData.standardized_reference_patterns[:,1:]
    correction_values = numpy.array(list(zip(*ReferenceData.correction_values)))
    if G.calculateUncertaintiesInConcentrations == True:
        if type(G.referenceFileUncertainties) != type(None): #Just mimicing the above lines. 
            correction_values_relative_uncertainties = numpy.array(list(zip(*ReferenceData.correction_values_relative_uncertainties)))
    #This function has inputs that are very general so that it could be easily understood and used in various 
    #circumstances, the function first gets the size of the data array and then uses that to index the loops
    #that find matching numbers in the abscissas and then keep those respective rows
    def ArrayRowReducer(reducedabscissa,abscissa,data): 
        reducedabscissa_length = len(reducedabscissa)
        data_width = len(data[0,:]) 
        data_height = len(data[:,0])
        matching_abscissa = []
        matching_abscissa_data = numpy.zeros([1,data_width])   
        for reducedabscissacounter in range(0,reducedabscissa_length): #array-indexed for loop
            for abscissacounter in range(0,data_height):#array-indexed for loop
                if reducedabscissa[reducedabscissacounter] == abscissa[abscissacounter]: #gets index for chosen mass fragment numbers within correction values/matching mass fragemtns arrays
                    matching_abscissa.append(abscissa[abscissacounter])
                    matching_abscissa_data = numpy.vstack([matching_abscissa_data, data[abscissacounter,list(range(0,data_width))]])
        matching_abscissa_data = numpy.delete(matching_abscissa_data,(0),axis = 0)
        return matching_abscissa_data, matching_abscissa
    #This small function just goes through every element of the correction array and inverses it; you can do
    #this more simply, but there are zeros here and we cannot have inf as our value, so the else statement simply
    #skips the zeros and inverse all others
    def ArrayElementsInverser(matching_correction_values):  #TODO: use numpy.divide for this (it is used for a similar task elsewhere in the program already, just search in program for example)
        reciprocal_matching_correction_values = copy.deepcopy(matching_correction_values)#we mean element reciprical 1/X. Not inverse Matrix
        for x in range(len(matching_correction_values[:,0])): #array-indexed for loop, these two loops go through all the values in the array
            for y in range(len(matching_correction_values[0,:])):#array-indexed for loop
                if matching_correction_values[x][y] != 0: #when a number is zero using **-1 gives a divide by zero error- so all these are skipped
                    reciprocal_matching_correction_values[x][y] = matching_correction_values[x][y]**float(-1)
        return reciprocal_matching_correction_values
    #here the main function, Populate_reciprocal_matching_correction_values, calls all of its sub-functions 
    ReferenceData.matching_correction_values, ReferenceData.matching_abscissa = ArrayRowReducer(mass_fragment_numbers,ReferenceData.referenceabscissa,correction_values)
    ReferenceData.monitored_reference_intensities, ReferenceData.matching_abscissa = ArrayRowReducer(mass_fragment_numbers,ReferenceData.referenceabscissa,referenceDataArray)
    if G.calculateUncertaintiesInConcentrations == True:
        if type(G.referenceFileUncertainties) != type(None): #Just mimicing the above lines.
            ReferenceData.reciprocal_matching_correction_values_relative_uncertainties, ReferenceData.matching_abscissa = ArrayRowReducer(mass_fragment_numbers,ReferenceData.referenceabscissa,correction_values_relative_uncertainties)
    ReferenceData.reciprocal_matching_correction_values = ArrayElementsInverser(ReferenceData.matching_correction_values)
    return ReferenceData
    
    
#This function will take the reference data and eliminate any molecules from the data that do not contain any mass fragment
#data this is done so that there will not be errors in the code later (as other molecules may also have 0 signal 
#relative to CO) It does this by looking at the matching mass fragments and deleting any columns which contain only
#zeros, and then also deletes that molecule form the molecules array and the correction values array.
def UnnecessaryMoleculesDeleter(ReferenceData):
    try:
        width = len(ReferenceData.monitored_reference_intensities[0,:])
    except:
        print("ERROR: There was an error in attempting to get the width of monitored_reference_intensities. This usually means that the Experimental data file is missing its comment line. Try adding a line with just the symbol # to the top of the experimental data file if there is not a line like this already. You may add any words you would like after the # symbol.")
    height = len(ReferenceData.monitored_reference_intensities[:,0])
    place_holder = 0
    for columncounter in range(width):#array-indexed for loop
        column = []
        for rowcounter in range(height):#array-indexed for loop
            column.append(ReferenceData.monitored_reference_intensities[rowcounter,columncounter-place_holder])
            if rowcounter == height-1: #at the end of the loop
                if sum(column) == 0:#if there are no relative intensities for the chosen mass fragments of this molecules, all its data is deleted from the arrays
                    ReferenceData.monitored_reference_intensities = numpy.delete(ReferenceData.monitored_reference_intensities,(columncounter-place_holder),axis = 1)
                    ReferenceData.reciprocal_matching_correction_values = numpy.delete(ReferenceData.reciprocal_matching_correction_values,(columncounter-place_holder),axis = 1)
                    ReferenceData.molecules = numpy.delete(ReferenceData.molecules,(columncounter-place_holder))
                    if G.calculateUncertaintiesInConcentrations == True:
                        if type(G.referenceFileUncertainties) != type(None): #Just mimicing the above lines.
                            ReferenceData.reciprocal_matching_correction_values_relative_uncertainties = numpy.delete(ReferenceData.reciprocal_matching_correction_values_relative_uncertainties,(columncounter-place_holder),axis = 1)
                    place_holder = place_holder + 1
    return ReferenceData

    
#this little function lets you choose your own times range, the inputs are the start of the range,
#the end of the range, a 'yes' or 'no' (timerangelimit), and the times and collected arrays
#the collected data will be shortened to the length of the new chosen times abscissa
def TimesChooser (ExperimentData,timeRangeStart,timeRangeFinish):
    if timeRangeStart > timeRangeFinish:
        print("ERROR: MSRESOLVE TimesChooser has received a timeRangeStart greater than timeRangeFinish. Fix your input and run again."); sys.exit
    #due to an uncommon numpy error "DeprecationWarning: in the future the special handling of scalars will be removed from delete and raise an error"  we will check if the rawsignals_absolute_uncertainties is really present or not. If it is present, we will need to operate on it as well.
    rawsignals_absolute_uncertainties_present = False #just initiailzing.
    if hasattr(ExperimentData,'rawsignals_absolute_uncertainties'):
        #check if the shape is really correct. Because it can otherwise be an array with just a single None object in it.
        if numpy.shape(ExperimentData.rawsignals_absolute_uncertainties) == numpy.shape(ExperimentData.workingData):
            rawsignals_absolute_uncertainties_present = True
        #else: #don't need to do the else explicitly since we intialized this as false.
            #rawsignals_absolute_uncertainties_present = False
    place_holder = 0 #enables indexing when parts of the array are being deleted
    deletePoint = False #just initializing.
    for timescounter in range(len(ExperimentData.times)): #array indexed for loop
        if ExperimentData.times[timescounter-place_holder] < timeRangeStart: #all rows that are before the time range are deleted from the collected data and times abscissa
            deletePoint = True
        if ExperimentData.times[timescounter-place_holder] > timeRangeFinish: #once the time is greater than the time range finish, all values after are deleted
            deletePoint = True
        if deletePoint == True:
            ExperimentData.times = numpy.delete(ExperimentData.times,timescounter-place_holder) #place holder subtracts from the for loop so that the correct index is maintained
            ExperimentData.workingData = numpy.delete(ExperimentData.workingData,timescounter-place_holder,axis = 0)
            if rawsignals_absolute_uncertainties_present == True: #copying the line from above.
                ExperimentData.rawsignals_absolute_uncertainties = numpy.delete(ExperimentData.rawsignals_absolute_uncertainties,timescounter-place_holder,axis = 0) 
            deletePoint = False #set back to the default after deleting.
            place_holder = place_holder + 1 #the place holder increased by one with every deleted row to maintain array indexing
    if len(ExperimentData.times) == 0:
        print("ERROR: TimesChooser received a timeRangeStart and timeRangeFinish that resulted in no times remaining for the Experimental data. Change the userinput and run your program again."); sys.exit()
    return None

''' ScaleDown takes an array and scales every value by the same factor so that
the largest value is below a chosen size.
Arguments:
a1DArray(required): the array to alter 
multiplier(optional): gives user option to set factor manually 
ScalesOf10(default = False): controls if the factor is also a factor of ten
Cap(default = 100): sets the chosen size '''

def ScaleDown(a1DArray, multiplier = None, Cap = 100, ScalesOf10 = False):
    #finds the factor to adjust max to 100
    if multiplier == None:
        #finds the largest entry
        maxNumber = float(max(a1DArray))
        # Confirm that the array needs scaling
        if maxNumber < Cap:
            multiplier = 1
            return a1DArray, multiplier
        #calculate multiplier
        multiplier = Cap/maxNumber
        # if neccessary, scale multiplier 
        if ScalesOf10:
             multiplier = math.pow(10, (math.floor(math.log(multiplier, 10))))
    # If given a multiplier, only need to make sure it's a float
    else:
        multiplier = float(multiplier)
    #applies factor to all entries
    for i in range(0,len(a1DArray)):
        a1DArray[i] = (a1DArray[i]*multiplier)
    #Returns the array  
    return a1DArray, multiplier

''' ScaleUp takes an array and scales every value by the same factor so that
the smallest value is above a chosen size.
Arguments:
a1DArray(required): the array to alter 
multiplier(optional): gives user option to set factor manually 
ScalesOf10(default = False): controls if the factor is also a factor of ten
Base(default = 1): sets the chosen size '''

def ScaleUp(a1DArray, multiplier = None, Base = 1, ScalesOf10 = False):
    #makes the collected array into a 1D numpy array
    a1DArray = numpy.array(a1DArray)
    #finds the factor to adjust min to 1
    if multiplier == None:
        #finds the smallest entry
        minNumber = float(numpy.min(a1DArray[numpy.nonzero(a1DArray)]>0))  #This first gets the nonzero values, then takes the ones greater than 0, then finds the minimum.
        # Confirm that the array needs scaling
        if minNumber > Base:
            multiplier = 1
            return a1DArray, multiplier
        # calculate multiplier
        multiplier = Base/minNumber
        # if neccessary, scale multiplier 
        if ScalesOf10:
            multiplier = math.pow(10, (math.ceil(math.log(multiplier, 10))))
    # If given a multiplier, only need to make sure it's a float
    else:        
        multiplier = float(multiplier)
    #applies factor to all entries
    for i in range(0,len(a1DArray)):
        a1DArray[i] = (a1DArray[i]*multiplier)
    #Returns the array       
    return a1DArray, multiplier

def ScaleRawData(data, scaleRawDataOption, scaleRawDataFactor):
    scaledData = copy.deepcopy(data)                                
    if scaleRawDataOption == 'auto':
        #automatically scales the smallest point to 1
        scaledData, multiplier = ScaleRawDataAuto(data)
    elif scaleRawDataOption == 'manual':
        #because a multiplier is given, it doesn't matter if ScaleUp or ScaleDown is used
        scaledData, multiplier = ScaleUp(data, multiplier = scaleRawDataFactor)
    return scaledData, multiplier

''' ScaleRawDataAuto operates in a similar way to ScaleUp or Down, except that 
    it will always scale the smallest value to 1, regardless of the scaling direction'''
def ScaleRawDataAuto(data): 
    a1Dholder=[]
    a1Dholder=data[0]
    minNum=numpy.min(a1Dholder[numpy.nonzero(a1Dholder)])

    for index in range(0, len(data)):
        a1Dholder =data[index]
        
        minNumHolder=numpy.min(a1Dholder[numpy.nonzero(a1Dholder)]>0) #This first gets the nonzero values, then takes the ones greater than 0, then finds the minimum.
        if ((minNumHolder>0) and (minNumHolder < minNum)):
            minNum = minNumHolder
            
    multiplier= 1/minNum

    for dataPointIndex in range(0, len(data)):
        for referenceIndex in range(0, len(data[0])):
            data[dataPointIndex, referenceIndex]=data[dataPointIndex, referenceIndex]*multiplier
    return data, multiplier

'''
Standardize is a simple algorithim that reads in one numpy array
(oneArray) and scales the whole array so that the
greatest number is 100 within the array
'''

def StandardizeTo100(a1DArray,n):
    a1DArrayStandardized = a1DArray.copy() #making a copy so original array is not changed.
    maxNumber=float(numpy.max(a1DArrayStandardized))
    multiplier=100.0/maxNumber
    for index2 in range(0,len(a1DArrayStandardized)):
        a1DArrayStandardized[index2]=(a1DArrayStandardized[index2]*multiplier)
    return a1DArrayStandardized



'''
StandardizeReferencePattern uses StandardizeTo100 to modify reference values so that each column is scaled to
100. NOTE THAT THIS FUNCTION ASSUMES THAT THE FIRST COLUMN in reference contains the 
mass framgnet numbers. 
Parameters: 
standardizedReference -  a name chosen for the numpy array that contains reference values
num_of_molecues-  an integer describing the number of  molecues that contributed to the reference file. 
If the num_of_molecues is fed in as 0 0 or not provided it will be calculated automatically.
'''
def StandardizeReferencePattern(referenceUnstandardized,num_of_molecules=0):
    if num_of_molecules==0:
        num_of_molecules = len(referenceUnstandardized[0]) -1 #subtract 1 for the abscissa.
    # preallocate new array for standardized values
    standardizedReference = copy.deepcopy(referenceUnstandardized)
    # standardize
    for moleculeIndex in range(1,num_of_molecules+1): #Note that we start at an index of 1 in order to skip the mass fragments (so they don't get 'standardized')
        standardizedReference[0:,moleculeIndex]=StandardizeTo100(referenceUnstandardized[0:,moleculeIndex],1) #TODO: Low priority. Change this to use amax so the loop isn't executed.
    return standardizedReference

'''This function is just a helper function to create a ReferenceData object from a file.'''
def createReferenceDataObject(referenceDataFileName = "ReferenceData.csv", referenceDataForm = "xyyy", AllMID_ObjectsDict={}):
    [provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form]=readReferenceFile(referenceDataFileName,referenceDataForm)
    ReferenceData = MSReference(provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName=referenceFileName, form=form, AllMID_ObjectsDict=AllMID_ObjectsDict)
    return ReferenceData
    
'''The following two functions are currently not used in the program,
but have been saved in case they are needed in the future'''
def CanBeFloat(value):
  return (type(value) == int or type(value) == float)

def ImportWorkingData(preProcessedDataOutputName):

    dataFrame = pandas.read_csv('%s' %preProcessedDataOutputName, header = None)
    # While we use ExportXYYYData() to write this data file out
    # we will need to remove the last column because
    # ExportXYYYData() ends every row with a comma
    # pandas interprets the empty value after the trailing comma
    # as a 'nan'. Thus the last column composed of all nans.
    # Make sure that's the case and then get rid of them
    if all(numpy.isnan(dataFrame.iloc[0:,-1])):
        dataFrame = dataFrame.iloc[:,0:-1]

    ''' generate mass fragment list'''
    #select only the 1st row down, all columns except for the first
    dfmass = dataFrame.iloc[0][1:]
    #convert to matrix
    masses = dfmass.values
    #sort through the matrix and remove labels
    #masses = numpy.delete(masses, -1)
    for i in range(0,len(masses)):
        masses[i] = str(masses[i]).replace('m','')
    #convert the matrix to floats if they aren't already 
    mass_fragment_numbers = masses.astype(float)
    
    '''generate time list'''
    #select column of times
    dftimes = dataFrame.iloc[1:][0]
    #convert to matrix
    times = dftimes.values
    #save with type float
    fulltimes = times.astype(float)
    
    '''collect preprocessed data'''
    #select matrix of signals
    dfpreprocessed = dataFrame.iloc[1:,1:]
    #convert to matrix
    preprocessed = dfpreprocessed.values
    #save  with type float
    preprocessedData = preprocessed.astype(float)

    
    '''create data set to work on'''
    workingData = preprocessedData
    
    return [workingData, mass_fragment_numbers, fulltimes]

def ImportAnalyzedData(concentrationsOutputName):

    dataFrame = pandas.read_csv('%s' %concentrationsOutputName, header = None)

    '''collect preprocessed data'''
    #select matrix of signals
    dfanalyzed = dataFrame.iloc[1:,0:]
    #convert to matrix
    analyzed = dfanalyzed.values
    analyzed = numpy.delete(analyzed, -1, 1)
    #save  with type float
    analyzedData = analyzed.astype(float)
    
    return analyzedData



'''
Performs some manipulations related to the reference pattern
'''
def ReferenceInputPreProcessing(ReferenceData, verbose=True):
    #We will skip the regular createReferencePatternWithTuningCorrection if using the gas mixture feature.
    try:
        if len(G.UserChoices['measuredReferenceYorN']['tuningCorrectorGasMixtureMoleculeNames']) == 0:
            ReferenceData = createReferencePatternWithTuningCorrection(ReferenceData, verbose=verbose)
    except: #This try and except is so that iterative can work correctly.
        ReferenceData = createReferencePatternWithTuningCorrection(ReferenceData, verbose=verbose)
    #TODO: the minimal reference value can cause inaccuracies if interpolating between multiple reference patterns if one pattern has a value rounded to 0 and the other does not
    #TODO: option 1: this issue can be fixed by moving this to after interpolation
    #TODO: option 2: Or we can below assign to preprocessed_reference_pattern rather than standardized_reference_patterns and then use that in data analysis (Note that interpolate would continue to use standardized_reference_patterns as well as preprocess the output)
    if G.minimalReferenceValue == 'yes':
        if G.implicitSLScorrection == True: #This feature requires us to have unfiltered reference patterns to do an implicit/recursive correction at the end.
            G.currentReferenceDataUnfiltered = copy.deepcopy(ReferenceData) #Make a copy before any filtering occurs. There is an implied return in global variable.
        ReferenceData.standardized_reference_patterns = ReferenceThresholdFilter(ReferenceData.standardized_reference_patterns,G.referenceValueThreshold)
        ReferenceData.ExportCollector('ReferenceThresholdFilter')
    
    #As the program is currently written, this function is called to act upon already threshold filtered standardized reference patterns which could cause innaccuracy.  
    #One could move this function prior to threshold filtering however then correction values would not be correctly calculated for interpolated reference patterns
    #We are not sure there are any other reasons we can't move this function call. However, there may be some care needed when using uncertainties.
    ReferenceData.correction_values, ReferenceData.correction_values_relative_uncertainties = CorrectionValuesObtain(ReferenceData)
    if G.implicitSLScorrection == True: #This feature requires us to have unfiltered reference patterns to do an implicit/recursive correction at the end. There is an implied return in global variable.
            G.currentReferenceDataUnfiltered.correction_values, G.currentReferenceDataUnfiltered.correction_values_relative_uncertainties = CorrectionValuesObtain(G.currentReferenceDataUnfiltered)
    
    #Only print if not called from interpolating reference objects
    if verbose:
        print('CorrectionValuesObtain')

    return ReferenceData

'''
GenerateReferenceDataAndFormsList takes in the list of referenceFileNamesList and the
list of forms.  A list is generated containing MSReference objects created based
on the referenceFileName and the corresponding form
It allows MSRESOLVE to be backwards compatible with previous user input files
'''
def GenerateReferenceDataList(referenceFileNamesList,referenceFormsList,AllMID_ObjectsDict={}):
    #referenceFormsList can take values of 'xyyy' or 'xyxy' and must be a string
    ##If referenceFileNamesList is a string or if form is a string then make them lists
    if isinstance(referenceFileNamesList,str):
        referenceFileNamesList = [referenceFileNamesList]
    if isinstance(referenceFormsList,str):
        referenceFormsList = [referenceFormsList]
    #If referenceFileNamesList and forms are lists of 1 then create a list of the single MSReference object
    #This allows MSRESOLVE to be backwards compatible with previous user input files while still incorporating the reference pattern time chooser feature
    if len(referenceFormsList) == 1 and len(referenceFileNamesList) == 1:
        [provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form]=readReferenceFile(referenceFileNamesList[0],referenceFormsList[0])
        ReferenceDataList = [MSReference(provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName=referenceFileName, form=form, AllMID_ObjectsDict=AllMID_ObjectsDict)]
        #save each global variable into the class objects
        ReferenceDataList[0].ExportAtEachStep = G.ExportAtEachStep
        ReferenceDataList[0].iterationSuffix = G.iterationSuffix
        if G.calculateUncertaintiesInConcentrations == True:
            if type(G.referenceFileUncertainties) != type(None):
                if type(G.referenceFileUncertainties) == type(float(5)) or  type(G.referenceFileUncertainties) == type(int(5)) :
                    #TODO: Low priority. The below results in "nan" values. It could be better to change it to make zeros using a numpy "where" statement.
                    G.referenceFileUncertainties = float(G.referenceFileUncertainties) #Maks sure we have a float.
                    #Get what we need.
                    provided_reference_patterns = ReferenceDataList[0].provided_reference_patterns
                    provided_reference_patterns_without_masses = ReferenceDataList[0].provided_reference_patterns[:,1:] #[:,0] is mass fragments, so we slice to remove those. 
                    #Make our variables ready.
                    absolute_standard_uncertainties = provided_reference_patterns*1.0 #First we make a copy.
                    absolute_standard_uncertainties_without_masses = (provided_reference_patterns_without_masses/provided_reference_patterns_without_masses)*G.referenceFileUncertainties #We do the division to get an array of ones. Then we multiply to get an array of the same size.
                    #now we populate our copy's non-mass area.
                    absolute_standard_uncertainties[:,1:] = absolute_standard_uncertainties_without_masses                                        
                    ReferenceDataList[0].absolute_standard_uncertainties = absolute_standard_uncertainties
                    #We can't convert to relative uncertainties yet because the file may not be standardized yet.
                if type(G.referenceFileUncertainties) == type('string'):
                    if '.csv' in referenceFileNamesList[0]:
                        provided_reference_patterns_absolute_uncertainties, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form = readReferenceFile(referenceFileNamesList[0][:-4]+"_absolute_uncertainties.csv",referenceFormsList[0])
                        ReferenceDataList[0].absolute_standard_uncertainties = provided_reference_patterns_absolute_uncertainties #Just initializing the variable before filling it properly.
                    if '.tsv' in referenceFileNamesList[0]:
                        provided_reference_patterns_absolute_uncertainties, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form = readReferenceFile(referenceFileNamesList[0][:-4]+"_absolute_uncertainties.tsv",referenceFormsList[0])
                        ReferenceDataList[0].absolute_standard_uncertainties = provided_reference_patterns_absolute_uncertainties #Just initializing the variable before filling it properly.
                    maximum_absolute_intensities = numpy.amax(ReferenceDataList[0].provided_reference_patterns[:,1:], axis = 0) #Find the maximum intensity for each molecule.
                    ReferenceDataList[0].absolute_standard_uncertainties[:,1:] = 100*ReferenceDataList[0].absolute_standard_uncertainties[:,1:]/maximum_absolute_intensities
                    #TODO: low priority, remove nan values and/or populate them with zero using numpy divide.
        #The below process **normally** happened later, but because of TuningCorrectorGasMixture, this was needed earlier so has been copied to occur earlier also. It should not cause any harm.        
        if G.calculateUncertaintiesInConcentrations == True: 
            if type(G.referenceFileUncertainties) != type(None):
                ReferenceDataList[0].update_relative_standard_uncertainties()
        return ReferenceDataList
    #Otherwise we have multiple reference files and forms
    #If just one form is used, make a list of forms that is the same length as referenceFileNamesList
    if len(referenceFormsList) == 1:
        #Generate a copy of referenceFileNamesList to be overwritten with forms
        listOfForms = copy.copy(referenceFileNamesList)
        #replace each value with the given form
        for i in range(len(referenceFileNamesList)):
            listOfForms[i] = referenceFormsList[0]
    #If list of forms is the same length of referenceFileNamesList then each form should correspond to the referenceFile of the same index
    elif len(referenceFormsList) == len(referenceFileNamesList):
        #So just set listOfForms equal to forms
        listOfForms = referenceFormsList
    #Initialize ReferenceDataList so it can be appended to
    ReferenceDataList = []
    #For loop to generate each MSReferenceObject and append it to a list
    for i in range(len(referenceFileNamesList)):
        [provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form]=readReferenceFile(referenceFileNamesList[i],listOfForms[i])
        ReferenceDataList.append(MSReference(provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName=referenceFileName, form=form, AllMID_ObjectsDict=AllMID_ObjectsDict))
        #save each global variable into the class objects 
        ReferenceDataList[i].ExportAtEachStep = G.ExportAtEachStep
        ReferenceDataList[i].iterationSuffix = G.iterationSuffix
        if G.calculateUncertaintiesInConcentrations == True:
            if type(G.referenceFileUncertainties) != type(None):
                if type(G.referenceFileUncertainties) == type(float(5)) or  type(G.referenceFileUncertainties) == type(int(5)) :
                    #TODO: The below results in "nan" values. It would be better to change it to make zeros using a numpy "where" statement.
                    G.referenceFileUncertainties = float(G.referenceFileUncertainties) #Maks sure we have a float.
                    #Get what we need.
                    provided_reference_patterns = ReferenceDataList[i].provided_reference_patterns
                    provided_reference_patterns_without_masses = ReferenceDataList[i].provided_reference_patterns[:,1:] #[:,0] is mass fragments, so we slice to remove those. 
                    #Make our variables ready.
                    absolute_standard_uncertainties = provided_reference_patterns*1.0 #First we make a copy.
                    absolute_standard_uncertainties_without_masses = (provided_reference_patterns_without_masses/provided_reference_patterns_without_masses)*G.referenceFileUncertainties #We do the division to get an array of ones. Then we multiply to get an array of the same size.
                    #now we populate our copy's non-mass area.
                    absolute_standard_uncertainties[:,1:] = absolute_standard_uncertainties_without_masses                                        
                    ReferenceDataList[i].absolute_standard_uncertainties = absolute_standard_uncertainties
                    #We can't convert to relative uncertainties yet because the file may not be standardized yet.
                if type(G.referenceFileUncertainties) == type('string'):
                    if '.csv'in referenceFileNamesList[0]:
                        provided_reference_patterns_absolute_uncertainties, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form = readReferenceFile(referenceFileNamesList[0][:-4]+"_absolute_uncertainties.csv",referenceFormsList[i])
                    if '.tsv'in referenceFileNamesList[0]:
                        provided_reference_patterns_absolute_uncertainties, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form = readReferenceFile(referenceFileNamesList[0][:-4]+"_absolute_uncertainties.tsv",referenceFormsList[i])
                    ReferenceDataList[i].absolute_standard_uncertainties = provided_reference_patterns_absolute_uncertainties #Just initializing the variable before filling it properly.
                    maximum_absolute_intensities = numpy.amax(ReferenceDataList[i].provided_reference_patterns[:,1:], axis = 0) #Find the maximum intensity for each molecule.
                    ReferenceDataList[i].absolute_standard_uncertainties[:,1:] = ReferenceDataList[i].absolute_standard_uncertainties[:,1:]/maximum_absolute_intensities
                    #TODO: low priority, remove nan values and/or populate them with zero using numpy divide.
        #The below process **normally** happened later, but because of TuningCorrectorGasMixture, this was needed earlier so has been copied to occur earlier also. It should not cause any harm.        
        if G.calculateUncertaintiesInConcentrations == True: 
            if type(G.referenceFileUncertainties) != type(None):
                ReferenceDataList[i].update_relative_standard_uncertainties()
    return ReferenceDataList

'''
InterpolateReferencePatterns is a function used in the reference pattern time chooser feature when a gap occurs between two time ranges
It inputs two reference files, the current time, the start time of the gap, and the end time of the gap
It interpolates row by row and returns the new interpolated data
As the program is currently written, this function is called to act upon already threshold filtered standardized reference patterns which could cause innaccuracy.  This prevents discontinuities in the interpolated intensities, in prinicpal it would sometimes be more accurate to do threshold filtering after interpolation
'''
def InterpolateReferencePatterns(firstReferenceObject,secondReferenceObject,time,gapStart,gapEnd):
    #Since we probably need values from firstReferenceObject for another interpolation, create a deepcopy of firstReferenceObject
    newReferenceObject = copy.deepcopy(firstReferenceObject)
    
    #loop through the provided reference intensities and linearly interpolate row by row
    for i in range(len(firstReferenceObject.standardized_reference_patterns)):
        #Overwrite provided_reference_patterns by interpolating the standardized_reference_patterns
        #[i,:] for every column in the ith row
        newReferenceObject.provided_reference_patterns[i,:] = DataFunctions.analyticalLinearInterpolator(firstReferenceObject.standardized_reference_patterns[i,:],secondReferenceObject.standardized_reference_patterns[i,:],time,gapStart,gapEnd)
        newReferenceObject.standardized_reference_patterns =StandardizeReferencePattern(newReferenceObject.provided_reference_patterns)
    return newReferenceObject

'''
CheckCurrentTimeRange is a function used for the Reference Pattern Time Chooser feature.  It looks at the current time in data analysis and determines
which reference pattern needs to be used based on user input.  If the time is in the current reference pattern's time range, the function does nothing.
If the time is in between two time ranges, the function calls InterpolateReferencePatterns where the two patterns are linearly interpolated.
If the time is at the beginning of the next time range, it will change the currentReferenceData to the nextReferenceData
#TODO: This is not a good algorithm as written. Instead of assuming that points can only be between current and next, the function should consider the entire referencePatternTimeRanges to see where the currentTime falls.  Right now it  only considers likely possibilities rather than all possibilities.
'''
def SelectReferencePattern(currentReferencePatternIndex, referencePatternTimeRanges, currentTime, firstReferenceObject, secondReferenceObject, ReferenceDataList):
    #NOTE: If time ranges overlap, SelectReferencePattern will select the reference pattern with the same list index as the first time range involved in the overlap
    #Print a warning if user has not filled time ranges from data analysis start and stop time
    if (currentTime > referencePatternTimeRanges[-1][1]) or (currentTime < referencePatternTimeRanges[0][0]):
        print("WARNING: User has chosen to use Reference Pattern Time Chooser.  \nUser needs to input reference pattern time ranges that fill the entirety of the data analysis time range. \nUser has not and the program is about to crash.")
        #If in the current time range, continue on with the for loop
    if currentTime >= referencePatternTimeRanges[currentReferencePatternIndex][0] and currentTime <= referencePatternTimeRanges[currentReferencePatternIndex][1]:
        pass
        #Return the original reference pattern and the index
        return (firstReferenceObject, currentReferencePatternIndex)
    #Otherwise we are not in the current time range so look for a gap or see if we are in the next time range
    else:
        #If we are out of the first time range but not yet in the second time range, we are in a gap
        if currentTime > referencePatternTimeRanges[currentReferencePatternIndex][1] and currentTime < referencePatternTimeRanges[currentReferencePatternIndex+1][0]:
            #Reads in current reference pattern, next reference pattern, current time, gap start time, and gap end time
            currentReferenceData = InterpolateReferencePatterns(firstReferenceObject,secondReferenceObject,currentTime,referencePatternTimeRanges[currentReferencePatternIndex][1],referencePatternTimeRanges[currentReferencePatternIndex+1][0])
            #Prepare the current reference data
            extractReferencePatternFromDataOptionHere = "no" #should be extracted already.
            currentReferenceData = PrepareReferenceObjectsAndCorrectionValues(currentReferenceData, ExperimentData.mass_fragment_numbers, ExperimentData, extractReferencePatternFromDataOptionHere, G.rpcMoleculesToChange, G.rpcMoleculesToChangeMF, G.rpcTimeRanges, verbose=False)
            if G.iterativeAnalysis: #If using iterative analysis, interpolate the subtracted signals' matching correction factors between the two reference objects
                currentReferenceData.SSreciprocal_matching_correction_values = DataFunctions.analyticalLinearInterpolator(firstReferenceObject.SSreciprocal_matching_correction_values,secondReferenceObject.SSreciprocal_matching_correction_values,currentTime,referencePatternTimeRanges[currentReferencePatternIndex][1],referencePatternTimeRanges[currentReferencePatternIndex+1][0])
        #If we are out of the first time range, not in a gap, and not in the last time range, then we are in the next time range
        elif currentTime >= referencePatternTimeRanges[currentReferencePatternIndex+1][0]:
            #Increase the index
            currentReferencePatternIndex = currentReferencePatternIndex + 1
            #Change the reference data accordingly
            currentReferenceData = ReferenceDataList[currentReferencePatternIndex]
        return (currentReferenceData, currentReferencePatternIndex)



# DataInput just asks the user for the data file names and indexes them into arrays 
#The input prompts are given once the instructions are printed out (lines12 and 15)
#then the data is taken by sending that input to the genfromtxt function, which 
#puts that data into arrays. It also gets the the mass fragments that are needed in 
# the collected data array and makes a one dimensional array out of that
#The first step is getting all of the data into arrays inside this script, for other 
#unctions to use, this includes being able to analyze excel sheets of the form xyyy 
#and xyxy. (Mass fragments and then data) This is all in the function DataInput ()
def DataInputPreProcessing(ExperimentData): 
    #records time of all reference data preprocessing
    #ExperimentData.ExportCollector("PreProcessing ReferenceData")
    
    #We do only scale the raw data if we are not doing iterative, or are in the first iteration.
    if G.iterativeAnalysis == False or G.iterationNumber == 1:
        #Scaling Raw Data
        #This if statement only applies ScaleRawData if scaling is set to automatic
        # or if there is a manual factor that is not 1.
        #We skip this function if manual factor was 1 because that would not change the data.
        if (G.scaleRawDataOption == 'auto' or (G.scaleRawDataOption == 'manual' and G.scaleRawDataFactor != 1)):
            ExperimentData.workingData, multiplier  = ScaleRawData(ExperimentData.workingData, G.scaleRawDataOption, G.scaleRawDataFactor)
            G.scaleRawDataFactor = multiplier                                 
            ExperimentData.ExportCollector("ScaleRawData")

    #displays graph of raw data
    if G.grapher == 'yes':
        print("Raw Signal Graph")
        Draw(ExperimentData.times, ExperimentData.workingData, ExperimentData.mass_fragment_numbers, 'no', 'Amp', graphFileName = 'rawData', fileSuffix = G.iterationSuffix, label="Raw Signal Graph", stopAtGraphs=G.stopAtGraphs, figureNumber=G.lastFigureNumber+1)
        G.lastFigureNumber = G.lastFigureNumber+1
        
    if len(G.backgroundMassFragment) != 0:
        SlopeEliminator (ExperimentData,G.backgroundMassFragment,G.backgroundSlopes,G.backgroundIntercepts)
        print('Linear Baseline Correction, Manual, Complete')
        ExperimentData.ExportCollector("SlopeEliminator")

    if G.linearBaselineCorrectionSemiAutomatic   == 'yes': #the data edit sheet is used here, to determine to run this function or not
        LinearBaselineCorrectorSemiAutomatic(ExperimentData, G.baselineType, G.massesToBackgroundCorrect, G.earlyBaselineTimes, G.lateBaselineTimes)
        print('Linear Baseline Correction, Semiautomatic, Complete')
        ExperimentData.ExportCollector("LinearBaselineCorrectorSemiAutomatic")
        
    if G.lowerBoundThresholdChooser == 'yes':#calls above function if option is desired in data edit file
        LowerBoundThresholdFilter (ExperimentData, G.massesToLowerBoundThresholdFilter, G.lowerBoundThresholdPercentage, G.lowerBoundThresholdAbsolute)
        print('LowerBound Threshold Filter Complete')
        ExperimentData.ExportCollector("LowerBoundThresholdFilter")
        
    #displays mid-preprocessing graph    
    if G.grapher == 'yes':
        print("Pre-marginalChangeRestrictor Graph")
        Draw(ExperimentData.times, ExperimentData.workingData, ExperimentData.mass_fragment_numbers, 'no', 'Amp', graphFileName ='midProcessingGraph', fileSuffix = G.iterationSuffix, label="Pre-marginalChangeRestrictor Graph", stopAtGraphs=G.stopAtGraphs, figureNumber=G.lastFigureNumber+1)
        G.lastFigureNumber = G.lastFigureNumber+1

    #Need to add uncertainties somewhere before interpolation and before datasmoother for the case that it is going to be provided as an input file.
    if G.calculateUncertaintiesInConcentrations:
        if type(G.collectedFileUncertainties) == type("String"):
            if G.collectedFileUncertainties.lower() == 'auto':
                if G.UserChoices['dataSmootherYorN']['on'].lower() == "yes":
                    G.UserChoices['uncertainties']['collectedFileUncertainties_radiusType'] = G.dataSmootherChoice
                    UncertaintiesWindowsTimeRadius = G.dataSmootherTimeRadius
                    UncertaintiesWindowsPointRadius = G.dataSmootherPointRadius
                    UncertaintiesType = "standardDeviation"
                else:
                    G.UserChoices['uncertainties']['collectedFileUncertainties_radiusType'] = 'pointrange' #This is the hardcoded default choice for 'auto'
                    UncertaintiesWindowsPointRadius = 5 # This is the hardcoded default choice for 'auto'
                    UncertaintiesWindowsTimeRadius = 5 #not relevant for 'auto'
                    UncertaintiesType = "aggregateError"                       
                UncertaintiesFromData, AverageResidualsFromData = DataFunctions.UncertaintiesFromLocalWindows(ExperimentData.workingData, ExperimentData.times, ExperimentData.mass_fragment_numbers, UncertaintiesWindowsChoice=G.UserChoices['uncertainties']['collectedFileUncertainties_radiusType'], UncertaintiesWindowsPointRadius=UncertaintiesWindowsPointRadius,  UncertaintiesWindowsTimeRadius=UncertaintiesWindowsTimeRadius,UncertaintiesType=UncertaintiesType)
                ExperimentData.rawsignals_absolute_uncertainties = UncertaintiesFromData
                ExperimentData.rawsignals_average_residuals = AverageResidualsFromData
            if G.collectedFileUncertainties.lower() == 'none' :
                ExperimentData.rawsignals_absolute_uncertainties = None
                ExperimentData.rawsignals_average_residuals = None
                G.collectedFileUncertainties = None
        if type(G.collectedFileUncertainties) == type(1): #If it's an integer then we will use that as the pointradius/timeradius.
                if G.UserChoices['dataSmootherYorN']['on'].lower() == "yes": #If smoothing, we don't use the residuals for errors.
                    UncertaintiesType = "standardDeviation"
                else:
                    UncertaintiesType = "aggregateError"            
                UncertaintiesFromData, AverageResidualsFromData = DataFunctions.UncertaintiesFromLocalWindows(ExperimentData.workingData, ExperimentData.times, ExperimentData.mass_fragment_numbers, UncertaintiesWindowsChoice=G.UserChoices['uncertainties']['collectedFileUncertainties_radiusType'],UncertaintiesWindowsPointRadius=G.collectedFileUncertainties, UncertaintiesWindowsTimeRadius=G.collectedFileUncertainties, UncertaintiesType=UncertaintiesType)
                ExperimentData.rawsignals_absolute_uncertainties = UncertaintiesFromData
                ExperimentData.rawsignals_average_residuals = AverageResidualsFromData
        elif type(G.collectedFileUncertainties) == type([1]) or type(G.collectedFileUncertainties) == type(numpy.array([1])): #if it's a list or a numpy array, it should be the same length as the number of mass fragments and will be populated for all times.
                ExperimentData.rawsignals_absolute_uncertainties = numpy.ones(numpy.shape(ExperimentData.workingData)) #just initializing as an array of ones.
                uncertaintiesPerTime = numpy.array(G.collectedFileUncertainties)
                ExperimentData.rawsignals_absolute_uncertainties = ExperimentData.rawsignals_absolute_uncertainties*uncertaintiesPerTime        
        elif type(G.collectedFileUncertainties) == type(None):
            ExperimentData.rawsignals_absolute_uncertainties = None
            ExperimentData.rawsignals_average_residuals = None

    if G.interpolateYorN == 'yes':
        [ExperimentData.workingData, ExperimentData.times] = DataFunctions.marginalChangeRestrictor(ExperimentData.workingData, ExperimentData.times, G.marginalChangeRestriction, G.ignorableDeltaYThreshold)
        if G.dataRangeSpecifierYorN == 'yes':#if the datafromcsv file does not exist(in the case that it is not chosen) then the function call cannot include it
            #Gathering data from the datarange csv
            ExperimentData.datafromcsv = genfromtxt( '%s' %G.csvFileName, delimiter=',',skip_header=1) 
            #In order for the datafromcsv file to be used for the data analysis, it must be the same size as the interpolated data. The interpolate accompanying arrays
	    #function matches the abscissa of the interpolated data with that of the accomanying data from csv file by interpolating the rest of the values between rows
	    #for each abscissa value that needs to be inserted
            ExperimentData.datafromcsv=DataFunctions.interpolateAccompanyingArrays(ExperimentData.times, ExperimentData.datafromcsv)
     
        print('Marginal Change Restrictor Finished')
        ExperimentData.ExportCollector("Marginal Change Restrictor")
        
    if G.timeRangeLimit == 'yes':
        print('Timechooser')
        TimesChooser(ExperimentData, G.timeRangeStart, G.timeRangeFinish)
        ExperimentData.ExportCollector("TimeChooser")
        
    if G.dataSmootherYorN == 'yes':
        print("DataSmoother")
        ExperimentData.workingData = DataFunctions.DataSmoother(ExperimentData.workingData, ExperimentData.times, ExperimentData.mass_fragment_numbers, G.dataSmootherChoice, G.dataSmootherTimeRadius, G.dataSmootherPointRadius, G.dataSmootherHeadersToConfineTo, G.polynomialOrder)
        ExperimentData.ExportCollector("DataSmoother")

    return ExperimentData

'''
PrepareReferenceObjectsAndCorrectionValues takes in ReferenceData to be prepared for data analysis 
'''
#TODO: remove the argument "extractReferencePatternFromDataOption" from this function. That feature has been moved outside of this function.
def PrepareReferenceObjectsAndCorrectionValues(ReferenceData, massesOfInterest=[], ExperimentData=None, extractReferencePatternFromDataOption='no', rpcMoleculesToChange=[], rpcMoleculesToChangeMF=[[]], rpcTimeRanges=[[]], verbose=True):
    # Reference Pattern Changer
    if extractReferencePatternFromDataOption == 'yes':
        ReferenceData.provided_reference_patterns, ReferenceData.provided_reference_patterns_absolute_uncertainties = ExtractReferencePatternFromData(ExperimentData, ReferenceData, rpcMoleculesToChange, rpcMoleculesToChangeMF, rpcTimeRanges)
        ReferenceData.ExportCollector('ExtractReferencePatternFromData',use_provided_reference_patterns = True)
        ReferenceData.ExportCollector('ExtractReferencePatternFromData_absolute_uncertainties', export_uncertainties= True)
        #Only print if not called from interpolating reference objects
        if verbose:
            print('ReferencePatternChanger complete')
    # Some initial preprocessing on the reference data
    ReferenceData = ReferenceInputPreProcessing(ReferenceData, verbose) #Note: if implicitSLScorrection is being used, G.currentReferenceDataUnfiltered gets created inside ReferenceInputPreProcessing
    # Set the ReferenceData.monitored_reference_intensities and
    # ReferenceData.reciprocal_matching_correction_values fields
    # based on the massesOfInterest, which is typically the ExperimentData.mass_fragment_numbers
    if len(massesOfInterest) > 0:
        ReferenceData = Populate_reciprocal_matching_correction_values(massesOfInterest,ReferenceData)
        if G.implicitSLScorrection == True: #if implicitSLS correction is being used, we need to do it for the unfiltered reference pattern also. 
            G.currentReferenceDataUnfiltered = Populate_reciprocal_matching_correction_values(massesOfInterest,G.currentReferenceDataUnfiltered)
        #Exports the matching correction value 
        ReferenceData.ExportCollector('ReciprocalMatchingCorrectionValues', export_matching_correction_value = True)              #exports matching correction value 
        #Only print if not called from interpolating reference objects
        if verbose:
            print("Matching Correction Values complete")
    # Remove reference species that have no mass fragment data
    # from the ReferenceData fields monitored_reference_intensities, reciprocal_matching_correction_values and molecules
    ## TODO: Consider changing this function to take the array directly i.e.
    ## (monitored_reference_intensities) so that it can potentially be applied to other arrays
    ## like ReferenceData.standardized_reference_patterns
    ReferenceData = UnnecessaryMoleculesDeleter(ReferenceData)
    if G.implicitSLScorrection == True: #if implicitSLS correction is being used, we need to do it for the unfiltered reference pattern also.
        G.currentReferenceDataUnfiltered = UnnecessaryMoleculesDeleter(G.currentReferenceDataUnfiltered)
    ReferenceData.ExportCollector('UnnecessaryMoleculesDeleter')
    # Export the reference data files that have been stored by ReferenceData.ExportCollector
    ReferenceData.ExportFragmentationPatterns(verbose)
    return ReferenceData

#Will take current reference data monitored intensities and will set any molecules to zero that do not meet their reference threshold.
def signalThresholdFilter(ReferenceDataObject, rawsignalsarrayline, minimumSignalRequired, minimumStandardizedReferenceHeightToBeSignificant):
    ReferenceDataObject = copy.deepcopy(ReferenceDataObject)    
    #first we flatten the rawSignalsArrayLine because we want it to be 1D.
    flattenedRawSignalsArrayLine = rawsignalsarrayline.flatten()
    #the below line creates a Boolean array which is true when the signal is below the minimum signal required.
    signalsNegligible = flattenedRawSignalsArrayLine < minimumSignalRequired
    #the below line creates a Boolean array which is true when the standardized reference intensity is above the threshold to be declared a significant peak
    signiciant_Peaks_locations_monitored_reference_intensities = ReferenceDataObject.monitored_reference_intensities > minimumStandardizedReferenceHeightToBeSignificant
    indexesOfMoleculesSetToZero = []
    #now we are going to loop across each molecules, if we find the case where the signals are negligible but the peak is significant, then we store that as a molecule that is not present.
    for moleculeIndex in range(len(signiciant_Peaks_locations_monitored_reference_intensities[0])): #just using the first fragment's row to get the number of molecules.
        moleculeFragmentationPattern = ReferenceDataObject.monitored_reference_intensities[:,moleculeIndex]
        moleculeFragmentationPatternSignificantPeakLocations = signiciant_Peaks_locations_monitored_reference_intensities[:,moleculeIndex]
        #Below is a Boolean array with values of true anytime there was a negligble signal for a significant peak
        negligibleSignalForSignificantPeak = signalsNegligible*moleculeFragmentationPatternSignificantPeakLocations
        #the sum of the Boolean array will be > 0 if there was any cases with negligible signals for significant peaks.
        if sum(negligibleSignalForSignificantPeak) > 0:
            ReferenceDataObject.monitored_reference_intensities[:,moleculeIndex] = ReferenceDataObject.monitored_reference_intensities[:,moleculeIndex]*0 #note that this changes the reference object directly.
            ReferenceDataObject.reciprocal_matching_correction_values[:,moleculeIndex] = ReferenceDataObject.reciprocal_matching_correction_values[:,moleculeIndex]*0 #note that this changes the reference object directly.
            indexesOfMoleculesSetToZero.append(moleculeIndex)
    return ReferenceDataObject


    
#exports the user input file so that it can be used in the next iteration
def ExportUserInputFile(fileName, userInputModule):
    
    #Creating an updated UI file
    globalsFE_saveFile = fileName 
    globalsFE_loadFile = fileName 
    #create export object
    globalsFE_object = ei.module_export_import(globalsFE_saveFile,globalsFE_loadFile,userInputModule)
    
    #save variables to the text file 
    globalsFE_object.save_params()

    
#this function is used to append any list to a file in an executable fashion
def AppendListToFile(listVariableName, List, FileName, entriesPerLine=float('Inf')):
    #open the file in an append fashion
    with open(FileName,'a+') as f:
        #write in the variable name and open the list
        f.write('\n%s = [' %listVariableName)
        #find the length of the list
        length = len(List)
        #loop through every value in the list 
        for listIndex in range(length):
            #Construct a string for each value in turn
           string = "'%s'," % str(List[listIndex])
           #if this is the last entry in the list, remove the trailing comma
           if listIndex == length - 1:
               string = string[:-1]
           #write the entry into the file 
           f.write(string)
           if entriesPerLine != float('Inf'): #This is infinity in python.
               #if the previous entry was the last one in that line, add a newline character
               if (listIndex + 1) % entriesPerLine == 0: 
                   f.write('\\'+'\n') #use \ character so that if any code reads this file later it will know the list is not over.
        #write in the closing bracket for the list
        f.write(']')
    return None
    
def StringSearch(string, keyword = '', antikeyword = ''):
    if keyword in string and not antikeyword in string:
        return True
    else:
        return False

'''
This is a helper function that removes the _iter_ from fileNames so the program can access reference files and collected files from the parent directory when running iterative analysis
'''    
def remove_iter_fromFileName(fileName):
    startingIndexOfStringToRemove = fileName.find('_iter_')
    if '.csv' in fileName:
        endingIndexOfStringToRemove = fileName.find('.csv')
    if '.tsv' in fileName:
        endingIndexOfStringToRemove = fileName.find('.tsv')
    originalFileName = fileName[0:startingIndexOfStringToRemove] + fileName[endingIndexOfStringToRemove:len(fileName)]
    return originalFileName
    
#This supporting function of IterativeAnalysisPreProcessing finds the highest suffix of any file that contains a given keyword. 
def FindHighestDirNumber(keyword):
    listIterDirectories =[]
    #Search all files/folders in the current directory
    for directoryname in os.listdir():
        #if one of them contains the keyword and is a directory i.e. no '.'
        if StringSearch(directoryname, keyword, '.'):
            #append it to the list
            listIterDirectories.append(directoryname)
        
    suffixlist = []
    #for all files/folders with the right keyword
    for directoryname in listIterDirectories:
        #append the last value of each name
        suffixlist.append(directoryname[-1])
    #return the highest of the last values
    if not len(suffixlist) == 0:
        return(max(suffixlist))
    return 1

#This supporting function of IterativeAnalysisPreProcessing confirms that a directory exists
def EnsureDirectory(dir_path):
    directory = dir_path
    #isolate the directory name
    #this line can be used if needed
    #directory = os.path.dirname(dir_path)
    #if the directory doesn't already exist
    if not os.path.exists(directory):
        #create the directory
        os.makedirs(directory)

def SpecificIterationName(iterativeAnalysis, iterationNumber):
    # Update: Clint, Sept 28 2018, this now returns an absolute path
    # Previously a relative path was returned but it was in windows format
    # i.e. ".\". That caused problems for linux.
    #if the user has entered an iteration name
    if iterativeAnalysis == False or iterativeAnalysis == True:
         #create the default directory
         iterationDirectoryName = os.path.join(
             os.curdir, "_iter_{}".format(str(iterationNumber)))
    else:
        #set that name to be the directory along with the correct number 
        iterationDirectoryName = os.path.join(os.curdir,
            str(iterativeAnalysis),
             "_iter_{}".format(str(iterationNumber)))
    return iterationDirectoryName
    
def IterativeDirectoryChange(iterativeAnalysis, iterationNumber):
    iterationDirectoryName = SpecificIterationName(iterativeAnalysis, iterationNumber)
    #confirm that the directory exists
    EnsureDirectory(iterationDirectoryName)
    #Change the working directory to the new directory name. 
    #'THIS IS A HIGHLY SIGNIFICANT LINE, because it redirects all of the output for the rest of the program run'
    os.chdir(iterationDirectoryName) #this will be changed back at the end of the program    

def IterationDirectoryPreparation(iterativeAnalysis, iterationNumber):
    #implied arguments for this function are G.referenceFileNamesList and G.collectedFileName
    IterativeDirectoryChange(iterativeAnalysis, iterationNumber)
    
    #naming for collected file
    #record the old file names 
    G.oldcollectedFileName = G.collectedFileName
    #construct the file names for the current run of the program
    #TODO FIXME, This syntax with -21 will not allow iterative to be compatible with more than 9 iterations
    collectedFileNameTemp = str(G.collectedFileName)[:-21] +  str(G.iterationSuffix) + str(G.collectedFileName)[-4:]      
    #copy the experimental and reference files into new names for this iterative run
    shutil.copy(G.collectedFileName, collectedFileNameTemp)
    
    #change the globals to reflect the renaming of the ref and exp files
    G.collectedFileName =  collectedFileNameTemp
    
    #construct file names for the next run of the program 
    #TODO FIXME, This syntax with -11 will not allow iterative to be compatible with more than 9 iterations
    G.nextExpFileName = G.collectedFileName[:-11] +  str('_remaining') + G.collectedFileName[-11:]
    #naming for reference files
    
    G.oldReferenceFileName = []
    G.nextRefFileName = []
    for RefIndex, RefName in enumerate(G.referenceFileNamesList): #a list
        #record the old file names 
        G.oldReferenceFileName.append(RefName)
        
        #construct the file names for the current run of the program
        #TODO FIXME, This syntax with -18 will not allow iterative to be compatible with more than 9 iterations
        referenceFileNameTemp = G.referenceFileNamesList[RefIndex][:-18] +  str(G.iterationSuffix) + G.referenceFileNamesList[RefIndex][-4:]
        
        #copy the experimental and reference files into new names for this iterative run
        shutil.copy(RefName, referenceFileNameTemp)
        
        #change the globals to reflect the renaming of the ref and exp files
        G.referenceFileNamesList[RefIndex] =  referenceFileNameTemp
        
        #construct file names for the next run of the program 
        #TODO FIXME, This syntax with -18 will not allow iterative to be compatible with more than 9 iterations
        G.nextRefFileName.append(RefName[:-18] + '_unused_iter_%s' %G.iterationNumber + RefName[-4:])
    return None
    #implied returns: G.oldReferenceFileName, G.oldcollectedFileName, G.referenceFileNamesList,G.collectedFileName, G.nextRefFileName, G. nextExpFileName, G.iterationNumber 

def IterationFirstDirectoryPreparation(iterativeAnalysis,iterationNumber):
    #implied arguments for this function are G.referenceFileNamesList and G.collectedFileName
    #this global value is set so that each export statement can label the output files correctly
    G.iterationNumber = iterationNumber
    
    iterationDirectoryName = SpecificIterationName(iterativeAnalysis, iterationNumber)
       
    #confirm that the directory exists
    EnsureDirectory(iterationDirectoryName)
    
    #Change the working directory to the new directory name. 
    'THIS IS A HIGHLY SIGNIFICANT LINE, because it redirects all of the output for the rest of the program run'
    os.chdir(iterationDirectoryName) #this will be changed back at the end of the program
    
    #copy the first UserInputFile into the first iteration directory
    ExportUserInputFile("UserInput_iter_1.py", G)
    #append the variable list to the user input file
    AppendListToFile("__var_list__", G.__var_list__, "UserInput_iter_1.py", float('Inf'))
    
    #record the old file names 
    G.oldcollectedFileName = G.collectedFileName
    #construct the file names for the first run of the program
    G.collectedFileName = G.collectedFileName[:-4] +  str(G.iterationSuffix) + G.collectedFileName[-4:]
    #construct file names for the second run of the program 
    #TODO FIXME, This syntax with -11 will not allow iterative to be compatible with more than 9 iterations
    G.nextExpFileName = G.collectedFileName[:-11] + '_remaining_iter_1' + G.collectedFileName[-4:]
    
    G.oldReferenceFileName = []
    for RefIndex, RefName in enumerate(G.referenceFileNamesList): #a list
        G.oldReferenceFileName.append(RefName)
        #construct the file names for the first run of the program
        G.referenceFileNamesList[RefIndex] = G.referenceFileNamesList[RefIndex][:-4] +  str(G.iterationSuffix) + G.referenceFileNamesList[RefIndex][-4:]
        #construct file names for the second run of the program 
        G.nextRefFileName.append(RefName[:-4] + '_unused_iter_1' + RefName[-4:])
    
    return None 
    #implied returns: G.oldReferenceFileName, G.oldcollectedFileName, G.referenceFileNamesList,G.collectedFileName, G.nextRefFileName, G. nextExpFileName, G.iterationNumber 

#The IterativeAnalysisDirectory and Variable Population function is used to shrink the size of the program analysis and redirect the output. 
def IADirandVarPopulation(iterativeAnalysis, chosenMassFragments, chosenMolecules, ExperimentData, ExperimentDataFullCopy, ReferenceDataList, ReferenceDataListFullCopy):
    #implied arguments: G.dataSimulation, G.referenceFileNamesList, G.collectedFileName, G.nextRefFileName, G.oldReferenceFileName, G.chosenMoleculesNames, G.iterationNumber
    #override data simulation to yes if it was not selected
    if G.dataSimulation != 'yes':
        print("Iterative analysis cannot find the remaining signals in the experiment without signal simulation being run.")
        print("User selection to skip signal simulation has been overridden. ")
        G.dataSimulation = 'yes'
    #Warn the user if they are trying to run an iteration that has no molecules to solve. (This would cause a complex error further on in the program if allowed to run.)    
    if len(ReferenceDataList[0].molecules) == 0:
        print("Warning Message: There are inadequate molecules to perform another iteration. Please confirm that there are still remaining molecules to solve.")
        sys.exit()

    ReferenceDataSS = []
    ReferenceDataSSreciprocal_matching_correction_valuesList = [] #This list will hold the correction factors used for each individual reference pattern
    for RefObjectIndex, RefObject in enumerate(ReferenceDataList): #a list
        #Creating a correction values matrix list for signal simulation at the end of the program
        ReferenceDataSS.append(copy.deepcopy(RefObject))
        ReferenceDataSS[RefObjectIndex] = ReferenceInputPreProcessing(ReferenceDataSS[RefObjectIndex])
        ReferenceDataSS[RefObjectIndex] = Populate_reciprocal_matching_correction_values(ExperimentDataFullCopy.mass_fragment_numbers,ReferenceDataSS[RefObjectIndex])
        ReferenceDataSSreciprocal_matching_correction_valuesList.append(ReferenceDataSS[RefObjectIndex].reciprocal_matching_correction_values)
        RefObject.SSreciprocal_matching_correction_values = ReferenceDataSS[RefObjectIndex].reciprocal_matching_correction_values
        
    #Selecting unused Reference Data
    unusedMolecules = []
    for molecule in ReferenceDataListFullCopy[0].molecules:
        if not molecule in G.chosenMoleculesNames:
            unusedMolecules.append(molecule)
    
    for RefObjectIndex, RefObject in enumerate(ReferenceDataList): #a list
        #Export current Reference Data  
        #Reference data is trimmed prior to this function
        ExportXYYYData(G.referenceFileNamesList[RefObjectIndex], RefObject.provided_reference_patterns, RefObject.molecules, abscissaHeader = 'M/Z')
    
    #Export current Experimental Data
    #Experimental data is trimmed prior to this function, but it still needs to be exported  
    ExportXYYYData(G.collectedFileName, ExperimentData.workingData, ExperimentData.mass_fragment_numbers,
              abscissaHeader = ExperimentData.abscissaHeader, dataType = 'preProcessed', rowIndex = ExperimentData.times)
   
    for RefObjectIndex, RefObject in enumerate(ReferenceDataList): #a list
        #export reference data for next iteration
        if G.iterationNumber == 1: #first iteration files aren't in standard locations
            referenceFilePath = os.path.normpath(
                os.path.join(os.curdir,
                    os.pardir,
                    str(G.oldReferenceFileName[RefObjectIndex])))
            referenceDataToExport = DataFunctions.TrimReferenceFileByMolecules(unusedMolecules,
                referenceFilePath)
            referenceDataToExport.to_csv(G.nextRefFileName[RefObjectIndex], header = False, index = False)
        else: #not first iteration
        #generate unused reference data
            referenceDataToExport = DataFunctions.TrimReferenceFileByMolecules(unusedMolecules, G.oldReferenceFileName[RefObjectIndex])
            
            referenceDataToExport.to_csv(G.nextRefFileName[RefObjectIndex], header = False, index = False)
    
    return ReferenceDataSSreciprocal_matching_correction_valuesList, unusedMolecules

#sortArrayByHeader is used in exportSimulatedSignalsSoFar to sort simulated data in order by mass fragment
def sortArrayByHeader(headerArray,dataArray):
    for headerIndex in range(len(headerArray)): #Loop through header
        headerArray[headerIndex] = headerArray[headerIndex][1:] #remove the 'm' on each string
    headerArray = headerArray.astype(float) #convert array to array of floats for sorting
    indiciesArray = numpy.argsort(headerArray) #return an array of indicies determined by sorting
    sortedArray = copy.copy(dataArray) #create a copy of the data array
    sortedHeader = copy.copy(headerArray)
    
    for index in range(0,len(headerArray)): #loop through the header array
        sortedHeader[index] = headerArray[indiciesArray[index]]
        sortedArray[:,index] = dataArray[:,indiciesArray[index]]
    
    sortedHeader = sortedHeader.astype(str)
    for headerIndex in range(len(headerArray)):
        sortedHeader[headerIndex] = 'm' + sortedHeader[headerIndex]
        
    return sortedArray, sortedHeader

def exportSimulatedSignalsSoFar(simulatedSignalsOutputName,iterationDirectoryName,iterationNumber):
    #read the simulated raw signals from the iteration's csv file
    simulatedSignalsOutputNameIterative = simulatedSignalsOutputName[:-4] + iterationDirectoryName + simulatedSignalsOutputName[-4:] #get the filename with the iterative suffix on it
    simulatedSignalsFromIterativeAbsoluteFileName = os.path.join(os.getcwd(),iterationDirectoryName,simulatedSignalsOutputNameIterative) #get the absolute path
    
    simulatedSignalsFromIterative = numpy.genfromtxt(simulatedSignalsFromIterativeAbsoluteFileName,delimiter=',',dtype=None, encoding='latin1').astype(str) #Get the data from the iteration's simultaed raw signals file (use astype(str) to keep b' from showing up before each entry)

    if iterationNumber == 1: #if in the first iteration we do not have SimulatedRawSignalsSoFarIterative.csv so we create it by saving the first iteration's simulated signals to it
        #Create the array to export
        numpy.savetxt('SimulatedRawSignalsSoFarIterative.csv',simulatedSignalsFromIterative,fmt='%s',delimiter=',', encoding='latin1')
    elif iterationNumber != 1: #We are at a later iteration so we need to read simulatedRawSignalsSoFarIterative.csv to see what signals are there
        simulatedSignalsSoFar = numpy.genfromtxt('SimulatedRawSignalsSoFarIterative.csv',delimiter=',',dtype=None, encoding='latin1') #read the simulatedRawSignalsSoFarIterative.csv file
        #Break up the simulated signals so far data into arrays
        SS_SoFarHeader = simulatedSignalsSoFar[0,1:].astype(str) #get the header
        SS_SoFarAbscissa = simulatedSignalsSoFar[1:,0].astype(float) #Get the abscissa
        SS_SoFarData = simulatedSignalsSoFar[1:,1:].astype(float) #get the data
        SS_SoFarAbscissaHeader = simulatedSignalsSoFar[0,0] #get the abscissa header
        
        #Similarly break up the iteration's simulated signals into arrays
        SS_IterativeHeader = simulatedSignalsFromIterative[0,1:].astype(str) #get the header
        SS_IterativeAbscissa = simulatedSignalsFromIterative[1:,0].astype(float) #Get the abscissa
        SS_IterativeData = simulatedSignalsFromIterative[1:,1:].astype(float) #get the data
        SS_IterativeAbscissaHeader = simulatedSignalsFromIterative[0,0] #get the abscissa header
        SS_SoFarAbscissa = numpy.reshape(SS_SoFarAbscissa,(len(SS_SoFarAbscissa),1)) #reshape the abscissa so it's at least in 2d with only 1 column
        
        #Loop through mass fragments in the current in this iteration and see if it is in simulatedRawSignalsSoFarIterative
        for massFragmentIndex in range(len(SS_IterativeHeader)):
            if SS_IterativeHeader[massFragmentIndex] in SS_SoFarHeader: #if the mass fragment is in both simulatedSignalsSoFar and the current iteration's simulated signals, then we need to add the columns
                currentMFIndex = numpy.where(numpy.array(SS_IterativeHeader[massFragmentIndex]) == numpy.array(SS_SoFarHeader))[0][0] #Get the column index in SS_SoFar, use [0][0] since numpy.where returns an an array in an array
                SS_SoFarData[:,currentMFIndex] = SS_SoFarData[:,currentMFIndex] + SS_IterativeData[:,massFragmentIndex] #add the signals to the appropriate column
            else: #otherwise the mass fragment is not in simulatedSignalsSoFar and needs to be appended to the data
                SS_SoFarHeader = numpy.append(SS_SoFarHeader,SS_IterativeHeader[massFragmentIndex]) #append the mass fragment to the header
                columnToAppend = numpy.atleast_2d(SS_IterativeData[:,massFragmentIndex])
                columnToAppend = columnToAppend.transpose()
                SS_SoFarData = numpy.hstack((SS_SoFarData,columnToAppend))
        sortedData, sortedHeaders = sortArrayByHeader(SS_SoFarHeader,SS_SoFarData) #Sort the data by the mass fragment number
        stackedData = numpy.vstack((sortedHeaders,sortedData)) #stack the sorted headers on the sorted data
        stackedAbscissa = numpy.vstack((SS_SoFarAbscissaHeader,SS_SoFarAbscissa)) #stack the abscissa header with the abscissa
        simulatedDataToExport = numpy.hstack((stackedAbscissa,stackedData))
        numpy.savetxt('SimulatedRawSignalsSoFarIterative.csv',simulatedDataToExport,fmt='%s',delimiter=',', encoding='latin1')

def IterativePrepareNextIterationInputFiles(ExperimentDataFullCopy):
    #Has implied argument of Globals file G.
    
    #update the suffix number and create the next user input file
    G.nextIterationSuffix = '_iter_%s' %str(G.iterationNumber + 1)
    G.nextIterationNumber = G.iterationNumber + 1
    nextUserInputFileName = 'UserInput%s.py' %G.nextIterationSuffix 
  
    #create the next iteration directory and change the cwd into the next iteration directory
    IterativeDirectoryChange(G.iterativeAnalysis, G.nextIterationNumber)
    #export the current input specifications, *specifically* so that we can get the new userinput module ready and start populating it without overwriting the origina userinput module.
    ExportUserInputFile(nextUserInputFileName, G)
    G.nextUserInputModule = False #first, just this equal to boolean because we can't do a deep copy. Annoyingly, I only got this to work by making the variable nextUserInputModule inside a globals module.
    sys.path.append(os.getcwd())  #Surprisingly, I was not able to import the module in the next line unless I added the current working directory to the system path.
    exec('import ' + nextUserInputFileName[:-3])
    exec('G.nextUserInputModule = '+ nextUserInputFileName[:-3])
    #Now G.nextUserInputModule contains the next user input module (currently like a draft module), which we can populate below and then export (writing over the draft file).

    G.nextUserInputModule.__var_list__ = G.__var_list__ 
    #save the new file name for the next user input file 
    G.nextUserInputModule.collectedFileName = G.nextExpFileName 
    G.nextUserInputModule.referenceFileNamesList = G.nextRefFileName
    #updating the selected molecules for the next user input file
    G.nextUserInputModule.chosenMoleculesNames = G.unusedMolecules
    #Updating the selected masses for the next user input file
    chosenMasses = list(ExperimentDataFullCopy.mass_fragment_numbers)
    G.nextUserInputModule.chosenMassFragments = [int(x) for x in chosenMasses]
    #Turn off Data Smoothing 
    G.nextUserInputModule.dataSmootherYorN = 'no'
    #Turn off manual baseline correction
    G.nextUserInputModule.backgroundMassFragment = []
    G.nextUserInputModule.backgroundSlopes = []
    G.nextUserInputModule.backgroundIntercepts = []
    #Turn off semiauto baseline correction
    G.nextUserInputModule.linearBaselineCorrectionSemiAutomatic = 'no'
    #turnoff reference pattern extractor.
    G.nextUserInputModule.extractReferencePatternFromDataOption = 'no'
    #turnoff any scaling of raw data since that should not occur after first iteration.
    G.nextUserInputModule.scaleRawDataOption = 'manual'
    G.nextUserInputModule.scaleRawDataFactor = G.scaleRawDataFactor
    
    #Now going to overwrite parallelized variables with their original versions if they were set to length of chosen molecules.
    delimitedStringOfVariablesToUnparallelize = 'moleculeLikelihoods, sensitivityValues, referenceValueThreshold, referenceSignificantFragmentThresholds'
    listOfVariablesToUnparallelize = delimitedStringOfVariablesToUnparallelize.split(", ") #Note that we are using ", " as the delimiter, not just ","
    for variable in listOfVariablesToUnparallelize:
        G.nextUserInputModule.__dict__[variable]=G.beforeParsedGDict[variable]
    
    #export the user input specifications 
    ExportUserInputFile(nextUserInputFileName, G.nextUserInputModule)
    #append the variable list to the user input file
    AppendListToFile("__var_list__", G.nextUserInputModule.__var_list__, nextUserInputFileName, float('Inf'))
    
    if G.nextUserInputModule.iterativeAnalysis == True:
        iterationDirectoryName = '_iter_%s' %(str(G.nextUserInputModule.iterationNumber))
    if not G.nextUserInputModule.iterativeAnalysis == True:
        iterationDirectoryName = '%s_iter_%s' %(G.nextUserInputModule.iterativeAnalysis, str(G.nextUserInputModule.iterationNumber))
    #copy the experimental signals to the next iteration
    copyFromPath = os.path.join(os.curdir, os.pardir,
            str(iterationDirectoryName),
            str(G.nextUserInputModule.collectedFileName))
    shutil.copy(copyFromPath, os.getcwd())
    for RefIndex, RefName in enumerate(G.nextUserInputModule.referenceFileNamesList): #a list
        #copy the next reference file from the previous iteration folder to the next iteration folder
        copyFromPath = os.path.join(os.curdir,
            os.pardir,
            str(iterationDirectoryName),
            str(G.nextUserInputModule.referenceFileNamesList[RefIndex]))
        shutil.copy(copyFromPath, os.getcwd())
    #returning to the parent directory
    os.chdir('..')

def IterativeAnalysisPostProcessing(ExperimentData, simulateddata, mass_fragment_numbers,ExperimentDataFullCopy, times, concdata, molecules):
    #implied arguments: G.iterationSuffix, G.nextRefFileName, G.nextExpFileName, G.iterativeAnalysis, G.unusedMolecules, G.iterationNumber
    #remove the signals that have already been solved for from the data set
    ExperimentData.workingData = DataFunctions.RemoveSignals(ExperimentDataFullCopy.workingData, ExperimentDataFullCopy.mass_fragment_numbers, simulateddata, mass_fragment_numbers)
    
    #Export the remaining experimental signals
    ExportXYYYData(G.nextExpFileName, ExperimentDataFullCopy.workingData, ExperimentDataFullCopy.mass_fragment_numbers, abscissaHeader = ExperimentData.abscissaHeader, dataType = 'Experiment', rowIndex = ExperimentData.times)
    #revert to the parent directory
    os.chdir('..')

    #Adding the Additional concentrations to the overall concentration results
    moleculeConcLabels = ['%s Concentration Relative to CO' % molecule for molecule in molecules] 
    DataFunctions.AppendColumnsToCSV(G.TotalConcentrationsOutputName, concdata, moleculeConcLabels, times, ["Time"])
    
    if G.iterativeAnalysis == True:
        iterationDirectoryName = '_iter_%s' %(str(G.iterationNumber))
    if not G.iterativeAnalysis == True:
        iterationDirectoryName = '%s_iter_%s' %(G.iterativeAnalysis, str(G.iterationNumber))

    #This function will take in the iteration directory name, iteration number, ExperimentDataFullCopy.abscissaHeader,simulated data, mass fragment numbers, and the times array and will add or append the simulated signals to a file called SimulatedRawSignalsSoFarIterative.csv
    exportSimulatedSignalsSoFar(G.simulatedSignalsOutputName,iterationDirectoryName,G.iterationNumber) #subtract 1 from the iteration number since the iteration number has already been changed
    IterativePrepareNextIterationInputFiles(ExperimentDataFullCopy)
    return None
     #implied returns: G.referenceFileNamesList, G.collectedFileName, G.nextRefFileName, G.chosenMoleculesNames, G.iterationSuffix
###############################################################################
#########################  Functions to read data files #######################
###############################################################################
#These functions read in the experimental data file and the reference file. The
#returned variables can then be used to initialize the respective classes.

#a small helper function to check if an extension exists in a filename and to return the delimiter based on that.
def getDelimiterFromExtension(filename):
    if ".tsv" in filename:
        delimiter = '\t'
    elif ".tab" in filename:
        delimiter = '\t'
    elif ".txt" in filename:
        delimiter = '\t'        
    elif ".skv" in filename:
        delimiter = ';'
    elif ".csv" in filename:
        delimiter = ',' #it could be something else, but we will assume that a csv
    else:
        delimiter = '\t' #for MSRESOLVE, this is now the default delimiter.
    return delimiter
        
def readDataFile(collectedFileName):

 #read the csv file into a dataframe.  dataFrame means "dataframe" and is a pandas object.
    dataFrame = pandas.read_csv('%s' %collectedFileName, header=None)
    ''' generate mass fragment list'''
    #select only the 2nd row down, all columns except for the first. 
		#"iloc" is a pandas dataframe function. All it does is select a portion of the data.
    dfmass = dataFrame.iloc[1][1:]
    #convert to matrix
    masses = dfmass.values
    #sort through the matrix and remove labels
    for i in range(0,len(masses)):
        if type(masses[i]) == type("string"):
            masses[i] = masses[i].replace('mass','')
            masses[i] = masses[i].replace('m','')
    #convert the matrix to doubles
    mass_fragment_numbers = masses.astype(numpy.double)
            
    '''generate time list'''
    # set abscissa header (time or Temp, etc.)
    abscissaHeader = dataFrame.iloc[1][0]
    #select column of times
    dftimes = dataFrame.iloc[2:][0]
    #convert to matrix
    times = dftimes.values
    #save as class object with type float
    times = times.astype(float)
    #if the user wants to analyze one point, the data is doubled in length
    #to prevent future index problems
    if len(times) == 1:
        times = numpy.append(times,times)
        times[1] = times[0]*1.1
   
    '''generate collected data'''
    #select matrix of raw signals
    dfcollected = dataFrame.iloc[2:,1:]
    #convert to matrix
    collected = dfcollected.values
    #save as class object with type float
    rawCollectedData = collected.astype(float)
    #if the user wants to analyze one point, the data is doubled in length
    #to prevent future index problems
    if len(rawCollectedData) == 1:
        rawCollectedData = numpy.vstack((rawCollectedData,rawCollectedData))

        
    return mass_fragment_numbers, abscissaHeader, times, rawCollectedData, collectedFileName

#readReferenceFile is a helper function that reads the reference file in a certain form and returns the
#variables and data that are used to initialize the class. It can read files both in XYYY and XYXY form.
def readReferenceFile(referenceFileName, form):        
    #This function converts the XYXY data to an XYYY format
    form = form.lower()
    def FromXYXYtoXYYY(provided_reference_patterns):
        #print("Warning: FromXYXYtoXYYY for converting data patterns has not been tested in a long time. A unit test should be created and checked prior to use. Then this warning updated (this warning appears in two parts of the code." )
        masslists = [] #future lists must be must empty here to append in the for loops
        relativeintensitieslists = [] #future list
        #this loops gathers all the mass fragment numbers for each molecule in one list of arrays, while a second
        #list is made, gathering the relative intensities so that they were indexed the same as their mass fragment
        #numbers in the other list
        #this for loop grabs stuff from the reference array, whose orientation and identity is shown in the flow chart arrays document
        for referenceBy2Index in range(0,len(provided_reference_patterns[0,:]),2):#array-indexed for loop, only gets every other value, as half the indexes are mass lists, and the other half are relative intensity
            masslists.append(provided_reference_patterns[:,referenceBy2Index])#these are lists of arrays
            relativeintensitieslists.append(provided_reference_patterns[:,referenceBy2Index+1])#the relative intensities are after every counter, so there is a +1 (it is array indexed so since the first column is a mass list all the +1's are relative intensities)
        masslist = [] #future list
        #This for loop gets all of the mass fragments from the first index of the list, basically by not adding the 
        #'nan's or empty spaces after the numbers
        provided_mass_fragments = provided_reference_patterns[:,0]
        for referenceIndex in range(len(provided_mass_fragments)): #array-indexed for loop
            if str(masslists[0][referenceIndex]) != 'nan': #we do not want nan's in our array, the genfromtxt function calls empty boxes in excel (might be in .csv as well)'nan'.
                masslist.append(masslists[0][referenceIndex])
        #this second nested for loop gathers all the other mass fragment numbers that have not already been added to
        #the masslist, basically obtaining all the masses in the reference data and then after the loop they are sorted
        #using .sort, then an empty array of zeros is made to accommodate the output array
        for masslistIndex in range(1,len(masslists)):#array-indexed for loop, starts at one because it's checking against all arrays besides itself
            for referenceIndex in range(len(provided_mass_fragments)):#array-indexed for loop
                if str(masslists[masslistIndex][referenceIndex]) != 'nan':
                    if sum(masslists[masslistIndex][referenceIndex] == numpy.array(masslist)) == 0:#if the value being looked at is not equal to anything in our masslist already
                        masslist.append(masslists[masslistIndex][referenceIndex])
        masslist.sort()#puts the list in order
        reference_holder = numpy.zeros([len(masslist),int(len(provided_reference_patterns[0,:])/2)+1])#makes an array that is full of zeros to hold our future reference array
        reference_holder[:,0:1] = numpy.vstack(numpy.array(masslist))#This puts the mass list in the first column of our new reference array
        #Finally, the for loop below makes a list each revolution, comparing each list of mass fragments (for each molecule)
        #and adding the relative intensities (from the identically indexed array) when they numbers were equal, and otherwise
        #adding a zero in its place. It then adds this list to the array (using numpy.vstack and numpy.array)
        for massListsIndex in range(len(masslists)):#array-indexed for loop
            relativeintensitieslist = [] #empties the list every loop
            for massListIndex in range(len(masslist)):
                placeholder = 0 #after the next for loop finishes, this is reset
                for specificMassListIndex in range(len(masslists[massListsIndex])):#array-indexed for loop, each column of .csv file being checked
                    if masslists[massListsIndex][specificMassListIndex] == masslist[massListIndex]:#This is when the index for the correct mass fragment is found
                        relativeintensitieslist.append(relativeintensitieslists[massListsIndex][specificMassListIndex])#relative intensities lists index the same way
                        placeholder = 1 #so that the next if statement will know that this has happened
                    if specificMassListIndex == len(masslists[massListsIndex])-1 and placeholder == 0:#If it comes to the end of the for loop, and there's no match, then the relative intensity is zero
                        relativeintensitieslist.append(0)
                if massListIndex == len(masslist)-1:#Once the larger for loop is done the 
                    reference_holder[:,(massListsIndex+1):(massListsIndex+2)] = numpy.vstack(numpy.array(relativeintensitieslist)) #the list is made into an array and then stacked (transposed)
        provided_reference_patterns = reference_holder
        return provided_reference_patterns
    
    #read the csv file into a dataframe
    if '.csv' in referenceFileName:
        dataFrame = pandas.read_csv('%s' %referenceFileName, header = None)
    elif '.tsv' in referenceFileName:
        try: #no easy way to assess utf16 vs utf8, so try both.
            dataFrame = pandas.read_csv('%s' %referenceFileName, header = None, delimiter = '\t', encoding = 'utf8') #need to specify encoding for cases of tab delimited files.
        except: #no easy way to assess utf16 vs utf8, so try both.
            dataFrame = pandas.read_csv('%s' %referenceFileName, header = None, delimiter = '\t', encoding = 'utf16') #need to use utf16 for some cases of tab delimited files.
    
    if form == 'xyyy':
        for rowIndex in range(len(dataFrame)): #Loop through each row and check the abscissa value
            try: #Try to convert the abscissa title to a float
                float(dataFrame.iloc[rowIndex][0]) #if successful, then this rowIndex is the first index of provided reference intensities
                dfreference = dataFrame.iloc[rowIndex:][:] #remove the rows of headers
                reference = dfreference.values #convert to matrix
                provided_reference_patterns = reference.astype(float) #convert the matrix to floats
                provided_reference_patterns = DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(provided_reference_patterns,startingRowIndex=1) #clear row of zeros
                break #exit the for loop since the first non-header row has been reached.
            except: #Otherwise the row consists of other information
                if (dataFrame.iloc[rowIndex][0] == 'SourceOfFragmentationPatterns') or (dataFrame.iloc[rowIndex][0] == 'Source:'): #if the abscissa titles the source (both old and new reference files)
                    dfSourceOfFragmentationPatterns = dataFrame.iloc[rowIndex][1:] #select the row of names
                    SourceOfFragmentationPatterns = dfSourceOfFragmentationPatterns.values #convert to matrix
                    SourceOfFragmentationPatterns = SourceOfFragmentationPatterns.astype(str) #save as class object with type string
                elif dataFrame.iloc[rowIndex][0] == 'sourceOfIonizationData':
                    dfsourceOfIonizationData = dataFrame.iloc[rowIndex][1:] #Select the row of names
                    sourceOfIonizationData = dfsourceOfIonizationData.values #convert to matrix
                    sourceOfIonizationData = sourceOfIonizationData.astype(str) #save as class object with type string
                elif dataFrame.iloc[rowIndex][0] == 'Molecules': #if the abscissa titles the molecule names
                    dfmolecules = dataFrame.iloc[rowIndex][1:] #select the row of names
                    molecules = dfmolecules.values #convert to matrix
                    molecules = molecules.astype(str) #save as class object with type string
                    molecules = list(molecules)
                    for moleculeIndex in range(len(molecules)):
                        molecules[moleculeIndex] = molecules[moleculeIndex].strip()#remove leading and trailing whitespaces.
                elif dataFrame.iloc[rowIndex][0] == 'Electron Numbers': #if the abscissa titles the electron numbers
                    dfelectronnumbers = dataFrame.iloc[rowIndex][1:] #select the row of names
                    electronnumbers = dfelectronnumbers.values #convert to matrix
                    electronnumbers = electronnumbers.astype(int) #save as class object with type int
                elif dataFrame.iloc[rowIndex][0] == 'Molecular Mass': #if the abscissa titles the molecular weights
                    dfmolecularWeights = dataFrame.iloc[rowIndex][1:] #select row of names
                    molecularWeights = dfmolecularWeights.values #convert to matrix
                    molecularWeights = molecularWeights.astype(numpy.double) #save as class object with type double
                elif dataFrame.iloc[rowIndex][0] == 'moleculeIonizationType':
                    dfmoleculeIonizationType = dataFrame.iloc[rowIndex][1:] #select row of names
                    moleculeIonizationType = dfmoleculeIonizationType.values #convert to matrix
                    moleculeIonizationType = moleculeIonizationType.astype(str) #save as class object with type string
                elif (dataFrame.iloc[rowIndex][0] == 'relativeIonizationEfficiencies') or (dataFrame.iloc[rowIndex][0] == 'knownIonizationFactorsRelativeToN2'):
                    dfrelativeIonizationEfficiencies = dataFrame.iloc[rowIndex][1:] #select row of names
                    relativeIonizationEfficiencies = dfrelativeIonizationEfficiencies.values #convert to matrix
                    for index in range(len(relativeIonizationEfficiencies)):
                        try: #try to convert to a float
                            relativeIonizationEfficiencies[index] = float(relativeIonizationEfficiencies[index])
                        except: #if not possible, the value is probably None or 'unknown' so leave as a string
                            pass
#                    relativeIonizationEfficiencies = relativeIonizationEfficiencies.astype(float) #save as class object with type float

        try: #if using an older reference file, it will not have ionization factors so the elif statement never gets entered meaning knownIonizationFactors does not exist
            relativeIonizationEfficiencies #Try calling this variable, if it exists there will be no error
        except: #if it does not exist, populate it with unknown
            relativeIonizationEfficiencies = ['unknown'] #initialize as a list of len(1)
            relativeIonizationEfficiencies = relativeIonizationEfficiencies*len(molecules)   #duplicate the list for # of molecules.
            relativeIonizationEfficiencies = numpy.array(relativeIonizationEfficiencies) #convert to matrix
            
        try: #if using an old reference file, it will not have ionization types so the elif statement never gets entered meaning moleculeIonizationType does not exist
            moleculeIonizationType #Try calling this variable, if it exists there will be no error
        except: #if it does not exist, populate it with unknown
            moleculeIonizationType = ['unknown'] #initialize as a list of len(1)
            moleculeIonizationType = moleculeIonizationType *len(molecules)  #duplicate the list for # of molecules.
            moleculeIonizationType = numpy.array(moleculeIonizationType) #convert to matrix
            
        try: #If using an older reference file, it will not have SourceOfIonizationInfo so the elif statement never gets entered meaning the variable does not exist
            sourceOfIonizationData #try calling the variable, if it exists there will be no error
        except: #If it does not exist, populate with empty strings
            sourceOfIonizationData = ['unknown'] #initialize as a list of len(1)
            sourceOfIonizationData = sourceOfIonizationData * len(molecules) #duplicate the list for # of molecules.
            sourceOfIonizationData = numpy.array(sourceOfIonizationData) #convert to matrix
          
                
#        ''' generate reference matrix'''
#        #remove top 4 rows
#        dfreference = dataFrame.iloc[4:][:]
#        #convert to matrix
#        reference = dfreference.values
#        #convert the matrix to floats
#        provided_reference_patterns = reference.astype(float)
#        #clear rows of zeros
#        provided_reference_patterns=DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(provided_reference_patterns,startingRowIndex=1)
#    
#        '''generate electron number list'''
#        #select row of electron numbers
#        dfelectronnumbers = dataFrame.iloc[2][1:]
#        #convert to matrix
#        electronnumbers = dfelectronnumbers.values
#        #save as class object with type int
#        electronnumbers = electronnumbers.astype(int)
#   
#        '''generate list of molecule names'''
#        #select row of names
#        dfmolecules = dataFrame.iloc[1][1:]
#        #convert to matrix
#        molecules = dfmolecules.values
#        #save as class object with type string
#        molecules = molecules.astype(str)
#        
#        '''generate list of molecular weights'''
#        #select row of names
#        dfmolecularWeights = dataFrame.iloc[3][1:]
#        #convert to matrix
#        molecularWeights = dfmolecularWeights.values
#        #save as class object with type float
#        molecularWeights = molecularWeights.astype(float)
#        
#        '''generate list of source information'''
#        #select row of names
#        dfsourceInfo = dataFrame.iloc[0][1:]
#        #convert to matrix
#        sourceInfo = dfsourceInfo.values
#        #save as class object with type string
#        sourceInfo = sourceInfo.astype(str)
        
        '''list of massfragments monitored is not part of reference file'''
        mass_fragment_numbers_monitored = None
        
    elif form == 'xyxy':
        breakLoop = False 
        for rowIndex in range(len(dataFrame)): #Loop through each row and check the abscissa value
            try: #rowIndex == 4: #Try to convert the abscissa title to a float
                float(dataFrame.iloc[rowIndex][0]) #if successful, then this rowIndex is the first index of provided reference intensities
                dfreference = dataFrame.iloc[rowIndex:][:] #remove the rows of headers
                reference = dfreference.values #convert to matrix
                provided_reference_patterns = reference.astype(float) #convert the matrix to floats
                #print("Warning: FromXYXYtoXYYY for converting data patterns has not been tested in a long time. A unit test should be created and checked prior to use. Then this warning updated (this warning appears in two parts of the code." )
                provided_reference_patterns = FromXYXYtoXYYY(provided_reference_patterns) #convert reference from XYXY to XYYY
                provided_reference_patterns = DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(provided_reference_patterns,startingRowIndex=1) #clear row of zeros
                break #exit the for loop
                
            except: #rowIndex < 4: #Otherwise the row consists of other information
                if (dataFrame.iloc[rowIndex][0] == 'SourceOfFragmentationPatterns') or (dataFrame.iloc[rowIndex][0] == 'Source:'): #if the abscissa titles the source (both old and new reference files)
                    dfSourceOfFragmentationPatterns = dataFrame.iloc[rowIndex][1::2] #select the row of names
                    SourceOfFragmentationPatterns = dfSourceOfFragmentationPatterns.values #convert to matrix
                    SourceOfFragmentationPatterns = SourceOfFragmentationPatterns.astype(str) #save as class object with type string
                elif dataFrame.iloc[rowIndex][0] == 'sourceOfIonizationData':
                    dfsourceOfIonizationData = dataFrame.iloc[rowIndex][1::2] #Select the row of names
                    sourceOfIonizationData = dfsourceOfIonizationData.values #convert to matrix
                    sourceOfIonizationData = sourceOfIonizationData.astype(str) #save as class object with type string
                elif dataFrame.iloc[rowIndex][0] == 'Molecules': #if the abscissa titles the molecule names
                    dfmolecules = dataFrame.iloc[rowIndex][1::2] #select the row of names
                    molecules = dfmolecules.values #convert to matrix
                    molecules = molecules.astype(str) #save as class object with type string
                    molecules = list(molecules)
                    for moleculeIndex in range(len(molecules)):
                        molecules[moleculeIndex] = molecules[moleculeIndex].strip()#remove leading and trailing whitespaces.
                elif dataFrame.iloc[rowIndex][0] == 'Electron Numbers': #if the abscissa titles the electron numbers
                    dfelectronnumbers = dataFrame.iloc[rowIndex][1::2] #select the row of names
                    electronnumbers = dfelectronnumbers.values #convert to matrix
                    electronnumbers = electronnumbers.astype(int) #save as class object with type int
                elif dataFrame.iloc[rowIndex][0] == 'Molecular Mass': #if the abscissa titles the molecular weights
                    dfmolecularWeights = dataFrame.iloc[rowIndex][1::2] #select row of names
                    molecularWeights = dfmolecularWeights.values #convert to matrix
                    molecularWeights = molecularWeights.astype(numpy.double) #save as class object with type double
                elif dataFrame.iloc[rowIndex][0] == 'moleculeIonizationType':
                    dfmoleculeIonizationType = dataFrame.iloc[rowIndex][1::2] #select row of names
                    moleculeIonizationType = dfmoleculeIonizationType.values #convert to matrix
                    moleculeIonizationType = moleculeIonizationType.astype(str) #save as class object with type string
                elif (dataFrame.iloc[rowIndex][0] == 'relativeIonizationEfficiencies') or (dataFrame.iloc[rowIndex][0] == 'knownIonizationFactorsRelativeToN2'):
                    dfrelativeIonizationEfficiencies = dataFrame.iloc[rowIndex][1::2] #select row of names
                    relativeIonizationEfficiencies = dfrelativeIonizationEfficiencies.values #convert to matrix
                    for index in range(len(relativeIonizationEfficiencies)):
                        try: #try to convert to a float
                            relativeIonizationEfficiencies[index] = float(relativeIonizationEfficiencies[index])
                        except: #if not possible, the value is probably None or 'unknown' so leave as a string
                            pass
#                    relativeIonizationEfficiencies = relativeIonizationEfficiencies.astype(float) #save as class object with type float

        try: #if using an older reference file, it will not have ionization factors so the elif statement never gets entered meaning knownIonizationFactors does not exist
            relativeIonizationEfficiencies #Try calling this variable, if it exists there will be no error
        except: #if it does not exist, populate it with unknown
            relativeIonizationEfficiencies = ['unknown'] #initialize as a list of len(1)
            relativeIonizationEfficiencies = relativeIonizationEfficiencies*len(molecules)   #duplicate the list for # of molecules.
            relativeIonizationEfficiencies = numpy.array(relativeIonizationEfficiencies) #convert to matrix
            
        try: #if using an old reference file, it will not have ionization types so the elif statement never gets entered meaning moleculeIonizationType does not exist
            moleculeIonizationType #Try calling this variable, if it exists there will be no error
        except: #if it does not exist, populate it with unknown
            moleculeIonizationType = ['unknown'] #initialize as a list of len(1)
            moleculeIonizationType = moleculeIonizationType *len(molecules)  #duplicate the list for # of molecules.
            moleculeIonizationType = numpy.array(moleculeIonizationType) #convert to matrix
            
        try: #If using an older reference file, it will not have SourceOfIonizationInfo so the elif statement never gets entered meaning the variable does not exist
            sourceOfIonizationData #try calling the variable, if it exists there will be no error
        except: #If it does not exist, populate with empty strings
            sourceOfIonizationData = ['unknown'] #initialize as a list of len(1)
            sourceOfIonizationData = sourceOfIonizationData * len(molecules) #duplicate the list for # of molecules.
            sourceOfIonizationData = numpy.array(sourceOfIonizationData) #convert to matrix
          
                
#        ''' generate reference matrix'''
#        #remove top 4 rows
#        dfreference = dataFrame.iloc[4:][:]
#        #convert to matrix
#        reference = dfreference.values
#        #convert the matrix to floats
#        provided_reference_patterns = reference.astype(float)
#        #convert reference from XYXY to XYYY
#        provided_reference_patterns=FromXYXYtoXYYY(provided_reference_patterns)
#        #clear rows of zeros
#        provided_reference_patterns=DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(provided_reference_patterns,startingRowIndex=1)
#
#        '''generate electron number list'''
#        #select row of electron numbers
#        dfelectronnumbers = dataFrame.iloc[2][1::2]
#        #convert to matrix
#        electronnumbers = dfelectronnumbers.values
#        #save as class object with type int
#        electronnumbers = electronnumbers.astype(int)
#   
#        '''generate list of molecule names'''
#        #select row of names
#        dfmolecules = dataFrame.iloc[1][1::2]
#        #convert to matrix
#        molecules = dfmolecules.values
#        #save as class object with type string
#        molecules = molecules.astype(str)
#        
#        '''generate list of molecular weights'''
#        #select row of names
#        dfmolecularWeights = dataFrame.iloc[3][1::2]
#        #convert to matrix
#        molecularWeights = dfmolecularWeights.values
#        #save as class object with type float
#        molecularWeights = molecularWeights.astype(float)
#        
#        '''generate list of source information'''
#        #select row of names
#        dfsourceInfo = dataFrame.iloc[0][1::2]
#        #convert to matrix
#        sourceInfo = dfsourceInfo.values
#        #save as class object with type string
#        sourceInfo = sourceInfo.astype(str)
        
        '''list of massfragments monitored is not part of reference file'''
        mass_fragment_numbers_monitored = None
    else:
        print("Form of reference pattern not recognized. Only XYYY or XYXY allowed."); sys.exit()
    return provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form

'''
getMoleculesFromReferenceData is a function that takes in the reference filename and just returns the molecules present
This function is written for XYYY data and does not work with XYXY data. Use of this function is discouraged. The readReferenceFile function should be used to pull the molecule names instead. 
'''
#def getMoleculesFromReferenceData(ReferenceFileName):
#    #Read the csv file
#    #TODO Change to use numpy.gen_from_text instead of pandas
#    #TODO: add an optional argument for the form "xxxy" vs "xyxy".
#    ReferenceInfo = pandas.read_csv(ReferenceFileName,header=0)
#    #Convert the reference info into an array
#    ReferenceInfoArray = numpy.array(ReferenceInfo)
#    #Get the names of molecules in reference data
#    molecules = ReferenceInfoArray[0,1:] #First row, all but the first column
#    #Strip the white space around the molecules
#    for i in range(0,len(molecules)):
#        molecules[i] = molecules[i].strip()
#    return molecules

'''
getMassFragmentsFromCollectedData is a function that takes in the collected filename and returns the mass fragments present in the data
'''
def getMassFragmentsFromCollectedData(CollectedFileName):
    #Read the csv file
    #TODO CHange to use numpy.gen_from_text instead of pandas
    DataInfo = pandas.read_csv(CollectedFileName,header=0)
    #Convert the data into an array
    DataInfoArray = numpy.array(DataInfo)
    #Get the names of mass fragments in collected data
    massFragments = DataInfoArray[0,1:] #First row, all but the first column
    #Remove the 'm' from each fragment and convert to float (i.e. 'm28' now becomes 28.)
    for i in range(0,len(massFragments)):
        massFragments[i] == str(massFragments[i])
        massFragments[i] = massFragments[i].replace("mass",'')
        massFragments[i] = massFragments[i].replace("m",'')
        massFragments[i] = massFragments[i].strip()
        massFragments[i] = float(massFragments[i])
    return massFragments

'''
getIE_Data reads in the ionization data file name and generates a dictionary filled with MolecularIonizationData objects
'''
def getIE_Data(IonizationDataFileName):
    ionizationData = numpy.genfromtxt(IonizationDataFileName,dtype=None,delimiter=',',encoding=None)
    AllMID_ObjectsDict = {} #initialize the MID Dictionary
    
    #For loop to get the column indexes of information in the ionization data csv file
    for colIndex in range(len(ionizationData[0])): #Loop across each column in the ionizationData
        if ionizationData[0][colIndex] == 'Name': #If header is name, then this column index points to molecule names
            NameIndex = colIndex
        elif ionizationData[0][colIndex] == 'RS': #If header is RS, then this column index points to the sensitivity factors
            RSIndex = colIndex
        elif ionizationData[0][colIndex] == 'Electron Number': #If header is Electron Number, then this column index points to molecules' electron numbers
            ENumberIndex = colIndex
        elif ionizationData[0][colIndex] == 'Type': #If header is Type, then this column index points to the molecules' ionization types
            TypeIndex = colIndex
        elif ionizationData[0][colIndex] == 'sourceOfIonizationData': #If header is source, then this column index points to the source of information
            SourceIndex = colIndex

    for rowIndex in range(1,len(ionizationData)): #loop across each row (not including the header)
        moleculeName = ionizationData[rowIndex][NameIndex] #get the molecule name
        RS_Value = ionizationData[rowIndex][RSIndex] #Get the ionization factor
        moleculeElectronNumber = ionizationData[rowIndex][ENumberIndex] #get the molecule's electron number
        moleculeIonizationType = ionizationData[rowIndex][TypeIndex] #get the molecule 
        moleculeIonizationTypeList = moleculeIonizationType.split(';')
        sourceOfIonizationData = ionizationData[rowIndex][SourceIndex]
        
        MID_ObjectName = moleculeName + '_IE' #The object name will be the moleculeName_IE
        
        if MID_ObjectName in AllMID_ObjectsDict: #If we already have this molecule in AllMID_ObjectsDict then we just need to add the RS_Value and source to the appropriate list
            AllMID_ObjectsDict[MID_ObjectName].addData(RS_Value,sourceOfIonizationData)
        else: #otherwise the object does not exist and needs to be created
            AllMID_ObjectsDict[MID_ObjectName] = MolecularIonizationData(moleculeName,RS_Value,moleculeElectronNumber,moleculeIonizationTypeList,sourceOfIonizationData) #Store MIDObject in AllMIDObjects Dictionary
    return AllMID_ObjectsDict

'''
populateAllMID_ObjectsDict tries to read the ionizationFileName and populates the dictionary
If unable to read the file an empty dictionary is returned
'''
def populateAllMID_ObjectsDict(ionizationDataFileName):
    try:
        AllMID_ObjectsDict = getIE_Data(ionizationDataFileName) #Read the ionization data and put the information into a dictionary
    except: #If the ionization file does not exist in the main directory, leave as an empty dictionary
        print('Warning: ionizationDataFileName was not found. Simple heuristics will be used for ionization efficiencies.')
        AllMID_ObjectsDict = {}
        
    return AllMID_ObjectsDict

###############################################################################
#########################  Classes: Data Storage  #############################
###############################################################################
#TODO add warning to user if their data contains any NaN values
# having NaN values may crash the program 
class MSData (object):
    #creates an MSData object that has the following sub-bojects:
		#self.times, 1D
		#self.mass_fragment_numbers , 1D and must be integers
		#self.rawCollectedData, a 2D array of the signals.
		
    def __init__(self, mass_fragment_numbers, abscissaHeader, times, rawCollectedData, collectedFileName=None):
        
        self.mass_fragment_numbers, self.abscissaHeader, self.times, self.rawCollectedData, self.collectedFileName=mass_fragment_numbers, abscissaHeader, times, rawCollectedData, collectedFileName
        #class object variable created to allow class to be used separately from the program. 
        self.ExportAtEachStep = ''
        
        '''create data set to work on'''
        self.workingData = self.rawCollectedData
        
        '''initilize variables to simplify future code'''
        self.datafromcsv = []
        #start the timer function
        self.previousTime = timeit.default_timer()
        #initalize debugging lists
        #These lists are appended in parallel so a variable of the same index from each list will be related
        self.runTimeAtExport = []
        self.labelToExport = []
        self.dataToExport = []
        self.experimentTimes = []
        
    def ExportCollector(self, callingFunction):
        #record current time
        currentTime = timeit.default_timer()
        #add net time to list of run times
        self.runTimeAtExport.append(currentTime - self.previousTime)
        #record current time for next function's use
        self.previousTime = currentTime
        #add the name of the calling function to mark its use
        self.labelToExport.append(callingFunction) 
        
        if self.ExportAtEachStep == 'yes':
            #record data of experiment
            self.dataToExport.append(self.workingData.copy())
            #record times from the data of the experiment
            self.experimentTimes.append(self.times.copy())
            
    def ExportMSData(self):
        print("\n Collected Export List:")
        for savePoint in range(len(self.runTimeAtExport)):
            print(self.labelToExport[savePoint])
            print(self.runTimeAtExport[savePoint])
            if self.ExportAtEachStep == 'yes':
                #inserting the data for a particular savePoint
                filename = 'Exported%s%s.csv'%(savePoint, self.labelToExport[savePoint]) #FIXME: for DataSmoother, and some others, the debug output has a "Time" header but the time is not exported.
                data = self.dataToExport[savePoint]
                abscissa = self.experimentTimes[savePoint]
                colIndex = ['m%s'% int(y) for y in self.mass_fragment_numbers]
                DataFunctions.MSDataWriterXYYY(filename, data, abscissa, colIndex, self.abscissaHeader)
                                        
class MSReference (object):
    def __init__(self, provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns=[], sourceOfIonizationData=[], relativeIonizationEfficiencies=[], moleculeIonizationType=[], mass_fragment_numbers_monitored=None, referenceFileName=None, form=None, AllMID_ObjectsDict={}):
        self.provided_reference_patterns, self.electronnumbers, self.molecules, self.molecularWeights, self.SourceOfFragmentationPatterns, self.sourceOfIonizationData, self.relativeIonizationEfficiencies, self.moleculeIonizationType, self.mass_fragment_numbers_monitored, self.referenceFileName, self.form = provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form
        #class object variable created to allow class to be used separately from the program. 
        self.ExportAtEachStep = ''
        self.iterationSuffix = ''
        if type(referenceFileName) != type(None):
            self.referenceFileNameExtension = self.referenceFileName.split(".")[1]
        #This loops through the molecules, and removes whitespaces from before and after the molecule's names.
        for moleculeIndex, moleculeName in enumerate(self.molecules):
            self.molecules[moleculeIndex] = moleculeName.strip()     
        #Initializing some list variables with defaults.
        if len(self.SourceOfFragmentationPatterns) == 0:
            self.SourceOfFragmentationPatterns = ['unknown']* len(self.molecules)
        if len(self.sourceOfIonizationData) == 0:
            self.sourceOfIonizationData = ['unknown']* len(self.molecules)       
            self.sourceOfIonizationData = numpy.array(self.sourceOfIonizationData, dtype="str")
        else:
            self.sourceOfIonizationData = numpy.array(self.sourceOfIonizationData, dtype="str")
        if len(self.relativeIonizationEfficiencies) == 0:
            self.relativeIonizationEfficiencies = ['unknown']* len(self.molecules)
            self.relativeIonizationEfficiencies = numpy.array(self.relativeIonizationEfficiencies, dtype="object") #can have Strings or floats, so use dtype="object"
        else:
            self.relativeIonizationEfficiencies = numpy.array(self.relativeIonizationEfficiencies, dtype="object") #can have Strings or floats, so use dtype="object"
        if len(self.moleculeIonizationType) == 0:
            self.moleculeIonizationType = ['unknown']* len(self.molecules)
        self.provided_mass_fragments = self.provided_reference_patterns[:,0]
        #clear ClearZeroRowsFromProvidedReferenceIntensities
        #self.ClearZeroRowsFromProvidedReferenceIntensities()  #commenting this out to avoid problems for now.
        #initialize the standardized_reference_patterns
        self.standardized_reference_patterns=StandardizeReferencePattern(self.provided_reference_patterns,len(self.molecules))
            
        #Initializing Export Collector Variables
        #start the timer function
        self.previousTime = timeit.default_timer()
        #initalize debugging lists
        #These lists are appended in parallel so a variable of the same index from each list will be related
        self.runTimeAtExport = []
        self.labelToExport = []
        self.dataToExport = []
        self.moleculesToExport = []
        self.exportSuffix = ''
        #self.experimentTimes = []       
        #Get ionization efficiencies and export their values and what method was used to obtain them
        self.AllMID_ObjectsDict = AllMID_ObjectsDict
        self.populateIonizationEfficiencies(self.AllMID_ObjectsDict)
        self.exportIonizationInfo()

    #a small helper function to check if an extension exists in a filename and to return the delimiter based on that.
    def getDelimiterFromExtension(self, filename=''):
        if filename == '':
            filename = self.referenceFileName
        if ".tsv" in filename:
            delimiter = '\t'
        elif ".tab" in filename:
            delimiter = '\t'
        elif ".txt" in filename:
            delimiter = '\t'        
        elif ".skv" in filename:
            delimiter = ';'
        elif ".csv" in filename:
            delimiter = ',' #it could be something else, but we will assume that a csv
        else:
            delimiter = '\t' #for MSRESOLVE, this is now the default delimiter.
        return delimiter

    #This function allows adding molecules to an existing reference patterns. When using TuningCorrector it is used to create MixedReference patterns.
    #Though these variable names are plural, they are expected to be lists of one. "molecules" is supposed to be a list of variable names.
    #provided_reference_patterns should be in an XYYY format.  If starting with XYXY data, is okay to feed a single "XY" at a time and to do so repeatedly in a loop.
    def addMolecules(self, provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns=[], sourceOfIonizationData=[], relativeIonizationEfficiencies=[], moleculeIonizationType=[], mass_fragment_numbers_monitored=None, referenceFileName=None, form=None, AllMID_ObjectsDict={}):      
        #In case somebody decides to try to feed single molecule's name directly.
        if type(molecules) == type("string"):
            molecules = [molecules]
        #Initializing some list variables with defaults.
        if len(SourceOfFragmentationPatterns) == 0:
           SourceOfFragmentationPatterns = ['unknown']* len(self.molecules)
        if len(sourceOfIonizationData) == 0:
           sourceOfIonizationData = ['unknown']* len(self.molecules)       
        if len(relativeIonizationEfficiencies) == 0:
           relativeIonizationEfficiencies = ['unknown']* len(self.molecules)
        if len(moleculeIonizationType) == 0:
           moleculeIonizationType = ['unknown']* len(self.molecules)
           
        #Now to concatenate to the existing reference objects. Can't use ".extend" because that modifies the original list while returning "None."
        self.electronnumbers = list(self.electronnumbers) + (list(electronnumbers))
        self.molecules = list(self.molecules) + (list(molecules))
        self.molecularWeights = list(self.molecularWeights) + (list(molecularWeights))
        self.SourceOfFragmentationPatterns = list(self.SourceOfFragmentationPatterns) + (list(SourceOfFragmentationPatterns))
        self.sourceOfIonizationData = list(self.sourceOfIonizationData) + (list(sourceOfIonizationData))
        self.relativeIonizationEfficiencies = list(self.relativeIonizationEfficiencies) + (list(relativeIonizationEfficiencies))
        self.moleculeIonizationType = list(self.moleculeIonizationType) + (list(moleculeIonizationType))
        
        #Will to standardize the molecules reference patterns before adding. This is important because the addXYYYtoXYYY function drops values that are small like 1E-8.
        standardized_reference_patterns=StandardizeReferencePattern(provided_reference_patterns,len(molecules))
        #Will add to both self.provided_reference_patterns  and self.standardized_reference_patterns
        self.provided_reference_patterns = DataFunctions.addXYYYtoXYYY(self.provided_reference_patterns, standardized_reference_patterns)
        self.standardized_reference_patterns = DataFunctions.addXYYYtoXYYY(self.standardized_reference_patterns, standardized_reference_patterns)
        #Rather than adding to the standarized reference patterns separately, will just go ahead and restandardize.
        self.standardized_reference_patterns=StandardizeReferencePattern(self.standardized_reference_patterns,len(self.molecules))
        return self

    #For removing molecules we will use an external function trimDataMoleculesToMatchChosenMolecules
    #Note that this actually makes a new reference object!
    def removeMolecules(self, listOfMoleculesToRemove):
        if type(listOfMoleculesToRemove) == type("string"):
            listOfMoleculesToRemove = [listOfMoleculesToRemove]

        moleculesToKeep = list(copy.deepcopy(self.molecules)) #make a copy of molecules.
        #now go through the original list of molecules and pop unwanted molecules from moleculesToKeep
        for moleculeName in self.molecules:
            if moleculeName in listOfMoleculesToRemove:
                moleculesToKeep.remove(moleculeName)
        trimmedRefererenceData = trimDataMoleculesToMatchChosenMolecules(self, moleculesToKeep)        
        return trimmedRefererenceData
        
    #TODO exportCollector should be updated to take in a string argument for the data type that it should record (patterns vs various intensities)
    #Additionally, it should take an optional variable to determine the headers that will be used.         
    #Basically, the logic in here is pretty bad!
    def ExportCollector(self, callingFunction, use_provided_reference_patterns = False, export_matching_correction_value = False, export_uncertainties = False, export_standard_uncertainties=False, export_tuning_uncertainties=False, export_relative_uncertainties=False):
        #record current time
        currentTime = timeit.default_timer()
        #Don't do anything unless ExportAtEachStep is on.
        if self.ExportAtEachStep == 'yes':
            #add net time to list of run times
            self.runTimeAtExport.append(currentTime - self.previousTime)
            #record current time for next function's use
            self.previousTime = currentTime
            #add the name of the calling function to mark its use
            self.labelToExport.append(callingFunction) 
            ##record molecules of experiment
            self.moleculesToExport.append(self.molecules.copy())
            #record data 
            if export_uncertainties:
                self.dataToExport.append(self.provided_reference_patterns_absolute_uncertainties.copy())
            elif export_matching_correction_value:
                self.dataToExport.append(self.reciprocal_matching_correction_values.copy())     #Export Correction values along with provided reference patterns
            elif export_tuning_uncertainties:
                self.dataToExport.append(self.standardized_reference_patterns_tuning_uncertainties.copy())
            elif export_standard_uncertainties:
                self.dataToExport.append(self.absolute_standard_uncertainties.copy())
            elif export_relative_uncertainties:
                self.dataToExport.append(self.relative_standard_uncertainties.copy())               
            elif use_provided_reference_patterns:
                self.dataToExport.append(self.provided_reference_patterns.copy())
            elif callingFunction == 'UnnecessaryMoleculesDeleter':
                self.dataToExport.append(self.monitored_reference_intensities.copy())
            elif not use_provided_reference_patterns:
                self.dataToExport.append(self.standardized_reference_patterns.copy())
            
    #this function exports only an XYYY file. To export a full referenceData file, use the function exportReferencePattern
    def ExportFragmentationPatterns(self, verbose=True):
        #Only print if not called from interpolating reference objects
        if verbose:
            print("\n Reference Debugging List:")
        for savePoint in range(len(self.runTimeAtExport)):
            #Only print if not called from interpolating reference objects
            if verbose:
                print(self.labelToExport[savePoint])
                print(self.runTimeAtExport[savePoint])
            if self.ExportAtEachStep == 'yes':
                #inserting the data for a particular savePoint
                delimiter = getDelimiterFromExtension(self.referenceFileNameExtension)
                filename = 'Exported%s%s.%s'%(savePoint, self.labelToExport[savePoint], self.referenceFileNameExtension)
                data = self.dataToExport[savePoint]
                colIndex = ['%s'% y for y in self.moleculesToExport[savePoint]]
                #colIndex = ['%s'% y for y in self.molecules]
                #print(numpy.shape(data),numpy.shape(colIndex))
                ExportXYYYData(filename,data,colIndex, fileSuffix = self.iterationSuffix)

    # This class function removes all rows of zeros from
    # the XYYY sorted provided_reference_patterns, and *also* provided_mass_fragments
    #The logic in the below funtion is badly written, in terms of efficiency. But it seems to work at present.
    #TODO: This is not a good practice, because provided_reference_patterns is getting changed, no longer "Provided".
    #TODO: (continued from previous line) It's more like "zero_trimmed" reference intensities after this.
    #TODO: it's better to clear zero rows from the standardized_reference_patterns
    def ClearZeroRowsFromProvidedReferenceIntensities(self):
        #initial a counter for the row index, which will be updated during the loop
        currentRowIndexAccountingForDeletions = 0
        #provided_reference_patternsOnly is not used, but is made for future use (see below)
        provided_reference_patternsOnly = self.provided_reference_patterns[:,1:]
        for intensitiesOnlyInRow in provided_reference_patternsOnly:
            #This line checks if there are any non-zeros in the row.
            numberOfNonzeros = numpy.count_nonzero(intensitiesOnlyInRow)
            if numberOfNonzeros == 0 :
                #If there are only zeros. we delete a row and adjust the row index to account for that deletion.
                self.provided_reference_patterns = numpy.delete(self.provided_reference_patterns, currentRowIndexAccountingForDeletions, axis=0 ) #axis = 0 specifies to delete rows (i.e. entire abscissa values at the integer of currentRowIndexAccountingForDeletions).
                self.provided_mass_fragments = numpy.delete(self.provided_mass_fragments, currentRowIndexAccountingForDeletions, axis=0 )
                try: #Under normal situations, this try is a bit like an implied if G.calculateuncertaintiesInConcentrations != None:
                    self.absolute_standard_uncertainties = numpy.delete(self.absolute_standard_uncertainties, currentRowIndexAccountingForDeletions, axis=0 )
                except:
                    pass                                                                                                                                     
                currentRowIndexAccountingForDeletions = currentRowIndexAccountingForDeletions -1
            #whether we deleted rows or not, we increase the counter of the rows.
            currentRowIndexAccountingForDeletions = currentRowIndexAccountingForDeletions + 1

    def ClearZeroRowsFromStandardizedReferenceIntensities(self):
        #initial a counter for the row index, which will be updated during the loop
        currentRowIndexAccountingForDeletions = 0
        #standardized_reference_patterns_unedited is not used, but is made for future use (see below)
        standardized_reference_patterns_unedited = self.standardized_reference_patterns[:,1:]
        self.standardized_reference_patterns_mass_fragments = self.standardized_reference_patterns[:,0]
        for intensitiesOnlyInRow in standardized_reference_patterns_unedited:
            #This line checks if there are any non-zeros in the row.
            numberOfNonzeros = numpy.count_nonzero(intensitiesOnlyInRow)
            if numberOfNonzeros == 0 :
                #If there are only zeros. we delete a row and adjust the row index to account for that deletion.
                self.standardized_reference_patterns = numpy.delete(self.standardized_reference_patterns, currentRowIndexAccountingForDeletions, axis=0 ) #axis = 0 specifies to delete rows (i.e. entire abscissa values at the integer of currentRowIndexAccountingForDeletions).
                self.standardized_reference_patterns_mass_fragments = numpy.delete(self.standardized_reference_patterns_mass_fragments, currentRowIndexAccountingForDeletions, axis=0 )
                try: #Under normal situations, this try is a bit like an implied if G.calculateuncertaintiesInConcentrations != None:
                    self.absolute_standard_uncertainties = numpy.delete(self.absolute_standard_uncertainties, currentRowIndexAccountingForDeletions, axis=0 )
                except:
                    print("Debugging line 2892") #This hardcoded print statement should not be removed.
                    pass                                                                                                                                     
                currentRowIndexAccountingForDeletions = currentRowIndexAccountingForDeletions -1
            #whether we deleted rows or not, we increase the counter of the rows.
            currentRowIndexAccountingForDeletions = currentRowIndexAccountingForDeletions + 1
        if hasattr(self, 'absolute_standard_uncertainties'):
            self.ExportCollector("ClearZeroRowsFromAbsoluteStandardUncertainties", export_standard_uncertainties=True)
            self.update_relative_standard_uncertainties()
            self.ExportCollector("ClearZeroRowsFromRelativeStandardUncertainties", export_relative_uncertainties=True)
        else:
            print("Warning: line 2897 was unable to update the relative uncertainties of a reference pattern.")
        self.ExportCollector("ClearZeroRowsFromStandardizedReferenceIntensities", use_provided_reference_patterns=False)
        
    def extendMassFragments(self, massFragmentsToExtendBy):
        #This will add in rows of zeros to the standardized_reference_patterns and also the uncertainties arrays.
        self.standardized_reference_patterns_mass_fragments = self.standardized_reference_patterns[:,0]*1.0
        self.standardized_reference_patterns = DataFunctions.extendXYYYtoZYYY(self.standardized_reference_patterns, massFragmentsToExtendBy)
        #Check if the absolute_standard_uncertainties is the right length. If it is, we will be extending it also.
        if hasattr(self, "absolute_standard_uncertainties"):
            if len(self.absolute_standard_uncertainties[:,0] ==  self.standardized_reference_patterns_mass_fragments):
                self.absolute_standard_uncertainties = DataFunctions.extendXYYYtoZYYY(self.absolute_standard_uncertainties, massFragmentsToExtendBy)
                #For the relative uncertainties, we will use the built in function.
                self.update_relative_standard_uncertainties()
            
    #This class function converts the XYXY data to an XYYY format
    def FromXYXYtoXYYY(self):
        masslists = [] #future lists must be must empty here to append in the for loops
        relativeintensitieslists = [] #future list
        #this loops gathers all the mass fragment numbers for each molecule in one list of arrays, while a second
        #list is made, gathering the relative intensities so that they were indexed the same as their mass fragment
        #numbers in the other list
        #this for loop grabs stuff from the reference array, whose orientation and identity is shown in the flow chart arrays document
        for referenceBy2Index in range(0,len(self.provided_reference_patterns[0,:]),2):#array-indexed for loop, only gets every other value, as half the indexes are mass lists, and the other half are relative intensity
            masslists.append(self.provided_reference_patterns[:,referenceBy2Index])#these are lists of arrays
            relativeintensitieslists.append(self.provided_reference_patterns[:,referenceBy2Index+1])#the relative intensities are after every counter, so there is a +1 (it is array indexed so since the first column is a mass list all the +1's are relative intensities)
        masslist = [] #future list
        #This for loop gets all of the mass fragments from the first index of the list, basically by not adding the 
        #'nan's or empty spaces after the numbers
        for referenceIndex in range(len(self.provided_reference_patterns[:,0])): #array-indexed for loop
            if str(masslists[0][referenceIndex]) != 'nan': #we do not want nan's in our array, the genfromtxt function calls empty boxes in excel (might be in .csv as well)'nan'.
                masslist.append(masslists[0][referenceIndex])
        #this second nested for loop gathers all the other mass fragment numbers that have not already been added to
        #the masslist, basically obtaining all the masses in the reference data and then after the loop they are sorted
        #using .sort, then an empty array of zeros is made to accommodate the output array
        for masslistIndex in range(1,len(masslists)):#array-indexed for loop, starts at one because it's checking against all arrays besides itself
            for referenceIndex in range(len(self.provided_reference_patterns[:,0])):#array-indexed for loop
                if str(masslists[masslistIndex][referenceIndex]) != 'nan':
                    if sum(masslists[masslistIndex][referenceIndex] == numpy.array(masslist)) == 0:#if the value being looked at is not equal to anything in our masslist already
                        masslist.append(masslists[masslistIndex][referenceIndex])
        masslist.sort()#puts the list in order
        reference_holder = numpy.zeros([len(masslist),len(self.provided_reference_patterns[0,:])/2+1])#makes an array that is full of zeros to hold our future reference array
        reference_holder[:,0:1] = numpy.vstack(numpy.array(masslist))#This puts the mass list in the first column of our new reference array
        #Finally, the for loop below makes a list each revolution, comparing each list of mass fragments (for each molecule)
        #and adding the relative intensities (from the identically indexed array) when they numbers were equal, and otherwise
        #adding a zero in its place. It then adds this list to the array (using numpy.vstack and numpy.array)
        for massListsIndex in range(len(masslists)):#array-indexed for loop
            relativeintensitieslist = [] #empties the list every loop
            for massListIndex in range(len(masslist)):
                placeholder = 0 #after the next for loop finishes, this is reset
                for specificMassListIndex in range(len(masslists[massListsIndex])):#array-indexed for loop, each column of .csv file being checked
                    if masslists[massListsIndex][specificMassListIndex] == masslist[massListIndex]:#This is when the index for the correct mass fragment is found
                        relativeintensitieslist.append(relativeintensitieslists[massListsIndex][specificMassListIndex])#relative intensities lists index the same way
                        placeholder = 1 #so that the next if statement will know that this has happened
                    if specificMassListIndex == len(masslists[massListsIndex])-1 and placeholder == 0:#If it comes to the end of the for loop, and there's no match, then the relative intensity is zero
                        relativeintensitieslist.append(0)
                if massListIndex == len(masslist)-1:#Once the larger for loop is done the 
                    reference_holder[:,(massListsIndex+1):(massListsIndex+2)] = numpy.vstack(numpy.array(relativeintensitieslist)) #the list is made into an array and then stacked (transposed)
        self.provided_reference_patterns = reference_holder

    def addSuffixToSourceOfFragmentationPatterns(self, suffixString):
        #this modifies the original object so that it is not necessary to return anything.
        self.SourceOfFragmentationPatterns = list(self.SourceOfFragmentationPatterns)
        for moleculeIndex, sourceString in enumerate(self.SourceOfFragmentationPatterns):
            self.SourceOfFragmentationPatterns[moleculeIndex] = sourceString + suffixString
        
    #populateIonizationEfficiencies is an MSReference function that populates a variable, relativeIonizationEfficiencies, that contains the ionization factors used in CorrectionValuesObtain
    #If the ionization factor is known and in the reference data, then that value is used
    #If the ionization factor is unknown the the function will look in the MID Dictionary and check if the molecule exists in the ionization data.  If it does the the ionization average of the ionization factors for that particular molecule in the data is used
    #If the ionization factor is unknown and the particular molecule does not exist in the MID Data, then the function checks the molecule's ionization type(s).  The function will take all molecules from the MID data that have the same type and will perform a linear fit on the data.  The ionization factor for this molecule is determined based on the linear fit and number of electrons
    #If the ionization factor is unknown, the molecule does not exist in the MID data, and the molecule's ionization type is unknown, then the function defaults to the Madix and Ko equation
    def populateIonizationEfficiencies(self, AllMID_ObjectsDict={}):
        #self.relativeIonizationEfficiencies = numpy.zeros(len(self.molecules)) #initialize an array the same length as the number of molecules that will be populated here and used in CorrectionValuesObtain
        #self.sourceOfIonizationData = copy.copy(self.molecules) #initialize an array to store which method was used to obtain a molecule's ionization factor
        self.sourceOfIonizationData = list(self.sourceOfIonizationData) #Making it a list here because if it's a numpy array the strings may get cut short.
        for moleculeIndex in range(len(self.molecules)): #loop through our initialized array
            if isinstance(self.relativeIonizationEfficiencies[moleculeIndex],float): #if the knownIonizationFactor is a float, then that is the value defined by the user
                self.relativeIonizationEfficiencies[moleculeIndex] = self.relativeIonizationEfficiencies[moleculeIndex]
                self.sourceOfIonizationData[moleculeIndex] = self.sourceOfIonizationData[moleculeIndex] #the molecule's factor was known
            else: #Ionization factor is not known so look at molecular ionization data from literatiure 
                #Initialize three lists
                MatchingMID_Objects = []
                MatchingMID_RS_Values = []
                MatchingMID_ElectronNumbers = []
                #Initialize a flag to overwrite if a molecule is in both self.molecules and the MID_ObjectDict
                matchingMolecule = False
                if self.relativeIonizationEfficiencies[moleculeIndex] == None or self.relativeIonizationEfficiencies[moleculeIndex] == 'unknown': #the ionization factor is not known so look in the AllMID_ObjectsDict
                    currentMoleculeKeyName = self.molecules[moleculeIndex] + '_IE' #holder variable to store the current molecule name + '_IE' to match the keys of AllMID_ObjectsDict (e.g. Acetaldehyde_IE)
                    if currentMoleculeKeyName in AllMID_ObjectsDict.keys(): #If the current molecule is in the dictionary, use its RS_Value
                        MatchingMID_Objects.append(currentMoleculeKeyName) #append the key
                        MatchingMID_RS_Values.append(numpy.mean(AllMID_ObjectsDict[currentMoleculeKeyName].RS_ValuesList)) #append the average RS_Value
                        MatchingMID_ElectronNumbers.append(AllMID_ObjectsDict[currentMoleculeKeyName].electronNumber) #append the electron number
                        matchingMolecule = True #set the flag to be true
                if matchingMolecule == True: #If the molecule matches a molecule in the MID dictionary, use the average RS_Value
                    self.relativeIonizationEfficiencies[moleculeIndex] = MatchingMID_RS_Values[0]
                    self.sourceOfIonizationData[moleculeIndex] = 'knownIonizationFactorFromProvidedCSV' #A molecule in the reference data is also in the ionization data
                elif matchingMolecule == False: #Otherwise matchingMolecule is False which means its not in the data from literature.  So we will approximate the ionization factor based on a linear fit of the data from literature that share the molecule's type or use the Madix and Ko equation
                    if type(self.moleculeIonizationType[moleculeIndex]) != type(None) and self.moleculeIonizationType[moleculeIndex] != 'unknown': #IF the user did not manually input the ionization factor and none of the molecules in the MID_Dict matched the current molecule
                        #Then get an estimate by performing a linear fit on the data in the MID Dictionary
			#TODO:The program currently only takes in one type but it is a desired feature to allow users to put in multiple types such as type1+type2 which would make a linear fit of the combined data between the two types
			#TODO continued:The user should also be able to put in type1;type2 and the program would find the ionization factor using a linear fit of data from type1 and using a linear fit of data from type2.  The largest of the two ionization factors would be used.
			#TODO continued:Then doing type1+type2;type3 would take the larger value between the linear fit of the combined type1 and type2 data or the value from the linear fit of type3 data
                        for key in AllMID_ObjectsDict: #Loop through the MID Dictionary
                            for MID_MoleculeType in AllMID_ObjectsDict[key].moleculeIonizationType: #Loop through the ionization types to get all the ionization types of a particular molecule (e.g. Ethanol is both an alcohol and a hydrogen non-metal-ide so its RS value(s) will be included if the user has a molecule that is either an alcohol or a hydrogen non-metal-ide)
                                #Use stringCompare to check if a molecule in the MID Dictionary matches a molecule in the reference data since casing and spacing may differ between the two (e.g. reference data may have carbon dioxide while MID Dictionary may have Carbon Dioxide)
                                if parse.stringCompare(self.moleculeIonizationType[moleculeIndex],MID_MoleculeType): #If the knownMoleculeType matches an MID object's molecule type
                                    MatchingMID_Objects.append(key) #Append the key
                                    MatchingMID_RS_Values.append(numpy.mean(AllMID_ObjectsDict[key].RS_ValuesList)) #Append the average of the RS values
                                    MatchingMID_ElectronNumbers.append(AllMID_ObjectsDict[key].electronNumber) #append the electron number
                        if len(MatchingMID_Objects) == 1: #If only one data point add (0,0) to the data
                            MatchingMID_Objects.append(MatchingMID_Objects[0]) #append the molecule type to the objects list
                            MatchingMID_RS_Values.append(0.0)
                            MatchingMID_ElectronNumbers.append(0.0)
                        if len(MatchingMID_Objects) > 1: #When we have more than one value in the data, find a linear fit
                            #Now we use polyfit, poly1d, and polyval to fit the data linearly and find an approximate ionization factor
                            #TODO: We think we should use a power law with y = mx^b (this implies an intercept of 0 and retains the type of curvature we see in the data)
                            #TODO continued: Link to do so: https://scipy-cookbook.readthedocs.io/items/FittingData.html
                            polynomialCoefficients = numpy.polyfit(MatchingMID_ElectronNumbers,MatchingMID_RS_Values,1) #Electron numbers as the independent var, RS_values as the dependent var, and 1 for 1st degree polynomial
                            poly1dObject = numpy.poly1d(polynomialCoefficients) #create the poly1d object
                            self.relativeIonizationEfficiencies[moleculeIndex] = numpy.polyval(poly1dObject,self.electronnumbers[moleculeIndex]) #use polyval to calculate the ionization factor based on the current molecule's electron number
                            self.sourceOfIonizationData[moleculeIndex] = 'evaluatedInterpolationTypeFit' #ionization factor was determined via a linear fit based on a molecule's ionization type
                if len(MatchingMID_Objects) == 0: #Otherwise use the original Madix and Ko equation
                    self.relativeIonizationEfficiencies[moleculeIndex] = (0.6*self.electronnumbers[moleculeIndex]/14)+0.4        
                    self.sourceOfIonizationData[moleculeIndex] = 'MadixAndKo' #ionization efficiency obtained via Madix and Ko equation

    #Export the ionization efficiencies used and their respective method used to obtain them (known factor, known molecule, known ionization type, or Madix and Ko)    
    def exportIonizationInfo(self):
        ionizationData = numpy.vstack((self.molecules,self.relativeIonizationEfficiencies,self.sourceOfIonizationData)) #make a 2d array containing molecule names (for the header), the ionization efficiencies, and which method was chosen
        ionizationDataAbsicca = numpy.array([['Molecule'],
                                             ['Ionization Efficiency'],
                                             ['Method to Obtain Ionization Efficiency']]) #create the abscissa headers for the csv file
        ionizationDataToExport = numpy.hstack((ionizationDataAbsicca,ionizationData)) #use hstack to obtain a 2d array with the first column being the abscissa headers
        #numpy.savetxt('ExportedIonizationEfficienciesSourcesTypes.csv',ionizationDataToExport,delimiter=',',fmt='%s') #export to a csv file
        numpy.savetxt('ExportedIonizationEfficienciesSourcesTypes.tsv',ionizationDataToExport,delimiter='\t',fmt='%s') #export to a tsv file
        
    def update_relative_standard_uncertainties(self):
        self.relative_standard_uncertainties = self.absolute_standard_uncertainties*1.0 #First make the array.
        #now populate the non-mass fragment parts by dividing.
        #Note that it's possible to get a divide by zero error for the zeros, which we don't want. So we fill those with 0 with the following syntax: np.divide(a, b, out=np.zeros(a.shape, dtype=float), where=b!=0) https://stackoverflow.com/questions/26248654/how-to-return-0-with-divide-by-zero
        a_array = self.relative_standard_uncertainties[:,1:]
        b_array = self.standardized_reference_patterns[:,1:]
        self.relative_standard_uncertainties[:,1:] = numpy.divide(a_array, b_array, out=numpy.zeros(a_array.shape, dtype=float), where=b_array!=0)
        #self.ExportCollector('StandardizeReferencePattern_relative_standard_uncertainties', export_relative_uncertainties= True)

    
    #This makes a full Reference data file, not just a pattern.
    def exportReferencePattern(self, referenceFileName):
        with open(referenceFileName, 'w') as the_file:          
            referenceFileHeader = ''
            referenceFileHeader += "Source:,"  + DataFunctions.arrayLikeToCSVstring(self.SourceOfFragmentationPatterns) + "\n"
            referenceFileHeader += "Molecules," + DataFunctions.arrayLikeToCSVstring(self.molecules) + "\n"
            referenceFileHeader += "Electron Numbers," + DataFunctions.arrayLikeToCSVstring(self.electronnumbers) + "\n"
            referenceFileHeader += "ionizationEfficiencies," + DataFunctions.arrayLikeToCSVstring(self.relativeIonizationEfficiencies) + "\n"
            referenceFileHeader += "ionizationEfficienciesSources," + DataFunctions.arrayLikeToCSVstring(self.sourceOfIonizationData) + "\n"
            referenceFileHeader += "Molecular Mass," + DataFunctions.arrayLikeToCSVstring(self.molecularWeights)# + "\n"
            if hasattr(self ,"standardized_reference_patterns"):
                self.standardized_reference_patterns=StandardizeReferencePattern(self.standardized_reference_patterns,len(self.molecules))
            else:
                self.standardized_reference_patterns=StandardizeReferencePattern(self.provided_reference_patterns,len(self.molecules))
            numpy.savetxt(referenceFileName, self.standardized_reference_patterns.copy(), delimiter=",", header = referenceFileHeader, comments='')
            
'''
The MolecularIonizationData class is used to generate a molecule's ionization factor based on its ionization type
'''        
class MolecularIonizationData (object):
    def __init__(self,moleculeName,RS_Value,electronNumber,moleculeIonizationType='unknown',sourceOfIonizationData='unknown'):
        #Store the MID variables
        self.moleculeName = moleculeName.strip() #TODO: consider making this a list to be consistent with the other objects.
        self.RS_ValuesList = parse.listCast(float(RS_Value)) #Since we can have slightly different RS_values for a molecule, make a list so a molecule with more than one RS_Value can contain all the info provided
        self.electronNumber = float(electronNumber) #TODO: consider making this a list to be consistent with the other objects.
        self.moleculeIonizationType = parse.listCast(moleculeIonizationType)
        self.sourceOfIonizationData = parse.listCast(sourceOfIonizationData) #Different RS values can come from different sources so make a list that will be parallel to RS_ValuesList containing the source of each RS Value at the same index
        
    def addData(self,RS_Value,sourceOfIonizationData):
        #if we have more than one RS_Value, then append to the list
        self.RS_ValuesList.append(float(RS_Value))
        #if we have more than one RS_Value, append the source
        self.sourceOfIonizationData.append(sourceOfIonizationData)																	
############################################################################################################################################
###############################################Algorithm Part 2: Analysing the Processed Data###############################################
############################################################################################################################################
#During trying to debug, the following variables are particularly useful:
#G.massesUsedInSolvingMoleculesForThisPoint keeps track of anything solved by SLSUnique.
#original_list_of_mass_fragments (array of mass fragments) vs. used_mass_fragments (array of 0 and 1 for which ones have been used)
#molecules_unedited (list of strings) vs.   solvedmolecules (array of 0 and 1 for which ones have been solved)
#reciprocal_matching_correction_values is not only a trimmed down version of correction_values, it is also transposed.
#uncertainties_dict contains all of the uncertainties including the dynamic variable versions that get changed during SLSUniqueFragments
#G.currentTimeIndex = timeIndex and there is also G.currentTimePoint #This is useful for debugging when in the for loop of data analysis.                                                                                                                                         
    
#this function compares the list of chosen mass fragments and those monitored and makes a raw signal
#array out of this data, which will be used with the inverse method to find percent signal and composition
#this function discards the mass fragments that were collected but are not present in the reference file
def RawSignalsArrayMaker(mass_fragment_numbers_monitored,mass_fragment_numbers,collected,counter,referenceabscissa):
    mass_fragment_length = len(mass_fragment_numbers)
    rawsignalsarrayline = numpy.zeros([1,1]) #zero used to stack onto the array
    for collectedcounter in range(len(mass_fragment_numbers_monitored)): #array-indexed for loop
        for massfragcounter in range(mass_fragment_length):#array-indexed for loop
            if mass_fragment_numbers[massfragcounter] == mass_fragment_numbers_monitored[collectedcounter]:#if there is a mass fragment number not contained in the mass fragment numbers (made by array builder) then it will not be added
                for referenceabscissacounter in range(len(referenceabscissa)):#array-indexed for loop
                    if referenceabscissa[referenceabscissacounter] == mass_fragment_numbers[massfragcounter]:#checks the reference for this mass fragment as well, before collected data is added
                        rawsignalsarrayline = numpy.vstack([rawsignalsarrayline,collected[counter,massfragcounter]])
    rawsignalsarrayline = numpy.delete(rawsignalsarrayline,(0),axis=0)#deletes zero used to start array building in the loops
    return rawsignalsarrayline


#CombinationMaker gets the combinations of matrices and solves each one by one and enters them into the list of answers-signals
#specifically to make square matrices
#itertools uses a combination function (below) and the function uses those to index drawing out of all the rows in an array
def CombinationMaker(reciprocal_matching_correction_values,rawsignalsarrayline,monitored_reference_intensities,mass_fragment_numbers):
    num_molecules = len(reciprocal_matching_correction_values[0,:])
    num_MassFragments = len(reciprocal_matching_correction_values[:,0])
    import itertools 
    combinations = list(itertools.combinations(list(range(num_MassFragments)),num_molecules)) 
    if len(combinations) == 0:#This function will not work without enough mass fragments, so the user must know the problem
        print('****************************************')
        print('Not enough matching mass fragments input')
        print("This means that at some point in the analysis, there were not enough masses in the reference file to apply the inverse method. It could mean you have too many overlapping masses for the molecules you are trying to resolve.  You can get around this by using the '#//Reference Mass Fragmentation Threshold//' feature to exclude tiny fragementation peaks. This would be done by setting the value to 'yes' for  minimalReferenceValue feature with referenceValueThreshold, such as referenceValueThreshold = 5.0 .  Alternatively, to be more targeted, if you know *which* fragmentation patterns could be overlapping, you could set those minor fragments to 0 in your reference pattern csv file. TODO: Print out the relevant masses here. This requires keeping track of when they are selected prior to combination maker, and possibly passing them as an additional argument.")
        print('****************************************')
    combinations_len = len(combinations) 
    correctionarray = numpy.zeros([1,num_molecules])
    intensityarray = numpy.zeros([1,num_molecules])
    rawsignalarray = numpy.zeros([1,1])
    correctionlist = []
    intensitylist = []
    rawsignallist = [] 
    massfragrow = 'yuh,'
    massfraglist = []
    #this loop gets two lists that contain al the arrays of possible answers and inverse correction values,first rows are made,
    #then they are stacked into an array, this array then, has the first row of zeros deleted, then it is added to an indice
    #of a list where it will later be accessed by the combination solver. three of these lists are output for signals,corrections
    #and relative intensities
    for combinationnum in range(combinations_len): #array-indexed for loop
        for moleculecounter in range(num_molecules):    #array-indexed for loop
            correctionrow = reciprocal_matching_correction_values[combinations[combinationnum][moleculecounter],:] 
            intensityrow = monitored_reference_intensities[combinations[combinationnum][moleculecounter],:]
            rawsignalrow = rawsignalsarrayline[combinations[combinationnum][moleculecounter],:]
            massfragrow = massfragrow + str(mass_fragment_numbers[combinations[combinationnum][moleculecounter]]) + ','
            correctionarray = numpy.vstack([correctionarray,correctionrow]) 
            intensityarray = numpy.vstack([intensityarray,intensityrow])
            rawsignalarray = numpy.vstack([rawsignalarray,rawsignalrow])
            if moleculecounter == num_molecules-1:#the end of the nested loop: the rows just made are entered into arrays
                correctionarray = numpy.delete(correctionarray,(0),axis=0)
                intensityarray = numpy.delete(intensityarray,(0),axis=0)
                rawsignalarray = numpy.delete(rawsignalarray,(0),axis=0)
                massfragrow1 = massfragrow.split(',',1)[1]
                massfraglist.append(massfragrow1)
                massfragrow = 'yuh,'
                correctionlist.append(correctionarray)
                intensitylist.append(intensityarray)
                rawsignallist.append(rawsignalarray)
                correctionarray = numpy.zeros([1,num_molecules])
                intensityarray = numpy.zeros([1,num_molecules])
                rawsignalarray = numpy.zeros([1,1])
    combinations_len = len(combinations)
    return [combinations_len,rawsignallist,correctionlist,intensitylist,massfraglist]


#This function simply solves each of the combinations, drawing the respective values out of the lists and uses numpy.linalg
def CombinationSolver(combinations_len,rawsignallist,correctionlist,molecules,massfraglist):
    compositions = []
    for combinationcounter in range (combinations_len):  #array-indexed for loop
        if numpy.linalg.det(correctionlist[combinationcounter]) != 0:#if the determinant is zero, then doing the linalg.solve function will stop the entire script- so you must use this method
            solutions = numpy.linalg.solve(correctionlist[combinationcounter], rawsignallist[combinationcounter])   #Linalg.solve will solve the matrix equation Ax=B. The first argument is A. The second argument is B. Comparing to the correction factor matrix equation in the madix and Ko excel sheet. A is a matrix of reciprocal correction values. X is a matrix of concentrations. B is a matrix of signals.
            composition = solutions 
            compositions.append(composition)
    return[compositions]
    

#compresses the data into an array of averages and standard deviations- then prints the results
def DataCompressor(signals,molecules,type):
    num_molecules = len(molecules)
    averagegroup = []
    stddevgroup = []
    average = []
    stddev = []
    for moleculecounter in range(num_molecules): #this part of the code is new for this version (3) and it gets avg and stddev
        for combinationnum in range(len(signals)):#array-indexed for loop
            averagegroup.append(signals[combinationnum][moleculecounter])#takes all of the different solutions and puts them in a group
            stddevgroup.append(signals[combinationnum][moleculecounter]) #does exactly what the line above did, with different names
            if combinationnum == len(signals)-1:#the end of the loop
                average.append(sum(averagegroup)/float(len(averagegroup))) #this actually determines the average of each number in the group
                stddev.append(numpy.std(stddevgroup)) #this gets the std dev, in a list
                averagegroup = []
                stddevgroup = []
    return [average,stddev]
    
    
#this function calls all the functions that make up the inverse method, so that there is no need to call them individually later
def InverseMethod(reciprocal_matching_correction_values,rawsignalsarrayline,monitored_reference_intensities,mass_fragment_numbers,molecules,type):
    [combinations_len,rawsignallist,correctionlist,intensitylist,massfraglist] = CombinationMaker (reciprocal_matching_correction_values,rawsignalsarrayline,monitored_reference_intensities,mass_fragment_numbers)
    [compositions] = CombinationSolver (combinations_len,rawsignallist,correctionlist,molecules,massfraglist)
    [averagecomp,stddevcomp] = DataCompressor(compositions,molecules,'composition')
    return averagecomp
    
    

#this function finds the significance of a specified value in an array to the array as a whole
def IndElemSignificanceCalculator(rowDataArray, specifiedColumnIndex, moleculesLikelihood, minThreshold=None, maxIntensityPossible=100):
    length = len(rowDataArray)
    #variable to hold the terms in the summation
    allSummationTerms = []
    #for each value in the array

    if minThreshold == 0: minThreshold=None #the code is written to have minThreshold as None or not None, but 0 counts as None  for this purpose.
    if type(minThreshold) != type(None): #assume the minThreshold is some kind of number, and use it to calculate and cap the indSummationTerm.
        indSummationTermCap = maxIntensityPossible/minThreshold - 1 #note that the 100 is the default maxIntensityPossible, with the assumption that standardized intensities  are capped at 100.
    
    for moleculecounter in range(length):
        if rowDataArray[moleculecounter] != 0: #if the value is not zero, calculate 
            if type(minThreshold) == type(None): #just calculate the indSummationTerm as normal.
                #calculates the unweighted ratio of each value, scaled by the likelihood of that molecule 
                indSummationTerm = abs((moleculesLikelihood[moleculecounter]*rowDataArray[moleculecounter])**float(-1)*(moleculesLikelihood[specifiedColumnIndex]*rowDataArray[specifiedColumnIndex]-1))
                allSummationTerms.append(indSummationTerm)       
            if type(minThreshold) != type(None): #assume the minThreshold is some kind of number, and use it to calculate and cap the indSummationTerm.
                indSummationTerm = abs((moleculesLikelihood[moleculecounter]*rowDataArray[moleculecounter])**float(-1)*(moleculesLikelihood[specifiedColumnIndex]*rowDataArray[specifiedColumnIndex]-1))
                if indSummationTerm > indSummationTermCap: indSummationTerm = indSummationTermCap
                allSummationTerms.append(indSummationTerm)
        if rowDataArray[moleculecounter] == 0: 
            if type(minThreshold) == type(None): 
                pass #the pass is like allSummationTerms.append(0) #if the intensity/value is zero, and we don't have a minThreshold to use, then the final value will be zero as well.
            if type(minThreshold) != type(None): #else assume the minThreshold is some kind of number.
                indSummationTerm = indSummationTermCap #use the cap directly if the term for rowDataArray[moleculecounter] == 0, because that means a denominator of 0 which is like infinity. 
                allSummationTerms.append(indSummationTerm)            
    #the following line can be replace with code such as "significance = (sum(allSummationTerms)**SumCoeffient)*(array[specifiedColumnIndex]**ValueCoefficent)"
    # if you would like to add coefficents to increase or decrease the weighting of each term
    significance = sum(allSummationTerms)*rowDataArray[specifiedColumnIndex]*moleculesLikelihood[specifiedColumnIndex]
    return significance

#This function compiles a list of the significances of each row to a particular column 
def ElemSignificanceCalculator(anArray,specifiedColumnIndex, moleculesLikelihood, minThreshold=None, maxIntensityPossible=100):
    #find the number of rows
    row_num = len(anArray)
    #empty list to store values
    sigValuesList = []
    #for each row...
    for rowcounter in range(row_num):
        # the "Significance" of that row to the column is calculated
        sigValue = IndElemSignificanceCalculator(anArray[rowcounter], specifiedColumnIndex, moleculesLikelihood, minThreshold=minThreshold, maxIntensityPossible=maxIntensityPossible)
        # the significance is stored in a list
        sigValuesList.append(sigValue)
        
    return sigValuesList
#this function sorts a list by its values, but returns the original indicies
#in the sorted order, rather than the sorted values
def ValueOrderedIndexSort(sigValuesList):
     #sort by smallest to largest, return indicies
     orderedSigIndices = numpy.argsort(sigValuesList)
     #flip indicies, so it is now sorted by largest to smallest
     orderedSigIndices = numpy.flipud(orderedSigIndices)
     #return the sorted indicies
     return orderedSigIndices

# this is a sub-function of DistinguishedArrayChooser. It serves to determine, in order, the most significant rows in a data set
# See Charles ProcDoc Summer 2018 page 3 for a mathamatical equation defining significance.     
def ImportantAbscissaIdentifier(anArray, moleculesLikelihood):
    #finding the size of the array
    row_num = len(anArray[:,0])
    column_num = len(anArray[0,:])
    #defining order variables to be passed on
    columnOrderList = [] #one order of rows for each column
    overallOrder = [] #order of rows that should be used for the whole array
    
    #search one column at a time
    for columncounter in range(column_num):
        
        #each row is looped through to find the significance of each point
        sigValuesList = ElemSignificanceCalculator(anArray, columncounter, moleculesLikelihood)
        #the point indicies are ordered by size
        OrderedSignIndices = ValueOrderedIndexSort(sigValuesList)
        
        #store the row indexs in order of most significant to least
        columnOrderList.append(OrderedSignIndices)
        #if I am trying to add the row for the first column, then I can do so 
        #without worry because it cannot possibly be used already
        if len(overallOrder) == 0:
            overallOrder.append(columnOrderList[columncounter][0])
        #if I am trying to add a row for a later column to the overall order
        else:
            #then I need to be careful that it isn't used before 
            #I do this by checking it against the added rows
            #for each row that could be added
            for row_counter in range(row_num):
                #if that row hasen't been used yet
                if not columnOrderList[columncounter][row_counter] in overallOrder:
                    #Then we append that row
                    overallOrder.append(columnOrderList[columncounter][row_counter])
                    #we are also done with that column, so we can go back up to the
                    #beginning of the function for the next column
                    break
    #return the order of rows for the array                
    return overallOrder
                

#List value checker simply confirms that a list matches the desired length 
#If the lengths don't match, list checker will edit them, so that they do match
#If non-obvious editing is required, a warning will also be printed to the user
def ListLengthChecker(aList, desiredLength, defaultNum):
    #true default: if user hasn't entered any values to the list
    if len(aList) == 0:#if the value is not in the data edit file
        aList = [defaultNum] * desiredLength
    #Perfect operation: user has provided the correct number of values to the list
    elif len(aList) == desiredLength:
        pass
    #if the list is one entry long, simply apply that value to all mssing spots
    elif len(aList) == 1:
        firstSensValue = aList[0]
        aList = [firstSensValue] * desiredLength
    #all other cases: warn user and simply use the first value for all entries
    else:
        firstSensValue = aList[0]
        aList = [firstSensValue] * desiredLength
        print("Warning, the distinguished inverse specifications that you have provided are of a different length than the number of molecules that you have provided.")
    return aList
    
#this function is going to be used by multiple sections of the code, including the updated sls method and a secondary inverse method
#this is a new way of selecting the most important rows, for each molecule, based on that molecules ratios with the other molecules 
#in that row and that molecules own value.
#uncertainties_dict should be an empty dictionary when uncertainties will not be used.
def DistinguishedArrayChooser(refMassFrags,correctionValues,rawSignals,moleculeLikelihoods,sensitivityValues, uncertainties_dict={}):
    #the shape of the referenceData is found 
    num_rows = len(refMassFrags[:,0])   #This is number of mass frags, if it's mass spec data.
    num_columns = len(refMassFrags[0,:]) #This is number of molecules, if it's mass spec data.
    
    #The acceptable threshold is determined by the SensitivityValue function
    sensitivityValues = ListLengthChecker(sensitivityValues, num_columns, 1)
   
    #the moleculesLikelihood is corrected if it wasn't entered by the user.
    moleculeLikelihoods = ListLengthChecker(moleculeLikelihoods, num_columns, 1)
    
    #all values below the specified relative intensity must be set the minThreshold value
    #This is because a subfunction attempts to divide by each value
    for columncounter in range(num_columns):
        for rowcounter in range(num_rows):
            if refMassFrags[rowcounter,columncounter] < sensitivityValues[columncounter]: 
                refMassFrags[rowcounter,columncounter] = 0 #sensitivityThresholdValue[0]
                
    #The correct order of the needed rows is determined by this function
    order = ImportantAbscissaIdentifier(refMassFrags,moleculeLikelihoods)
    
    #empty lists to store results i.e. shortened arrays
    shortRefMassFrags = []
    shortCorrectionValues = []
    shortCorrectionValues_relative_uncertainties = []
    shortRawSignals = []
    shortRawSignals_absolute_uncertainties = []
    
    #add the correct row to each list. I believe the row_num here is the mass fragment index..
    for row_num in order:
        shortRefMassFrags.append(refMassFrags[row_num])
        shortCorrectionValues.append(correctionValues[row_num])
        shortRawSignals.append(rawSignals[row_num])
        if len(uncertainties_dict) > 0:
            shortCorrectionValues_relative_uncertainties.append(uncertainties_dict['reciprocal_matching_correction_values_relative_uncertainties_one_time'][row_num])
            shortRawSignals_absolute_uncertainties.append(uncertainties_dict['rawsignals_absolute_uncertainties_one_time'][row_num])
    #This section stacks the chosen rows from lists into arrays
    shortRefMassFrags = numpy.asarray(shortRefMassFrags)
    shortCorrectionValues = numpy.asarray(shortCorrectionValues)
    shortRawSignals = numpy.asarray(shortRawSignals)
    
    #finding the size of the new array
    num_rows = len(shortRefMassFrags[:,0])
    num_columns = len(shortRefMassFrags[0,:])
    #This section replaces the minThreshold's that were chosen as threshold values with 0s
    for columncounter in range(num_columns):
        for rowcounter in range(num_rows):
            if shortRefMassFrags[rowcounter,columncounter] <= sensitivityValues[columncounter]: 
                shortRefMassFrags[rowcounter,columncounter] = 0
                
    #The shortened arrays are finally returned to the Inverse Method solver                
    return shortRefMassFrags,shortCorrectionValues,shortRawSignals, shortCorrectionValues_relative_uncertainties, shortRawSignals_absolute_uncertainties
    
    
#this function takes the data from important abscissa identifier and 
def InverseMethodDistinguished(monitored_reference_intensities,reciprocal_matching_correction_values,rawsignalsarrayline, uncertainties_dict={}):
    monitored_reference_intensities,reciprocal_matching_correction_values,rawsignalsarrayline, reciprocal_matching_correction_values_relative_uncertainties, rawsignals_absolute_uncertainties = DistinguishedArrayChooser(monitored_reference_intensities,reciprocal_matching_correction_values,rawsignalsarrayline, G.moleculeLikelihoods,G.sensitivityValues, uncertainties_dict)
    #The below try and except statemnt is meant to catch cases as described in the except statement.
    try:
        numpy.linalg.det(reciprocal_matching_correction_values)
    except:
        print("There is an error in a matrix operation evaluation: The number of feasible mass fragments to check is probably less than the number of molecules. This can happen if referenceValueThreshold is too strict, leaving not enough feasible fargments to consider. The program is probably about to crash.")               
    if numpy.linalg.det(reciprocal_matching_correction_values) != 0:#only solves if determinant is not equal to zero     
        #Now we will solve with the inverse way, and also for uncertainties if that module is present.
        #The uncertainties module works well for matrix inverse, but not for dot product. So we make a custom function.                                                                                                   
        if len(uncertainties_dict) > 0:
            try:
                import uncertainties
                from uncertainties import unumpy
                uncertaintiesModulePresent = True
                # Below function works but does not work correctly with uncertainties module. It gives the wrong answer when you use it, as of Feb 2nd 2020
                # def loopBasedDotProductComma(matrix_a, matrix_b, verbose = False):# explicit for loops 
                    # #the result will have the same number of rows as the 1st matrix, and the same number of columns as the 2nd matrix.
                    # #The number of columns of the 1st matrix must equal the number of rows of the 2nd matrix.    
                    # #print("here 10", matrix_a)
                    # num_a_rows = len(matrix_a)
                    # num_a_columns = len(matrix_a[0])
                    # num_b_rows = len(matrix_b)
                    # num_b_columns = len(matrix_b[0]) 
                    # numFinalRows = num_a_rows
                    # numFinalColumns = num_b_columns
                    # dotProductArray = numpy.zeros((numFinalRows,numFinalColumns), dtype="object") #Make it an object array so it can hold anything.
                    # #row_a and row_b below are actually indices, not rows. The wrod index is omitted for brevity.
                    # #We are going to fill the final matrix top row first, then second row next, etc.
                    # print(num_a_rows, num_a_columns, num_b_rows, num_b_columns,  numFinalRows, numFinalColumns, numpy.shape(dotProductArray))
                    # for finalColumn in range((numFinalColumns)):
                        # #We're going to loop down each column as well.
                        # for finalRow in range((numFinalRows)):        
                            # #Within the first row (each row) of the product matrix... 
                            # #    The first item gets filled from the first row of matrix_a and the corresponding column of matrix_b.
                            # #    The second item gets filled from the second row of matrix_a and the corresponding column of matrix_b. 
                            # #We're going to loop across a_columns (which is same number as b_rows) to make the temporary sum.
                            # temporarySum = 0
                            # for additionIndex in range((num_a_columns)):
                                # if verbose:
                                    # print("line 2839",  matrix_a[finalRow,additionIndex])
                                    # print(finalColumn, finalRow, matrix_a[finalRow,additionIndex], matrix_b[additionIndex,finalColumn])
                                # #NEEDED TO USE COMMAS. THAT WAS CAUSING PROBLEMS FOR MATRICES.
                                # temporarySum = temporarySum + matrix_a[finalRow,additionIndex]*matrix_b[additionIndex,finalColumn]
                            # if verbose:
                                # print(temporarySum)
                                # print(temporarySum.n)
                                # print(temporarySum.s)
                            # temporarySum =uncertainties.ufloat(temporarySum.n, temporarySum.s) #Was not working properly, so had to pull n and s out and cast. (nominal value and std_dev value)
                            # dotProductArray[finalRow,finalColumn]= temporarySum
                    # return dotProductArray 
            except:
                uncertaintiesModulePresent = False
                solutions = numpy.linalg.solve(reciprocal_matching_correction_values,rawsignalsarrayline) #Linalg.solve will solve the matrix equation Ax=B. The first argument is A. The second argument is B. Comparing to the correction factor matrix equation in the madix and Ko excel sheet. A is a matrix of reciprocal correction values. X is a matrix of concentrations. B is a matrix of signals. 
                uncertainties_dict['concentrations_absolute_uncertainties_one_time'] = numpy.array(solutions*0.0).transpose() #The uncertainties dictionary does have these the other way, so have to transpose back.
                uncertainties_dict['concentrations_relative_uncertainties_one_time'] = numpy.array(solutions*0.0).transpose() #The uncertainties dictionary does have these the other way, so have to transpose back.
                print("WARNING: InverseMethodDistinguished could not return uncertainties and is returning uncertainties of zero. This is for the datapoint with intensities of", rawsignalsarrayline, "and concentrations of", solutions)
            if uncertaintiesModulePresent == True:              
                reciprocal_matching_correction_values_relative_uncertainties = reciprocal_matching_correction_values_relative_uncertainties#Note that we cannot use uncertainties_dict['reciprocal_matching_correction_values_relative_uncertainties'], because we have truncated the array already based on distinguished.
                reciprocal_matching_correction_values_absolute_uncertainties = reciprocal_matching_correction_values*reciprocal_matching_correction_values_relative_uncertainties                                                        
                unumpy_reciprocal_matching_correction_values = unumpy.uarray(numpy.array(reciprocal_matching_correction_values), numpy.array(reciprocal_matching_correction_values_absolute_uncertainties)) #inverse doesn't work if they're not numpy arrays already. They were lists before this line.
                #Needs to be a umatrix, not uarray, if want to invert and do dot product.          
                #following https://pythonhosted.org/uncertainties/numpy_guide.html?highlight=inverse then adding dot product.
                uInverseInBetween = unumpy.ulinalg.pinv(unumpy_reciprocal_matching_correction_values)                                  
                #Up until this line, rawsignals_absolute_uncertainties has been 1D. We need to make it 2D and the same direction as rawsignalsarrayline which is already 2D.
                rawsignals_absolute_uncertainties_2D = numpy.atleast_2d(rawsignals_absolute_uncertainties).transpose()
                unumpy_rawsignalsarrayline = unumpy.uarray(numpy.array(rawsignalsarrayline), numpy.array(rawsignals_absolute_uncertainties_2D)) #When this has 0.0 we get a very small final uncertainties listed as 0.        
                uSolutions = uInverseInBetween@unumpy_rawsignalsarrayline #the @ forces matrix multiplication.
                solutions= unumpy.nominal_values(uSolutions)
                solutions_uncertainties =unumpy.std_devs(uSolutions)
                uncertainties_dict['concentrations_absolute_uncertainties_one_time'] = abs(numpy.array(solutions_uncertainties).transpose()) #The uncertainties dictionary does have these the other way, so have to transpose back.
                uncertainties_dict['concentrations_relative_uncertainties_one_time'] = abs(numpy.array(solutions_uncertainties/solutions).transpose()) #The uncertainties dictionary does have these the other way, so have to transpose back.
        solutions = numpy.linalg.solve(reciprocal_matching_correction_values,rawsignalsarrayline)   #Linalg.solve will solve the matrix equation Ax=B. The first argument is A. The second argument is B. Comparing to the correction factor matrix equation in the madix and Ko excel sheet. A is a matrix of reciprocal correction values. X is a matrix of concentrations. B is a matrix of signals.
    else:
        print('The Array Chosen is Singular')
        solutions = numpy.zeros(len(rawsignalsarrayline)) # the solutions are         if len(uncertainties_dict) > 0:
        if len(uncertainties_dict) > 0:
            print('The Array Chosen is Singular. Using a 0 for all uncertainties. In future, may try to use pseudo inverse uncertainties propagation, but it would likely underestimate the uncertainties.')
            #In future can probably use uInverseInBetween = unumpy.ulinalg.pinv(unumpy_reciprocal_matching_correction_values), assuming that function exists.
            uncertainties_dict['concentrations_absolute_uncertainties_one_time']=solutions*0.0
            uncertainties_dict['concentrations_relative_uncertainties_one_time']=solutions*0.0        
    return solutions

    
#this function is going to get the data from the array that is being made (holding all the solutions thus far) but using the molecules
#copy from the main sls method in order know the order the full array has the molecules in
def MoleculeRange(molecules,timeIndex,molecules_copy,scaledConcentrationsarray):
    moleculedata = []
    moleculedata1 = molecules_copy
    #on the first time, the signalsarray that is given is actually just the predicted values for all the molecules given
    #so that the list is already read and does not need to be modified
    if timeIndex == 0:
        moleculedata = scaledConcentrationsarray
    elif timeIndex == 1:#after the first time, so there is only one row on the array with the chosen signals
        moleculedata2 = scaledConcentrationsarray[1:]
        #once the two lists of strings are formatted so that the molecules that are the same are exactly the same, the two lists are compared so
        #that the values that were used for the molecules that need ranges are extracted and used in the next function
        for moleculecounter in range(len(molecules)):#array-indexed for loop
            for moleculedatacounter in range(len(moleculedata1)):#array-indexed for loop
                if molecules[moleculecounter] == moleculedata1[moleculedatacounter]:#get index for correction molecule
                    if moleculedata2[moleculedatacounter] < 0.00001:
                        moleculedata.append(0)
                    else:
                        moleculedata.append(moleculedata2[moleculedatacounter])
    elif timeIndex > 1:#after that, you must index for the row too
        moleculedata2 = scaledConcentrationsarray[timeIndex-1,1:]
        #once the two lists of strings are formatted so that the molecules that are the same are exactly the same, the two lists are compared so
        #that the values that were used for the molecules that need ranges are extracted and used in the next function
        for moleculecounter in range(len(molecules)):#array-indexed for loop
            for moleculedatacounter in range(len(moleculedata1)):#array-indexed for loop
                if molecules[moleculecounter] == moleculedata1[moleculedatacounter]:#get index for correction molecule
                    if moleculedata2[moleculedatacounter] < 0.00001:
                        moleculedata.append(0)
                    else:
                        moleculedata.append(moleculedata2[moleculedatacounter])
    return moleculedata

#TODO Ashi needs to check that dataRangeSpecifier is working
#this function here is called in the section at the bottom, right before BruteForce is called, because this function
#sets the ranges and increments of search for the brute force function below. It does so by using prompts within a for
#loop to get the upper and lower bounds and lastly the increments, all being added to a tuple which is then put into a 
#list with all the other tuples, to be used as the ranges for the individual molecules
def DataRangeSpecifier(molecules,timeIndex,molecules_copy,conversionfactor,datafromcsv,DataRangeSpecifierlist,scaledConcentrationsarray):
    [dataRangeSpecifierYorN,signalOrConcentrationRange,csvFile,moleculesRange,csvFileName,higherBound,lowerBound,increments,permutationNum] = DataRangeSpecifierlist
    
    #TODO: #FIXME: DataRangeSpecifier has a bunch of hard coded values like
    # 0.5, 2.0, and 0.01. These need to be fixed, eventually.
    
    #the user can choose to input a set of values for the entire times or to input a .csv file full of values so that 
    #brute only checks in that range
    if dataRangeSpecifierYorN == 'yes':
        holder = numpy.zeros(len(molecules))
        for moleculeIndex in range(len(molecules)):#array-indexed for loop
            for rangeIndex in  range(len(moleculesRange)):#array-indexed for loop
                if moleculesRange[rangeIndex] == molecules[moleculeIndex]:#array of zeros and ones made where ones are molecules included in data edit file
                    holder[moleculeIndex] = 1
        place_holder = 0
        for moleculeIndex in range(len(holder)):#array-indexed for loop
            if holder[moleculeIndex] == 1: #if there are ones, those molecules are deleted from the molecules array- so this function can make ranges for them
                molecules = numpy.delete(molecules,moleculeIndex-place_holder)
                place_holder = place_holder + 1
        if csvFile == 'no':# This part is for non csv file ranges- single ranges for the whole times
            specifications1 = []
            #The conversion factor is used so that the user's specified molecular concentration boundaries can be convertd into boundaries for the scaled to CO concentrations that are solved for
            if lowerBound != []:#if they are not empty (they should not be since the user input 'yes'), then they are set as the permanent ranges 
                if signalOrConcentrationRange == 'concentration':#if the user input concentration ranges, the ranges are divided by the conversion factor- because this function uses raw signals
                    lowerBound = lowerBound/float(conversionfactor)
                    higherBound = higherBound/float(conversionfactor)
                    increments = increments/float(conversionfactor)
                for moleculecounter in range(len(moleculesRange)):#array-indexed for loop, adds the input values to specifications1 - user inputs
                    specifications1.append(tuple([lowerBound[moleculecounter],higherBound[moleculecounter],increments[moleculecounter]]))
        else: #if the csv file is used
            specifications1 = []
            datafromcsv = numpy.array(datafromcsv)[timeIndex,1:]
            if signalOrConcentrationRange == 'concentration':#divides all concentrations by conversion factor
                datafromcsv = datafromcsv/float(conversionfactor)
            for moleculecounter in range(0,len(moleculesRange)*3,3):#gets three values for every molecule, so the moleculeIndex skips every three, the indexes below selected each of those three with +0,+1 and +2.
                specifications1.append(tuple([datafromcsv[moleculecounter],datafromcsv[moleculecounter+1],datafromcsv[moleculecounter+2]]))
        permutations1 = 1
        for specIndex in range(len(specifications1)):#array-indexed for loop, adds product of iteration numbers, to be used in the next little section
            permutations1 = permutations1*(specifications1[specIndex][1]-specifications1[specIndex][0])/float(specifications1[specIndex][2])
    else:#if there are no user inputs, then permutations will still divide a certain number, as such it must be one to leave the number unchanged
        permutations1 = 1
    #in this section the number of combinations will be limited; in this way there will either more or less values iterated
    #across for each molecule; it turns out that just multiplying whatever number of permutations that you would like to limit
    #the function to put to the power of the inverse of the 
    permutationNum = permutationNum/float(permutations1)
    iterationnumber = int(numpy.floor(permutationNum**(float(len(molecules))**(-1))))
    #the ranges are defined by the previous set of data after the moleculeIndex gets to one, this way they are known to be 
    #within a certain range of each other, thanks to the Interpolater. In this way, the first round is calculated by 
    #the inverse function if the values aren't defined by the user. If the previous value was zero, the function checks
    #5 values between 0 and 0.1. This will probably become a possible user input.
    rangeslist = []
    specifications = []
    moleculedata = MoleculeRange(molecules,timeIndex,molecules_copy,scaledConcentrationsarray)
    for moleculedatacounter in range(len(moleculedata)):#array-indexed for loop
        if moleculedata[moleculedatacounter] == 0:#if the value was zero the time before, the function uses the default range
            spread = numpy.linspace(0.0,0.01,iterationnumber)
            difference = spread[1]-spread[0]
            rangeslist.append((0,0.01,difference))
        else:
            spread = numpy.linspace(moleculedata[moleculedatacounter]*0.5,moleculedata[moleculedatacounter]*2.0, iterationnumber)
            difference = spread[1]-spread[0]
            rangeslist.append((moleculedata[moleculedatacounter]*0.5,moleculedata[moleculedatacounter]*2,difference))
        specifications = rangeslist
    #this is the second part of the user input option, as they may only select ranges for certain molecules, and this
    #will concatenate the results which were calculated based on the previous time's signals. 
    if dataRangeSpecifierYorN == 'yes':#answer in data edit file
        place_holder = 0
        place_holder2 = 0
        specifications2 = []
        for moleculeIndex in range(len(holder)):#array-indexed for loop
            if holder[moleculeIndex] == 0:#these concatenate the specification from user input and calculated from past answers
                specifications2.append(specifications[place_holder])
                place_holder = place_holder + 1
            else:
                specifications2.append(specifications1[place_holder2])
                place_holder2 = place_holder2 + 1
        specifications = specifications2
    return specifications
    
    
#this is the forward function for creating the different possible raw signals, but what it does is multiply the matrix
#of the correction values by that of the generated percentages, subtracting the actual raw signals from this simulated
#array of raw signals, thereby making the distance from zero the error, which is now contained in an array. These values
#are each squared, then summed in order to get the SSR, which the scipy.optimize.brute then gets the minimum of, thereby
# finding the percentages that are closest to the correct percentages
def SimulateSignalAndReturnObjFunc(concentrations,*otherArgumentsList):
    rawsignalsarrayline = otherArgumentsList[0]#the item in the list is the raw signal arrayline
    reciprocal_matching_correction_values = otherArgumentsList[1] #the second item is the matching correction values
    objectiveFunctionOption = otherArgumentsList[2]#the input argument contains both the objectiveFunctionOption (objectiveFunctionType) and the two arrays: raw signals array line and matching correction values
    xyyData = numpy.zeros([3,len(rawsignalsarrayline)]) #a three line array is made that will be entered into the function below
    xyyData[1:2,:] = numpy.hstack(rawsignalsarrayline) #this second line is the raw signals
    xyyData[2:3,:] = numpy.hstack(numpy.array(numpy.array(reciprocal_matching_correction_values)@numpy.array(numpy.vstack(concentrations)))) #the third row is the calculated signals
    objectiveFunctionDictionary = ObjectiveFunctionGenerator(xyyData,0)
    if objectiveFunctionOption == 'weightedSAR': #based on the choice given the output will be chosen from this called functions dictionary
        objective_function = objectiveFunctionDictionary['weightedSAR']
    if objectiveFunctionOption == 'weightedSSR':
        objective_function = objectiveFunctionDictionary['weightedSSR']
    if objectiveFunctionOption == 'ssr':
        objective_function = objectiveFunctionDictionary['ssr']
    if objectiveFunctionOption == 'sar':
        objective_function = objectiveFunctionDictionary['sar']
    return objective_function

    
###The ObjectiveFunctionGenerator function takes in either a 2D list or array of data in xyy format with the 0th index containing
#the event numbers, the 1st index containing the observed data, and the 2nd index containing the predicted data, and args.
#The args is a list containing first the absoluteError, second the minVal, and third the percentError. You must always provide 
#the absoluteError and minVal, but the percent Error is optional. The function calculates the weighting factors, residuals, 
#raw and weighted absolute
# residuals and square residuals, sum of the weighted absolute residuals, sum of the weighted square residuals, and sum 
#of the raw square residuals. It outputs a dictionary containing the sum of weighted absolute residuals, sum of weighted square
#residuals, and sum of raw square residuals.#####
def ObjectiveFunctionGenerator(xyyData, minVal, args=[1.,1.], maxValOption=0):
    dataArray = numpy.asarray(xyyData)
    percentError = args[1]
    weightingFactorArray = numpy.ones(len(dataArray[0])) #Initializing a weighting factor array same length as the data
    maxAbsValue = max(abs(dataArray[1]))
    for row in range(len(dataArray[1])): #Array-indexing for loop
        value = dataArray[1][row]
        if type(args[0]) == type(float(1.0)): #If the absolute error argument is a float, then there is only one value and 
        #it is used throughout the entire weighting process
            absoluteError = args[0]
        else:
            if len(args[0]) == 1: #If the absoluteError array only has one value, then that single value is used
            #throughout the entire weighting process
                absoluteError = args[0]
            elif len(args[0]) > 1: #If the absoluteError array has an absolute error for each data point, then this
            #tells the loop to move through the absolute errors concurrent to the 
                absoluteError = args[0][row]
        ##The below section calculates the weighting factor        
        if maxValOption == 0:
            if value <= minVal: #This avoids dividing by zero, since calculating the weighting factor involves dividing by
            #the value
                if minVal == 0.:
                    weightingFactor = 0.
                else:
                    weightingFactor = (1/minVal)*(1/max((value*percentError), absoluteError))
            else:
                weightingFactor = (1/abs(value))*(1/max((value*percentError), absoluteError))
        if maxValOption==1:
            if value <= minVal: #This avoids dividing by zero, since calculating the weighting factor involves dividing by
                    #the value
                if minVal == 0.:
                    weightingFactor = 0.
                else:
                    weightingFactor = (1/maxAbsValue)*(1/max((value*percentError), absoluteError))
            else:
                weightingFactor = (1/maxAbsValue)*(1/max((value*percentError), absoluteError))
        weightingFactorArray[row] = weightingFactor
    absResidualsArray = numpy.zeros([len(dataArray[0])]) #Initializing an absolute residuals array same length as the data
    for row in range(len(dataArray[1])): #Array-indexing for loop
        observedValue = dataArray[1][row]
        predictedValue = dataArray[2][row]
        residualVal = observedValue - predictedValue #Calculates residual
        absRes = abs(residualVal)
        absResidualsArray[row] = absRes
    sqResArray = absResidualsArray*absResidualsArray
    sar = sum(absResidualsArray)
    ssr = sum(sqResArray) #Calculates the sum of the square residuals (ssr)
    weightedAbsResArray = numpy.zeros([len(dataArray[0])]) #Initializing a weighted absolute residuals array same length as data
    for row in range(len(dataArray[1])): #Array-indexing for loop
        weightedAbsRes = weightingFactorArray[row]*absResidualsArray[row] #Weights absolute residuals
        weightedAbsResArray[row] = weightedAbsRes
    weightedSAR = sum(weightedAbsResArray) #Calculates the sum of the weighted absolute residuals
    weightedSqResArray = weightingFactorArray*sqResArray
    weightedSSR = sum(weightedSqResArray) #Calculates the sum of the weighted square residuals 
    objectiveFunctionDictionary = {'weightedSAR':weightedSAR, 'weightedSSR':weightedSSR, 'ssr':ssr, 'sar':sar} 
    return objectiveFunctionDictionary
    

#This is the function that actually calls the function brute and runs it based on the function above, with ranges as specified
#by the datarangespecifier function, it proceeds to output the answers as well as print them to the screen.
def OptimizingFinisher(molecules,specifications,reciprocal_matching_correction_values,rawsignalsarrayline,objectiveFunctionType,maxPermutations=100001, optimizer="BruteForce"):
    #The optimizer can be BruteForce (grid based) or can come from scipy.optimize.minimize
    #specifications is normally a list. If the optimizer is BruteForce, then the specifications variable contains information about the grid. Index 1 and 0 are max and min along an axis, index 2 is increment/spacing along that axis.
    #If the optimizer is not BruteForce, then specifications includes the initial guess in its first index.
    SimulateSignalAndReturnObjFuncArgs = []#a list is made in order to use inside the forward problem since we are not defining anything globally
    SimulateSignalAndReturnObjFuncArgs.append(rawsignalsarrayline)
    SimulateSignalAndReturnObjFuncArgs.append(reciprocal_matching_correction_values)
    SimulateSignalAndReturnObjFuncArgs.append(objectiveFunctionType)
    def SimulateSignalAndReturnObjFuncWrapper(concentrationGuess):
        objFuncValue = SimulateSignalAndReturnObjFunc(concentrationGuess,*SimulateSignalAndReturnObjFuncArgs)        
        return objFuncValue
    
    from scipy import optimize
    if optimizer=="BruteForce":
        range_by_subtraction = numpy.zeros(len(specifications))
        quotient = numpy.zeros(len(specifications))
        
        product = 1 #Here, product is initialized.
        for specificationnumber in range(len(specifications)):
            range_by_subtraction[specificationnumber] = specifications[specificationnumber][1]-specifications[specificationnumber][0] #subtraction between max and min along axis.
            quotient[specificationnumber] = range_by_subtraction[specificationnumber]/float(specifications[specificationnumber][2])  #range_by_subtraction is divided by increment/spacing desired along that axis.
            product = quotient[specificationnumber]*product #Here, product is multiplied by itself and gets bigger and bigger.

        if product < maxPermutations:
            resbrute = optimize.brute(SimulateSignalAndReturnObjFunc, specifications, args = SimulateSignalAndReturnObjFuncArgs, full_output=True, finish=None)#calls our forward function SignalSimulation

        else:
            print('Error: Too many Permutations')
            sys.exit()
        answers = resbrute[0]
    if optimizer!="BruteForce":
        initialGuess = specifications[0]
        answers=initialGuess*1.0 #The default value is what we were provided as an initial guess.
        SimulateSignalAndReturnObjFuncWrapper(initialGuess)
        optimizedResultObject = optimize.minimize(SimulateSignalAndReturnObjFuncWrapper, initialGuess, method=optimizer)
        optimizedAnswers = optimizedResultObject.x #This "x" is how to pull the results out from scipy optimize minimize  results object.
        answers = optimizedAnswers
        #answers = initialGuess
    return answers

def excludeEmptyMolecules(remaining_num_molecules, solutions, solvedmolecules, remaining_monitored_reference_intensities, remaining_reciprocal_correction_factors_SLS, remaining_reference_intensities_SLS, remaining_molecules_SLS, molecules_unedited, uncertainties_dict={}):   
    #initialize a variable for moleculeIndex before the loop across all molecules.
    moleculeIndexIncludingDeletions = 0
    for moleculeIndex in range(remaining_num_molecules):#array-indexed for loop. Ideally, we'll do SLS once for each molecule.
        referenceIntensitiesForThisMolecule = remaining_monitored_reference_intensities[:,moleculeIndex]  #Note that this must be monitored_reference_intensities, so it has same indexing as moleculeIndex even while remaining_reference_intensities_SLS gets shortened.
        if sum(referenceIntensitiesForThisMolecule) == 0.0:
            #Find the name of the molecule getting deleted.
            nameOfThisMolecule = remaining_molecules_SLS[moleculeIndexIncludingDeletions]
            #Find the index of where that molecule was originally, to update solvedmolecules accordingly, and also solutions accordingly.
            originalMolecularIndex = list(molecules_unedited).index(nameOfThisMolecule)
            #this is setting concentration to 0 for that molecule, in this iteration.
            solutions[originalMolecularIndex] = 0.0  #note that we use the actual moleculeIndex here.
            #No need to subtract any signals
            #update the used molecules list, and amounts remaining for other things.
            solvedmolecules[originalMolecularIndex] = 1 #note that we use the actual moleculeIndex here. 
            remaining_reciprocal_correction_factors_SLS = numpy.delete(remaining_reciprocal_correction_factors_SLS,(moleculeIndexIncludingDeletions),axis = 1)
            remaining_reference_intensities_SLS = numpy.delete(remaining_reference_intensities_SLS,(moleculeIndexIncludingDeletions),axis = 1)
            remaining_molecules_SLS = numpy.delete(remaining_molecules_SLS,(moleculeIndexIncludingDeletions))      
            remaining_num_molecules = remaining_num_molecules -1
            moleculeIndexIncludingDeletions = moleculeIndexIncludingDeletions - 1 
            if len(uncertainties_dict) > 0:
                if type(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']) != type(None):
                    uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'] = numpy.delete(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'],(moleculeIndexIncludingDeletions),axis = 1)
        #increase the index no matter what, to go to next part of loop.
        moleculeIndexIncludingDeletions = moleculeIndexIncludingDeletions + 1
    #connect to return variables after the loop:
    remaining_correction_factors_SLS_after_exclusion = remaining_reciprocal_correction_factors_SLS
    remaining_reference_intensities_SLS_after_exclusion = remaining_reference_intensities_SLS
    remaining_molecules_SLS_after_exclusion = remaining_molecules_SLS
    return remaining_num_molecules, solutions, solvedmolecules, remaining_correction_factors_SLS_after_exclusion, remaining_reference_intensities_SLS_after_exclusion, remaining_molecules_SLS_after_exclusion, uncertainties_dict
    
#this function is a path is sequential linear subtraction, which can be used alongside the inverse
#method or as opposed to it. Either way, it only produces one set of values, so it has no need for the 
#data compressor function and starts right after the correction values are obtained
#TODO: make some kind of unit test that tests a good order being chosen.
def SLSUniqueFragments(molecules,monitored_reference_intensities,reciprocal_matching_correction_values,rawsignalsarrayline, timeIndex, time, uncertainties_dict={}):
    #FIXME: I am using the currentReferenceData.mass_fragment_numbers_monitored but it needs to be passed in from Reference or Experimental datal.
    
    combinedString = '' #This line is for debugging purposes. 
    try:
        type(G.massesUsedInSolvingMoleculesForThisPoint)
    except: #This means somebody is trying to call the function directly and does not need G.massesUsedInSolvingMoleculesForThisPoint.  Probably due to using a unit test or best mass frag chooser.
        G.massesUsedInSolvingMoleculesForThisPoint = numpy.zeros(len(molecules))
    original_list_of_mass_fragments = copy.deepcopy(currentReferenceData.mass_fragment_numbers_monitored)
    # This is creating a local copy of 'monitored_reference_intensities' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reference_intensities_SLS = copy.deepcopy(monitored_reference_intensities)
    

    # This is creating a local copy of 'reciprocal_matching_correction_values' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reciprocal_correction_factors_SLS = copy.deepcopy(reciprocal_matching_correction_values)
    original_correction_factors = copy.deepcopy(reciprocal_matching_correction_values)
    #We keep track of the original number of mass fragments and molecules for indexing purposes.
    #In the loop, we'll update the number of remaining molecules and  mass fragments from remaining correction values.
    original_num_MassFragments = len(original_correction_factors[:,0])
    original_num_molecules = len(original_correction_factors[0,:]) 
    remaining_num_molecules = original_num_molecules #This is how to start the variable.

    # This is creating a local copy of 'rawsignalsarrayline' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_rawsignals_SLS = copy.deepcopy(rawsignalsarrayline)

    # This is creating a local copy of 'molecules' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_molecules_SLS = copy.deepcopy(molecules)
    molecules_unedited = copy.deepcopy(molecules) #old variable, but being kept to prevent need to change things.
    
    if len(uncertainties_dict) > 0: #This means that the uncertainties_dict argument has been passed in with values to use.
        uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'] = copy.deepcopy(uncertainties_dict['reciprocal_matching_correction_values_relative_uncertainties_one_time'])
        if (type(uncertainties_dict['rawsignals_absolute_uncertainties']) != type(None)) and (str(uncertainties_dict['rawsignals_absolute_uncertainties']).lower() != 'none'):
            uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'] = copy.deepcopy(uncertainties_dict['rawsignals_absolute_uncertainties_one_time'])
            #need to also make it 2D and then transpose.
            uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'] = numpy.atleast_2d(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS']).transpose()
        elif type(uncertainties_dict['rawsignals_absolute_uncertainties']) == type(None): #if it's a None type, we need to make a zeros array because we *will* need to populate it.
            uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'] = remaining_rawsignals_SLS*0.0 #This will give the correct size and shape.
        if not 'concentrations_relative_uncertainties_one_time' in uncertainties_dict:
            uncertainties_dict['concentrations_relative_uncertainties_one_time'] = numpy.zeros(len(molecules))
        if not 'concentrations_absolute_uncertainties_one_time' in uncertainties_dict:
            uncertainties_dict['concentrations_absolute_uncertainties_one_time'] = numpy.zeros(len(molecules))

    #initiailzing some lists for population (and truncation).
    nonZeroCorrectionValuesList = []
    signalsCorrespondingToNonZeroCorrectionValuesList = []
    solvedmolecules = numpy.zeros(len(remaining_molecules_SLS))
    #To keep track of the mass fragments used. (just an array of zeros to start, then fill with ones)
    used_mass_fragments= numpy.array(copy.deepcopy(original_list_of_mass_fragments))*0.0
    last_solvedmolecules = solvedmolecules*1.0  # This is just initializing for use later.
    last_used_mass_fragments = used_mass_fragments*1.0 #This is just initailizing for use later.
     
    #Not quite sure why things were done the below way, but it must work.
    solutions1 = numpy.zeros([1,len(remaining_reciprocal_correction_factors_SLS[0,:])])
    solutions = solutions1[0] #this is the variable where concentrations are placed and is the primary return.
    
    #First, remove any molecules where they have no signals due to threshold filtering etc.
    #initialize a variable for moleculeIndex before the loop across all molecules.
    moleculeIndexIncludingDeletions = 0
    for moleculeIndex in range(original_num_molecules):#array-indexed for loop. Ideally, we'll do SLS once for each molecule.
        referenceIntensitiesForThisMolecule = monitored_reference_intensities[:,moleculeIndex] 
        if sum(referenceIntensitiesForThisMolecule) == 0:
            #this is setting concentration to 0 for that molecule, in this iteration.
            solutions[moleculeIndex] = 0.0  #note that we use the actual moleculeIndex here.
            #No need to subtract any signals
            #update the used molecules list, and amounts remaining for other things.
            solvedmolecules[moleculeIndex] = 1 #note that we use the actual moleculeIndex here. 
            remaining_reciprocal_correction_factors_SLS = numpy.delete(remaining_reciprocal_correction_factors_SLS,(moleculeIndexIncludingDeletions),axis = 1)
            remaining_reference_intensities_SLS = numpy.delete(remaining_reference_intensities_SLS,(moleculeIndexIncludingDeletions),axis = 1)
            remaining_molecules_SLS = numpy.delete(remaining_molecules_SLS,(moleculeIndexIncludingDeletions))      
            moleculeIndexIncludingDeletions = moleculeIndexIncludingDeletions - 1 
            if len(uncertainties_dict) > 0:
                if type(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']) != type(None):
                    uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'] = numpy.delete(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'],(moleculeIndexIncludingDeletions),axis = 1)
        moleculeIndexIncludingDeletions = moleculeIndexIncludingDeletions + 1
    num_remaining_molecules_before_loop = len(remaining_reciprocal_correction_factors_SLS[0,:]) #could have also used a different array or way of doing this.
    
    listFor_remaining_num_molecules_during_loop = list(range(num_remaining_molecules_before_loop))
    deletionAdjustedIndex = 0
    #These nested for loops make the whole function run - they go through all the rows
    #and all the columns, but they do this as many times as there are rows, so that all
    #the values that can be found using this method will be. The values for remaining_num_MassFragments and
    #remaining_num_molecules are re-evaluted every cycle.
    for molNumIndex in listFor_remaining_num_molecules_during_loop:#array-indexed for loop. Ideally, we'll do SLS once for each molecule.     
    
        remaining_num_MassFragments = len(remaining_reciprocal_correction_factors_SLS[:,0])
        remaining_num_molecules = len(remaining_reciprocal_correction_factors_SLS[0,:])
        #FIXME: need to remove the globals here, including currentReferenceData.  these need to be passed in as arguments.
        
        if G.excludeMoleculesIfSignificantFragmentNotObserved == 'yes':
            numMoleculesExcluded = 0#reset this variable.
            remaining_num_molecules_before_excluding  = remaining_num_molecules #keep track of this variable.
            
            slsReferenceDataObject = copy.deepcopy(currentReferenceData)
            slsReferenceDataObject.monitored_reference_intensities = remaining_reference_intensities_SLS
            slsReferenceDataObject.reciprocal_matching_correction_values = remaining_reciprocal_correction_factors_SLS
            slsReferenceDataObject.molecules = remaining_molecules_SLS
            slsReferenceDataObject = signalThresholdFilter(slsReferenceDataObject, remaining_rawsignals_SLS, G.minimumSignalRequired, G.minimumStandardizedReferenceHeightToBeSignificant)
            #after done, update local variables from the object that has now been changed.
            remaining_monitored_reference_intensities = slsReferenceDataObject.monitored_reference_intensities *1.0  #These are the remaining_monitored_reference_intensities after exclusion.
            remaining_reference_intensities_SLS = slsReferenceDataObject.monitored_reference_intensities *1.0 #This one will get further shortened inside the loop.
            remaining_reciprocal_correction_factors_SLS = slsReferenceDataObject.reciprocal_matching_correction_values
            remaining_molecules_SLS = slsReferenceDataObject.molecules #should not have changed, since it just adds zeroes to patterns without removing molecules.
            #now we exclude (remove) molecules with no signals left in their reference patterns on any of the remaining mass fragments, which can happen due to the rawThresholdFilter.
            
            #Now we need to actually delete the molecules that have 0's thanks to the signalThreshold filter. The delting is what excludeEmptyMolecules does.
            remaining_num_molecules, solutions, solvedmolecules, remaining_reciprocal_correction_factors_SLS, remaining_reference_intensities_SLS, \
                                            remaining_molecules_SLS, uncertainties_dict = excludeEmptyMolecules(remaining_num_molecules, solutions,  \
                                            solvedmolecules, remaining_monitored_reference_intensities, remaining_reciprocal_correction_factors_SLS, \
                                            remaining_reference_intensities_SLS, remaining_molecules_SLS, molecules_unedited, uncertainties_dict)  
            numMoleculesExcluded = remaining_num_molecules_before_excluding - remaining_num_molecules #note that here remaining_num_molecules is after excluding.

            if numMoleculesExcluded > 0:#reset this variable.            
                #Revised the numbers since we've deleted some molecules.
                remaining_num_MassFragments = len(remaining_reciprocal_correction_factors_SLS[:,0])
                remaining_num_molecules = len(remaining_reciprocal_correction_factors_SLS[0,:])

                #Check if we should export to file what happened.
                #TODO: Below should probably be made  a function (occurs at another place below)
                specific_mass_fragment_used = 0 #Can't calculate from used_mass_fragments and last_used_mass_fragments(like in other part of code) because no specific unique fragments were actually used.
                #finding which molecules were solved by this method of excluding molecules without significant fragments.
                difference_solvedmolecules = list(solvedmolecules - last_solvedmolecules)
                #There may be more than one molecule solved, so need to make a list of them.
                specific_molecules_solved = [] #make the list blank and then will append to it.
                if 1 in difference_solvedmolecules:
                    #There could be more than 1 solved molecule, so we need to look through them one at a time to make a string.
                    for x in range(len(difference_solvedmolecules)):
                        if difference_solvedmolecules[x] == 1:
                            specific_molecules_solved.append(molecules_unedited[x])
                            G.massesUsedInSolvingMoleculesForThisPoint[x] = 0 #It is actually already initialized as zero, but this line is put here for clarity.
                
                if G.SLSUniqueExport == 'yes' and (G.answer == 'sls' or G.answer == 'autosolver'):
                    outputMoleculesAndChosenMassFragmentsFilename = 'ExportedSLSUniqueMoleculesAndChosenMassFragments.csv'
                    if G.iterativeAnalysis:
                        #then the filename will have a suffix attached
                        outputMoleculesOrderFileName = outputMoleculesAndChosenMassFragmentsFilename[:-4] + '_iter_%s' %G.iterationNumber + outputMoleculesAndChosenMassFragmentsFilename[-4:] 
                    with open(outputMoleculesAndChosenMassFragmentsFilename, 'a') as f:
                        f.write('%s,' %timeIndex)
                        f.write('%s,' %time)                    
                        f.write(str(specific_mass_fragment_used) + "," + str(list(specific_molecules_solved))[1:-1] + "\n")
                    
                    if currentReferenceData.referenceFileNameExtension == 'csv':
                        outputMoleculesOrderFileName = 'ExportedSLSUniqueMoleculesOrder.csv'
                    if currentReferenceData.referenceFileNameExtension == 'tsv':
                        outputMoleculesOrderFileName = 'ExportedSLSUniqueMoleculesOrder.tsv'
                    if G.iterativeAnalysis:
                        #then the filename will have a suffix attached
                        outputMoleculesOrderFileName = outputMoleculesOrderFileName[:-4] + '_iter_%s' %G.iterationNumber + outputMoleculesOrderFileName[-4:] 
                    with open(outputMoleculesOrderFileName, 'a') as f:
                        f.write('%s,' %timeIndex)
                        f.write('%s,' %time)
                        for x in range(len(solvedmolecules)):
                            f.write('%s,' %solvedmolecules[x])
                        f.write("SolvedMolecules \n")

                    if currentReferenceData.referenceFileNameExtension == 'csv':
                        outputMassFragmentsOrderFileName = 'ExportedSLSUniqueMassFragmentsUsed.csv'
                    if currentReferenceData.referenceFileNameExtension == 'tsv':
                        outputMassFragmentsOrderFileName = 'ExportedSLSUniqueMassFragmentsUsed.tsv'
                    if G.iterativeAnalysis:
                        #then the filename will have a suffix attached
                        outputMassFragmentsOrderFileName = outputMassFragmentsOrderFileName[:-4] + '_iter_%s' %G.iterationNumber + outputMassFragmentsOrderFileName[-4:]
                    with open(outputMassFragmentsOrderFileName, 'a') as f:
                        f.write('%s,' %timeIndex)
                        f.write('%s,' %time)                        
                        f.write(str(list(used_mass_fragments))[1:-1] + "\n") #the [1:-1] is to remove the list symbols during printing to file.
                last_solvedmolecules = solvedmolecules #some changes have occurred.
                last_used_mass_fragments = used_mass_fragments #no change has occurred.
            
            
        ####The below block of code is just to choose the next molecule to perform SLS on.###
        chosenMolecule = None
        tuplesOfUniqueFragmentsList = []
        
        if G.minimalReferenceValue == "yes": #We only will do some filtering things if it's requested.
        #Before going forward, we're going to make a variable called remaining_referenceSignificantFragmentThresholds, using a function.       
            def get_remaining_referenceSignificantFragmentThresholds(referenceSignificantFragmentThresholds, molecules_unedited, remaining_molecules_SLS):
                remaining_referenceSignificantFragmentThresholds = list(copy.deepcopy(referenceSignificantFragmentThresholds))
                molecules_unedited_to_reduce = list(copy.deepcopy(molecules_unedited))
                for moleculeIndex, moleculeName in enumerate(molecules_unedited):
                    if moleculeName in remaining_molecules_SLS:
                        pass
                    else:

                        indexToPop = molecules_unedited_to_reduce.index(moleculeName) #The list keeps shrinking, so we have to keep searching for the new/current index to pop at.
                        molecules_unedited_to_reduce.pop(indexToPop)
                        remaining_referenceSignificantFragmentThresholds.pop(indexToPop)
                return remaining_referenceSignificantFragmentThresholds
            remaining_referenceSignificantFragmentThresholds = get_remaining_referenceSignificantFragmentThresholds(G.referenceSignificantFragmentThresholds, molecules_unedited, remaining_molecules_SLS)
        
                              
        for massFragmentIndex_i in range(remaining_num_MassFragments):#array-indexed for loop (over all fragments)
            referenceIntensitiesAtThatMassFragment = remaining_reference_intensities_SLS[massFragmentIndex_i]
            ReciprocalCorrectionFactorsAtThatMassFragment = remaining_reciprocal_correction_factors_SLS[massFragmentIndex_i]
            if len(uncertainties_dict) > 0: #Just mimicing the above line to pass on the uncertainties.
                if type(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']) != type(None):
                    correctionFactorsUncertaintiesAtThatMassFragment = uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'][massFragmentIndex_i]
            signalsAtThatMassFragment = remaining_rawsignals_SLS[massFragmentIndex_i]
            if len(uncertainties_dict) > 0: #Just mimicing the above line to pass on the uncertainties.
                if (type(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS']) != type(None)) and (str(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS']).lower() != 'none'):
                    signalsAtThatMassFragmentForThisSLS_absolute_uncertainty = uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'][massFragmentIndex_i]
                else:
                    signalsAtThatMassFragmentForThisSLS_absolute_uncertainty = 0 #We will set it to zero here if there is no uncertainty provided.
            else:
                signalsAtThatMassFragmentForThisSLS_absolute_uncertainty = 0 #We will set it to zero here if there is no uncertainty provided.
            #this line checks if that mass fragment is unique to a particular molecule.
            if numpy.count_nonzero(referenceIntensitiesAtThatMassFragment) == 1:
                #the nonzero value will be the one at the maximum intensity for the reference pattern.
                valueOfUniqueStandardizedIntensity = max(referenceIntensitiesAtThatMassFragment)
                #the below nubby function will return the relevant index in array indexing.
                moleculeIndexOfUniqueIntensity = numpy.argmax(referenceIntensitiesAtThatMassFragment)
                #consideredMolecule = remaining_molecules_SLS[moleculeIndexOfUniqueIntensity] #This is for debugging purposes, this can be printed.
                #For debugging, also print remaining_reference_intensities_SLS and massFragmentIndex_i.
                #However, now we have a few lines of code to check if we are above the referenceSignificantFragmentThresholds.
                if G.minimalReferenceValue == "yes": #We only check for the remaining_referenceSignificantFragmentThresholds if this option has been chosen. 
                    if (max(remaining_reference_intensities_SLS[massFragmentIndex_i]) < remaining_referenceSignificantFragmentThresholds[moleculeIndexOfUniqueIntensity]): #This allows separate referenceSignificantFragmentThresholds for each molecule.
                        significantFragment = False  #Set to false if the fragment is too small.
                    else: 
                        significantFragment = True #This means the fragment is greater than or equal to the threshold for significance.
                if G.minimalReferenceValue != "yes": #if the option is not selected, then all fragments are considered significant.
                    significantFragment = True 
                if significantFragment == True:
                    #now make a tuple with the unique standardized intensity in the front so we can sort by that
                    ReciprocalCorrectionFactorOfUniqueIntensity = ReciprocalCorrectionFactorsAtThatMassFragment[moleculeIndexOfUniqueIntensity]
                    if len(uncertainties_dict) > 0: #Just mimicing the above line to pass on the uncertainties.
                        if type(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']) != type(None):
                            correctionFactorOfUniqueIntensity_relative_uncertainty = correctionFactorsUncertaintiesAtThatMassFragment[moleculeIndexOfUniqueIntensity]
                        else: 
                            correctionFactorOfUniqueIntensity_relative_uncertainty = 0
                    else: 
                        correctionFactorOfUniqueIntensity_relative_uncertainty = 0
                    #TODO: consider changing the primary weighting to valueOfUniqueStandardizedIntensity*signalsAtThatMassFragment
                    #and/or to a user specified argument.
                    #Note that for now, by having signals as the second slot, if two molecules each have 100% that they will sort by intensity of signals next.
                    
                    #This is for using error for primary WeightingSLS (and also prevents need to calculate it again later).
                    #now need to use the chosen mass to calculate concentration.
                    concentrationOfConsideredMolecule_relative_uncertainty = 1 #This is a default value in case the length is zero below such that it can't be calculated.
                    concentrationOfConsideredMolecule = ((float(signalsAtThatMassFragment))/float(ReciprocalCorrectionFactorOfUniqueIntensity))
                    if len(uncertainties_dict) > 0:
                            #We do not put a second "if" statement in here like for the other places, because this uncertainty will be calculated from either type (reference or collected), and if a type isn't present it is already set to 0 above.
                            #To calculate the Uncertainty:
                            #we use the following formula:  if Q = (a*b)/(c*d), Q_relative_uncertainty = sqrt( a_relative_uncertainty^2 + b_relative_uncertainty^2 + c_relative_uncertainty^2 + d_relative_uncertainty^2)
                            #Note that there is no difference between how we propagate for multiplication and division in this formula: they can be done in the same propagation step.
                            signalsAtThatMassFragmentForThisSLS_relative_uncertainty = signalsAtThatMassFragmentForThisSLS_absolute_uncertainty/signalsAtThatMassFragment #We are putting this as zero for now. TODO: Figure out if it should be non-zero instead of later step.
                            #We already have correctionFactorOfUniqueIntensity_relative_uncertaintyForThisSLS
                            concentrationOfConsideredMolecule_relative_uncertainty = (correctionFactorOfUniqueIntensity_relative_uncertainty**2 + signalsAtThatMassFragmentForThisSLS_relative_uncertainty**2)**0.5
                    #G.slsWeighting is a list-like object: [1,1,1,1] An slsWeighting value of 0 makes any of the three below factors into a coefficent of 1. A value of 1 turns it on, higher values make it even more strongly weighted.
                    #The first uses uncertainties weighting. The second solves for largest concentrations first. The third uses reference peak height. The fourth uses the signal intensity.  All can be on at the same time. 
                    moleculeWeightingSLS = (concentrationOfConsideredMolecule_relative_uncertainty**(-2*G.slsWeighting[0]))*(concentrationOfConsideredMolecule**G.slsWeighting[1])*(valueOfUniqueStandardizedIntensity**G.slsWeighting[2])*(signalsAtThatMassFragment**G.slsWeighting[3])
                    uniqueFragmentTuple = (moleculeWeightingSLS, signalsAtThatMassFragment, massFragmentIndex_i, moleculeIndexOfUniqueIntensity, ReciprocalCorrectionFactorOfUniqueIntensity, correctionFactorOfUniqueIntensity_relative_uncertainty, signalsAtThatMassFragmentForThisSLS_absolute_uncertainty, concentrationOfConsideredMolecule, concentrationOfConsideredMolecule_relative_uncertainty)
                    tuplesOfUniqueFragmentsList.append(uniqueFragmentTuple)
        #now we sort according to the biggest standardized intensities (signals as second spot), in descending order.
        tuplesOfUniqueFragmentsList.sort(reverse=True) # there is no return by list sort, the list object is directly modified.
        
        #now we simply take the first index which is the best to subtract first.
        if len(tuplesOfUniqueFragmentsList)>=1: #need at least one unique one, or can't perform SLS!
            tupleForThisSLS = tuplesOfUniqueFragmentsList[0]
            #for simplicity in reading the code, we will break out the different parts of the tuple
            signalsAtThatMassFragmentForThisSLS = tupleForThisSLS[1]
            massFragmentIndexForThisSLS = tupleForThisSLS[2]
            used_mass_fragments[massFragmentIndexForThisSLS]=1 #TODO: This should be returned so that the SLSUniqueOrder.csv can have an accompanying file of SLSUniqueOrderMassFragments
            moleculeIndexForThisSLS = tupleForThisSLS[3]
            reciprocalCorrectionFactorOfUniqueIntensityForThisSLS = tupleForThisSLS[4]
            correctionFactorOfUniqueIntensity_relative_uncertaintyForThisSLS = tupleForThisSLS[5]
            signalsAtThatMassFragmentForThisSLS_absolute_uncertainty = tupleForThisSLS[6]
            chosenMolecule = remaining_molecules_SLS[moleculeIndexForThisSLS] #This line is for debugging etc.
            chosenMolecule_original_molecular_index = list(molecules_unedited).index(chosenMolecule) #we take the chosenMolecule string, and search for it in the original list of molecules to get the original index.
            #TODO: make (or better yet, take in) a list called "moleculeSolvingOrder", append chosenMolecule to that, and return that from this function. Then we can export a file from main called moleculeSolvingOrder for each time point.
            #TODO continued: The reason to take in a list (default value blank list) is because SLSCommon may call SLSunique multiple times, so we need to append rather than just making a blank list each time.
	
            
            #now need to use the chosen mass to calculate concentration.
            concentrationOfMoleculeForThisSLS = ((float(signalsAtThatMassFragmentForThisSLS))/float(reciprocalCorrectionFactorOfUniqueIntensityForThisSLS))
            if len(uncertainties_dict) > 0:
                    #We do not put a second "if" statement in here like for the other places, because this uncertainty will be calculated from either type (reference or collected), and if a type isn't present it is already set to 0 above.
                    #To calculate the Uncertainty:
                    #we use the following formula:  if Q = (a*b)/(c*d), Q_relative_uncertainty = sqrt( a_relative_uncertainty^2 + b_relative_uncertainty^2 + c_relative_uncertainty^2 + d_relative_uncertainty^2)
                    #Note that there is no difference between how we propagate for multiplication and division in this formula: they can be done in the same propagation step.
                    signalsAtThatMassFragmentForThisSLS_relative_uncertainty = signalsAtThatMassFragmentForThisSLS_absolute_uncertainty/signalsAtThatMassFragmentForThisSLS #We are putting this as zero for now. TODO: Figure out if it should be non-zero instead of later step.
                    #We already have correctionFactorOfUniqueIntensity_relative_uncertaintyForThisSLS
                    concentrationOfMoleculeForThisSLS_relative_uncertainty = (correctionFactorOfUniqueIntensity_relative_uncertaintyForThisSLS**2 + signalsAtThatMassFragmentForThisSLS_relative_uncertainty**2)**0.5
                    uncertainties_dict['concentrations_relative_uncertainties_one_time'][chosenMolecule_original_molecular_index]=concentrationOfMoleculeForThisSLS_relative_uncertainty
                    uncertainties_dict['concentrations_absolute_uncertainties_one_time'][chosenMolecule_original_molecular_index]=abs(concentrationOfMoleculeForThisSLS_relative_uncertainty*concentrationOfMoleculeForThisSLS)
                      
            if G.slsUniquePositiveConcentrationsOnly == True:
                if concentrationOfMoleculeForThisSLS <= 0:
                    concentrationOfMoleculeForThisSLS = 0  #We can't just skip the below lines, because the molecule gets removed from arrays and info gets exported etc.
            
            ## These print statements will be preserved for debugging purposes.
            #print("Debugging","current moleculeChosen is", remaining_molecules_SLS[moleculeIndexForThisSLS], concentrationOfMoleculeForThisSLS)
            #print("Debugging","which is also", molecules[chosenMolecule_original_molecular_index], chosenMolecule_original_molecular_index)
            #print("Debugging","current moleculeChosen is", remaining_molecules_SLS[moleculeIndexForThisSLS] , concentrationOfMoleculeForThisSLS)
            #print("Debugging",remaining_molecules_SLS)
            #print("Debugging",tuplesOfUniqueFragmentsList)
            #print("Debugging","signals", remaining_rawsignals_SLS)
            #print("Debugging",remaining_reciprocal_correction_factors_SLS)
            #print("Debugging","predicted this_round signal:", concentrationOfMoleculeForThisSLS*reciprocalCorrectionFactorOfUniqueIntensityForThisSLS)
            #print("Debugging",original_list_of_mass_fragments)
            #print("Debugging",used_mass_fragments)
            
            
            ## These below lines will be preserved for debugging purposes.
            # #now we need to collect the list of masses/signals and correction factors that correspond to that molecule, i.e. moleculeIndexForThisSLS, which are nonzero.
            # for massFragmentIndex_jjj in range(remaining_num_MassFragments):
                # if remaining_reciprocal_correction_factors_SLS[massFragmentIndex_jjj,moleculeIndexForThisSLS] != 0:#If the value in the correction_values is not zero, it is kept
                    # nonZeroCorrectionValuesList.append(remaining_reciprocal_correction_factors_SLS[massFragmentIndex_jjj,moleculeIndexForThisSLS]) 
                    # signalsCorrespondingToNonZeroCorrectionValuesList.append(remaining_rawsignals_SLS[massFragmentIndex_jjj])  #This only appends signals where a correction value exists for *this* molecule, molecule_ii.
            #We need to propagate the error that is in the concentration into the simulated raw signals.
            
            
            ## These print statements will be preserved for debugging purposes.
            #print("Debugging", nonZeroCorrectionValuesList)
            #print("Debugging","here are the signals to subtract from", signalsCorrespondingToNonZeroCorrectionValuesList)
        
        
            #now we are going to solve for the signals after modifying the old code.
            #this for loop multiplies the relative amount that was just found by the correction values in order to determine the amount of signal that the 
            #certain molecule in question accounted for so that right after the for loop, this signal could be erased from the full signal, and thus the other
            #molecules relative amounts could be solved for
            #initializing an array to be populated with the amount to subtract.       
            solvedSignalsForSubtractionArray = numpy.zeros([remaining_num_MassFragments,1])                           
            solvedSignalsForSubtractionArray_relative_uncertainties = numpy.zeros([remaining_num_MassFragments,1])                           
            for massFragmentIndex_jj in range(remaining_num_MassFragments):#array-indexed for loop. #This is being called massFragmentIndex_jj to distinguish it from the outer loop.
                if remaining_reciprocal_correction_factors_SLS[massFragmentIndex_jj,moleculeIndexForThisSLS]!= 0: #used to find the raw signals that are made b the molecules that are being deleted
                    solvedSignalsForSubtractionArray[massFragmentIndex_jj] = remaining_reciprocal_correction_factors_SLS[massFragmentIndex_jj,moleculeIndexForThisSLS] * concentrationOfMoleculeForThisSLS
                    
                    if debuggingExportIndividualItem :  #This if statement is for SLS debugging purposes. 
                    
                        ObjectforExport1= str(solvedSignalsForSubtractionArray[massFragmentIndex_jj]) + " " + str(chosenMolecule) + str(moleculeIndexForThisSLS) + " " + str(original_list_of_mass_fragments[massFragmentIndexForThisSLS])
                 
                        Filename = " " + str(chosenMolecule)
                        Prefix = " solvedSignals4Sub "     #solvedSignals4Sub stands for solved signals for subtraction.
                        debuggingExportIndividualItem (Prefix, ObjectforExport1, Filename) 
                        combinedString = combinedString + Prefix + Filename + ' \n' + str(ObjectforExport1) + ' \n'

                        ObjectforExport2= str(remaining_reciprocal_correction_factors_SLS[massFragmentIndex_jj,moleculeIndexForThisSLS]) + " " + str(chosenMolecule) + str(moleculeIndexForThisSLS) + " " + str(original_list_of_mass_fragments[massFragmentIndexForThisSLS])
                 
                        Filename = " " + str(chosenMolecule)
                        Prefix = " RRCF "    #RRCF stands for Reciprocal Remaining Correction factor
                        debuggingExportIndividualItem (Prefix, ObjectforExport2, Filename) 
                        combinedString = combinedString + Prefix + Filename + ' \n' + str(ObjectforExport2) + ' \n'
                        
                        ObjectforExport3= str(concentrationOfMoleculeForThisSLS) + " " + str(chosenMolecule) + str(moleculeIndexForThisSLS) + " " + str(original_list_of_mass_fragments[massFragmentIndexForThisSLS])
                 
                        Filename = " " + str(chosenMolecule)
                        Prefix = " ConcentrationSLS "       #ConcentrationSLS stands for solved concetrations for SLS
                        debuggingExportIndividualItem (Prefix, ObjectforExport3, Filename) 
                        combinedString = combinedString + Prefix + Filename + ' \n' + str(ObjectforExport3) + ' \n'
                    
                    if len(uncertainties_dict) > 0: #Just propagating the error to make the concentration uncertainties, including the signal associated with massFragmentIndexForThisSLS.
                        solvedSignalsForSubtractionArray_relative_uncertainties[massFragmentIndex_jj] = (concentrationOfMoleculeForThisSLS_relative_uncertainty**2+(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'][massFragmentIndex_jj,moleculeIndexForThisSLS])**2)**0.5

            #The below line is the key line where the subtraction is done. 
            remaining_rawsignals_SLS = remaining_rawsignals_SLS - solvedSignalsForSubtractionArray        
            if len(uncertainties_dict) > 0: #Now need to propagating the error to make the concentration uncertainties. First need to multiply times signal to get absolute uncertainties.
                solvedSignalsForSubtractionArray_absolute_uncertainties = abs(solvedSignalsForSubtractionArray_relative_uncertainties * solvedSignalsForSubtractionArray)
                #now need to take the remaining_rawsignals_SLS uncertainties, and combine them after zeroing out the signal associated with massFragmentIndexForThisSLS.
                uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'][massFragmentIndexForThisSLS] = 0.0
                #To combine the uncertainties we do the   c_absolute_uncertainty = sqrt(a_absolute_uncertainty^2 + a_absolute_uncertainty^2).  For arrays, need to use numpy.square to square each element in the array. For arbitrary powers, we use numpy.power, which is slower and less accurate.
                sqrt_term = numpy.square(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS']) + numpy.square(solvedSignalsForSubtractionArray_absolute_uncertainties)
                #uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'] = numpy.power(sqrt_term,0.5) #apparently numpy.power is slow, so should use **.
                uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'] = sqrt_term**0.5
                
            #Since it's done, we'll update the solved molecules array etc.
            solutions[chosenMolecule_original_molecular_index] = concentrationOfMoleculeForThisSLS
                        
            if debuggingExportIndividualItem :   #This if statement is for SLS debugging purposes.
            
                ObjectforExport4= str(solutions[chosenMolecule_original_molecular_index]) + " " + str(chosenMolecule) + " " + str(moleculeIndexForThisSLS) 
             
                Filename = " " + str(chosenMolecule)
                Prefix = " Solutions "
                debuggingExportIndividualItem (Prefix, ObjectforExport4, Filename) 
                combinedString = combinedString + Prefix + Filename + ' \n' + str(ObjectforExport4) + ' \n'

                debuggingExportIndividualItem ( " SLS Step By Step ", combinedString) #Exports SLS in step by step format for debugging purposes. 
            
            solvedmolecules[chosenMolecule_original_molecular_index] = 1 #This updates a list that keeps track of which molecules have been used up.

       
            #now delete that molecule from the correction values array, etc.
            remaining_reciprocal_correction_factors_SLS = numpy.delete(remaining_reciprocal_correction_factors_SLS,(moleculeIndexForThisSLS),axis = 1)
            remaining_reference_intensities_SLS = numpy.delete(remaining_reference_intensities_SLS,(moleculeIndexForThisSLS),axis = 1)
            remaining_molecules_SLS = numpy.delete(remaining_molecules_SLS,moleculeIndexForThisSLS)
            #update this variable:
            remaining_num_molecules = remaining_num_molecules -1
            if len(uncertainties_dict) > 0:
                if type(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']) != type(None):
                    uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'] = numpy.delete(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'],(moleculeIndexForThisSLS),axis = 1)
                    uncertainties_dict['correction_values_relative_uncertainties_one_time'] = uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']*1.0
                    uncertainties_dict['reciprocal_matching_correction_values_relative_uncertainties_one_time'] = uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']*1.0
            #Reset these lists to start again.
            nonZeroCorrectionValuesList = []
            signalsCorrespondingToNonZeroCorrectionValuesList = []
        
            # This block of code is a printing statement to show the user what order the molecules are being solved in
            # This is a csv file so should be delimited with commas
            #TODO: Below should probably be made  a function (occurs at another place above)
            
            
            #Find the most recently used mass fragment and solved molecule.
            difference_used_mass_fragments = list((used_mass_fragments - last_used_mass_fragments)) #We subtract the two arrays looking for a "1" for which mass fragment was used. But there are special cases with no mass fragment used.
            if 1 in difference_used_mass_fragments:
                index_of_mass_fragment = difference_used_mass_fragments.index(1) 
                specific_mass_fragment_used = original_list_of_mass_fragments[index_of_mass_fragment] #Find subtract the arrays, make them a list, use index function to find the position of the new "1", call the mass fragment.
            else:
                specific_mass_fragment_used = 0 
                
            difference_solvedmolecules = list((solvedmolecules - last_solvedmolecules))
            if 1 in difference_used_mass_fragments:
                index_of_molecule = difference_solvedmolecules.index(1)
                specific_molecule_solved = molecules_unedited[ index_of_molecule ] #Find subtract the arrays, make them a list, use index function to find the position of the new "1", call the molecule.
                G.massesUsedInSolvingMoleculesForThisPoint[index_of_molecule] = specific_mass_fragment_used
            else:
                specific_molecule_solved = 'None'
            
            if G.SLSUniqueExport == 'yes':
                if currentReferenceData.referenceFileNameExtension == 'csv':
                    outputMoleculesAndChosenMassFragmentsFilename = 'ExportedSLSUniqueMoleculesAndChosenMassFragments.csv'
                if currentReferenceData.referenceFileNameExtension == 'tsv':
                    outputMoleculesAndChosenMassFragmentsFilename = 'ExportedSLSUniqueMoleculesAndChosenMassFragments.tsv'
                if G.iterativeAnalysis:
                    #then the filename will have a suffix attached
                    outputMoleculesOrderFileName = outputMoleculesAndChosenMassFragmentsFilename[:-4] + '_iter_%s' %G.iterationNumber + outputMoleculesAndChosenMassFragmentsFilename[-4:] 
                with open(outputMoleculesAndChosenMassFragmentsFilename, 'a') as f:
                    f.write('%s,' %timeIndex)
                    f.write('%s,' %time)                    
                    f.write(str(specific_mass_fragment_used) + "," + specific_molecule_solved + "\n")

                if currentReferenceData.referenceFileNameExtension == 'csv':
                    outputMoleculesOrderFileName = 'ExportedSLSUniqueMoleculesOrder.csv'
                if currentReferenceData.referenceFileNameExtension == 'tsv':
                    outputMoleculesOrderFileName = 'ExportedSLSUniqueMoleculesOrder.tsv'
                if G.iterativeAnalysis:
                    #then the filename will have a suffix attached
                    outputMoleculesOrderFileName = outputMoleculesOrderFileName[:-4] + '_iter_%s' %G.iterationNumber + outputMoleculesOrderFileName[-4:] 
                with open(outputMoleculesOrderFileName, 'a') as f:
                    f.write('%s,' %timeIndex)
                    f.write('%s,' %time)
                    for x in range(len(solvedmolecules)):
                        f.write('%s,' %solvedmolecules[x])
                    f.write("SolvedMolecules \n")
                
                outputMassFragmentsOrderFileName = 'ExportedSLSUniqueMassFragmentsUsed.csv'
                if G.iterativeAnalysis:
                    #then the filename will have a suffix attached
                    outputMassFragmentsOrderFileName = outputMassFragmentsOrderFileName[:-4] + '_iter_%s' %G.iterationNumber + outputMassFragmentsOrderFileName[-4:] 
                with open(outputMassFragmentsOrderFileName, 'a') as f:
                    f.write('%s,' %timeIndex)
                    f.write('%s,' %time)                    
                    f.write(str(list(used_mass_fragments))[1:-1] + "\n") #the [1:-1] is to remove the list symbols during printing to file.
            last_solvedmolecules = solvedmolecules*1.0 #Be careful not to simply make a pointer!
            last_used_mass_fragments = used_mass_fragments*1.0 #Be careful not to simply make a pointer!
        if remaining_num_molecules == 0:
            #TODO: Add in an if G.SLSUniqueExport == 'yes' here that exports the full line of solved molecules and masses at one time.
            break
        
    if remaining_reciprocal_correction_factors_SLS.size > 0:#if there are correction values left (i.e. not all the solutions have been found)
        #this for loop is used to delete any rows entirely composed of zeros, since the molecules percentages are found
        #based on the column index, the rest of the loop still works fine, this is just so that there are not singular
        #matrices made by combination maker. (originally there would be up to hundreds of simply useless matrices solved)
        place_holder = 0
        for correctionIndex in range(len(remaining_reciprocal_correction_factors_SLS[:,0])):#array-indexed for loop
            if any(remaining_reciprocal_correction_factors_SLS[correctionIndex-place_holder,:]) == 0:#all rows of only zeros are deleted
                remaining_reciprocal_correction_factors_SLS = numpy.delete(remaining_reciprocal_correction_factors_SLS,correctionIndex-place_holder,axis = 0)
                remaining_rawsignals_SLS = numpy.delete(remaining_rawsignals_SLS,correctionIndex-place_holder,axis = 0)
                remaining_reference_intensities_SLS = numpy.delete(remaining_reference_intensities_SLS,correctionIndex-place_holder,axis = 0)
                if len(uncertainties_dict) > 0: #Will just mimic the deletions done for the above lines in the uncertainties.
                    if (type(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']) != type(None)) and (str(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']).lower() != 'none'):
                        uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'] = numpy.delete(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'],correctionIndex-place_holder,axis = 0)
                    if (type(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS']) != type(None)) and (str(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS']).lower() != 'none'):
                        uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'] = numpy.delete(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'],correctionIndex-place_holder,axis = 0)
                place_holder = place_holder + 1#since the arrays are being deleted, this keeps the indexing correct
    if sum(solvedmolecules) == 0:#if none of the solutions have been found
        solutions = []
        solvedmolecules = []   
    return [remaining_molecules_SLS,remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,solutions,molecules_unedited,solvedmolecules, uncertainties_dict]
    
#this sls method cuts smaller, solvable arrays out of the large array and uses numpy.linalg.solve to find the signals
#relative to Co for those molecules, then via the SLSMethod function below the function works in conjunction with the 
#sls unique fragments in order to solve all the molecules in the array that can be solved via arrays of ize one, two or
#three (ultimately, any size array can be solved via this method (if the method went to infinity rather than stopping at
#3) provided that the array were square- for, in fact, none of these method (save for the brute method) can solve for the 
#signals if the first array is not square
def SLSCommonFragments(reciprocal_matching_correction_values,rawsignalsarrayline,monitored_reference_intensities,molecules,scaledConcentrationsarray,molecules_unedited,conversionfactor,datafromcsv,DataRangeSpecifierlist,objectiveFunctionType,counterforspecifications,maxPermutations = 100001, uncertainties_dict={}):

    #TODO: #FIXME: It seems 
    # like this function forces Brute towards the end, rather than checking 
    # what the user defined finisher is. This should be checked, eventually.

    # This is creating a local copy of the monitored reference intensities which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reference_intensities_SLS = copy.deepcopy(monitored_reference_intensities)

    # This is creating a local copy of the reciprocal_matching_correction_values which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reciprocal_correction_factors_SLS = copy.deepcopy(reciprocal_matching_correction_values)

    # This is creating a local copy of the rawsignalsarrayline which will become
    # truncated as the molecules are solved and masses are removed
    remaining_rawsignals_SLS = copy.deepcopy(rawsignalsarrayline)

    # This is creating a local copy of 'molecules' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_molecules_SLS = copy.deepcopy(molecules)
    
    sizedecide = 1 #right now, can cause crashing when given values greater than 1. This has to do with how much depth the nesting can have in SLSCommonFragments.
    reminder3 = 0
    #molecules_unedited = remaining_molecules_SLS
    selectedVariable2DArrayAsListholder = []
    selectedcorrectionsholder = []
    selectedmassfragsholder = []
    selectedrawsigsholder = []
    firstTimeCounter = 1
    #this for loop goes through the whole nested loop as many times as the user would like, with the default setting being 3,
    #It resets the selected variable holder, which finds each row with the same number of variables 
    for sizecounter in range(2,sizedecide+1):#array-indexed for loop
    #if the array is empty here then the values are set to zero using the else statement
    #the original indexeers will cause an error because they require a two dimensional array
        if remaining_reciprocal_correction_factors_SLS.shape[0] != 0:
            remaining_num_MassFragments = len(remaining_reciprocal_correction_factors_SLS[:,0])
            remaining_num_molecules = len(remaining_reciprocal_correction_factors_SLS[0,:])
        else:  #####TODO: Figure out why this else statement is here.
            remaining_num_MassFragments = len(remaining_reciprocal_correction_factors_SLS)
            remaining_num_molecules = len(remaining_reciprocal_correction_factors_SLS)
        selectedVariable2DArrayAsList = []
        selectedcorrections = []
        selectedmassfrags = []
        selectedrawsigs = []
        useablecorrections = []
        useablemassfrags = []
        useablerawsigs = []
        #this loop looks at each row for the correction value array, resets the selected variable value
        for rowcounter in range(remaining_num_MassFragments):#array-indexed for loop
            selectedvariables = []
            reminder2 = 0
            reminder4 = 0
            #this loop looks at each molecule's correction value for each row, only if the firstTimeCounter is equal to zero(will be
            #used later)
            for molecule_ii in range(remaining_num_molecules):#array-indexed for loop
                if firstTimeCounter == 1:#will only run if no row has been chosen already
                    #the code below this only happens if the firstTimeCounter is equal to one
                    #here a zero or a one is added to the selected variable, based on whether the correction value array has a
                    #value present or not, this way the different rows can be compared
                    if remaining_reciprocal_correction_factors_SLS[rowcounter,molecule_ii] != 0:#chooses all values that are non zero- inputs ones
                        selectedvariables.append(1)  
                    else:
                        selectedvariables.append(0) 
                    #once the selected variable has been completed, this if statement is called
                    if molecule_ii == remaining_num_molecules-1: #at the end of the row (last molecule to consider) based on array indexing
                        #this if statement gets all the rows with the needed 
                        if sum(selectedvariables) == sizecounter: #if the row has the correct number of variables
                            selectedVariable2DArrayAsList.append(selectedvariables)
                            selectedcorrections.append(remaining_reciprocal_correction_factors_SLS[rowcounter])
                            selectedmassfrags.append(remaining_reference_intensities_SLS[rowcounter])
                            selectedrawsigs.append(remaining_rawsignals_SLS[rowcounter])
            #If the loop has finished making the selected variables list, then this second loop starts, and checks each of these
            #lists to see if this variable has enough rows with the same columns being used to be used as our array to solve
            if rowcounter == remaining_num_MassFragments - 1:#after all the selected variables have been put in arrays
                #This is being called massFragmentIndex_jj to distinguish it from the outer loop.
                for massFragmentIndex_jj in range(remaining_num_MassFragments):#array-indexed for loop
                    currentSelectedVariablesOption = []
                    reminder2 = 0
                    for molecule_iii in range(remaining_num_molecules):#array-indexed for loop
                        if firstTimeCounter == 1:#will only run if no row has been chosen already
                            if remaining_reciprocal_correction_factors_SLS[massFragmentIndex_jj,molecule_iii] != 0:#makes the selected variable again
                                currentSelectedVariablesOption.append(1)  
                            else:
                                currentSelectedVariablesOption.append(0)
                            if molecule_iii == remaining_num_molecules-1: #once it has been made
                            #this for loop iterates across all of the values that are the same length as the variable we're 
                            #checking for, it then makes an array using all of the rows that match (it will at least match with
                            #itself, because itself is in the selected variables)
                                for selectedIndex in range(len(selectedVariable2DArrayAsList)):#array-indexed for loop. This loops as many times as there are COMBINATIONS from the first loop.
                                    if selectedVariable2DArrayAsList[selectedIndex] == currentSelectedVariablesOption:#adds itself to the list of matching rows,checks other rows of same size
                                        reminder2 = reminder2 + 1   #This increments the counter named reminder2 when a match is found.
                                        if reminder2 == 1:# If it's the first match, it keeps the values that will be used in solving; if this row ends up being used
                                            useablecorrections = numpy.array(selectedcorrections[selectedIndex])
                                            useablemassfrags = numpy.array(selectedmassfrags[selectedIndex])
                                            useablerawsigs = numpy.array(selectedrawsigs[selectedIndex])
                                        else:  #If it's not the first match, it uses a vstack to add to the existing matches. A match is "useable"
                                            useablecorrections = numpy.vstack([useablecorrections,selectedcorrections[selectedIndex]])
                                            useablemassfrags = numpy.vstack([useablemassfrags,selectedmassfrags[selectedIndex]])
                                            useablerawsigs = numpy.vstack([useablerawsigs,selectedrawsigs[selectedIndex]])
                                #if you are checking for arrays larger than size two then the size three array will check all the size 
                                #two arrays for lines that only contain values in the columns that have values for this size 3 row; 
                                #adding the matching rows to an array, the same one as before
                                if reminder4 == 0 and sizecounter > 2:#checks against past sizes for matching within one (has only values where the chosen row has values),reminder4
                                                                      #says that this has already been done so that the function will not add a bunch of extra rows and therefore
                                                                      #make the rest of the array unsolvable
                                    for rowIndex in range(len(selectedVariable2DArrayAsListholder)):#this goes through each layer of chosen rows
                                        for colIndex in range(len(selectedVariable2DArrayAsListholder[rowIndex])):#this goes through the rows of each chosen layer previous to the one being looked at
                                            #the discrepancy between the rows 
                                            if (len(currentSelectedVariablesOption) - sum(numpy.array(currentSelectedVariablesOption) == numpy.array(selectedVariable2DArrayAsListholder[rowIndex][colIndex]))) <= sizecounter - (rowIndex + 2):#if they are good, then it adds all the required fields from that row
                                                reminder2 = reminder2 + 1
                                                reminder4 = reminder4 + 1
                                                if reminder2 == 1: #if this is the first row then everything else is just stacked on it
                                                    useablecorrections = selectedcorrectionsholder[rowIndex][colIndex]
                                                    useablemassfrags = selectedmassfragsholder[rowIndex][colIndex]
                                                    useablerawsigs = selectedrawsigsholder[rowIndex][colIndex]
                                                else:#just stacked on first row
                                                    useablecorrections = numpy.vstack([useablecorrections,selectedcorrectionsholder[rowIndex][colIndex]])
                                                    useablemassfrags = numpy.vstack([useablemassfrags,selectedmassfragsholder[rowIndex][colIndex]])
                                                    useablerawsigs = numpy.vstack([useablerawsigs,selectedrawsigsholder[rowIndex][colIndex]])
                                #if there are a number of rows equal to or greater than the size that is being looked for then this set of
                                #rows can be used, and are used, in order to solve for a certain set of molecules (either two or three),
                                #then it deletes these values from the total correction value array and deletes the columns not being used
                                #in the synthesized arrays and solves for their signals relative to CO
                                if reminder2 >= sizecounter:#If there are enough rows
                                    place_holder = 0
                                    for deleter in range(len(remaining_reciprocal_correction_factors_SLS[:,0])):#array-indexed for loop
                                        place_holder2 = 0
                                        for checker in range(len(useablecorrections[:,0])):#array-indexed for loop
                                            if place_holder2 == 0:
                                                if all(remaining_reciprocal_correction_factors_SLS[deleter-place_holder,:] == useablecorrections[checker,:]):#gets index of correction values to delete- the ones chosen
                                                    remaining_reciprocal_correction_factors_SLS = numpy.delete(remaining_reciprocal_correction_factors_SLS,deleter-place_holder,axis = 0)
                                                    remaining_rawsignals_SLS = numpy.delete(remaining_rawsignals_SLS,deleter-place_holder,axis = 0)
                                                    remaining_reference_intensities_SLS = numpy.delete(remaining_reference_intensities_SLS,deleter-place_holder,axis = 0)
                                                    if len(uncertainties_dict) > 0: #will just mimic the above lines for the uncertainties.
                                                        if type(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']) != type(None):
                                                            uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'] = numpy.delete(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'],deleter-place_holder,axis = 0)
                                                        if type(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS']) != type(None):
                                                            uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'] = numpy.delete(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'],deleter-place_holder,axis = 0)
                                                    place_holder = place_holder + 1
                                                    place_holder2 = place_holder2 + 1
                                    place_holder = 0
                                    place_holder2 = 0
                                    solvedmolecules = currentSelectedVariablesOption #Is this just to get the array size? Should a copy command be used here?
                                    for selectedcounter in range(len(currentSelectedVariablesOption)):#array-indexed for loop
                                        if currentSelectedVariablesOption[selectedcounter - place_holder] == 0:#all of the zeros in the chosen rows are eliminated, because they will only make solving more difficult
                                            currentSelectedVariablesOption = numpy.delete(currentSelectedVariablesOption,selectedcounter - place_holder)
                                            useablecorrections = numpy.delete(useablecorrections,selectedcounter - place_holder,axis = 1)
                                            useablemassfrags = numpy.delete(useablemassfrags,selectedcounter - place_holder,axis = 1)
                                            place_holder = place_holder + 1 #saves place when deleting
                                        if solvedmolecules[selectedcounter] == 1:# in same loop, the values that are going to be used will be deleted from the arrays getting passed on- that are not solved yet
                                            remaining_reciprocal_correction_factors_SLS = numpy.delete(remaining_reciprocal_correction_factors_SLS,selectedcounter - place_holder2,axis = 1)
                                            remaining_reference_intensities_SLS = numpy.delete(remaining_reference_intensities_SLS,selectedcounter - place_holder2,axis = 1)
                                            remaining_molecules_SLS = numpy.delete(remaining_molecules_SLS,selectedcounter - place_holder2)
                                            place_holder2 = place_holder2 + 1 #saves place when deleting
                                    useablemassfrags,useablecorrections,useablerawsigs = DistinguishedArrayChooser (useablemassfrags,useablecorrections,useablerawsigs, G.moleculeLikelihoods, G.sensitivityValues)
                                    if numpy.linalg.det(useablecorrections) != 0: #solves if det is not zero
                                        #the counter below is equal to zero so that the molecule chooser will just choose the solutions given in the last line
                                        #the molecules, correction and raw signals arrays here are fragments of the whole, using the function in the same way that
                                        #it is used in the SLSMethod function
                                        place_holder3 = 0
                                        useablemolecules = []
                                        for molecules_unedited_counter in range(len(molecules_unedited)):
                                            for moleculescounter in range(len(remaining_molecules_SLS)):
                                                if molecules_unedited[molecules_unedited_counter] == remaining_molecules_SLS[moleculescounter]:
                                                    place_holder3 = place_holder3 + 1
                                                if moleculescounter == len(remaining_molecules_SLS) - 1:
                                                    if place_holder3 == 0:
                                                        useablemolecules.append(molecules_unedited[molecules_unedited_counter])
                                            place_holder3 = 0
                                        if len(useablemolecules) == 0:
                                            useablemolecules = molecules_unedited
                                        #this works in the same way as the brute force method for the remaining matrix at the end of the sls method, where when
                                        #the counter is equal to zero the moleculespecifier chooses all of the array line. This is because when the counter is zero
                                        #the input is always the predicted values according to the inverse method, while afterwards it is the actual signals array
                                        #in this way, for when the counter is zero, the program changes the value of the array to the predicted solutions
                                        if counterforspecifications == 0:
                                            scaledConcentrationsarray = numpy.linalg.solve(useablecorrections,useablerawsigs)
                                        specifications = DataRangeSpecifier(useablemolecules,counterforspecifications,molecules_unedited,conversionfactor,datafromcsv,DataRangeSpecifierlist,scaledConcentrationsarray)
                                        solutions = OptimizingFinisher(useablemolecules,specifications,useablecorrections,useablerawsigs,objectiveFunctionType,maxPermutations)
                                    else:
                                        print('The Array Chosen is Singular')
                                        solutions = numpy.zeros(len(remaining_rawsignals_SLS)) # the solutions are made into all zeros if the chosen array is singular
                                    firstTimeCounter = 0
                                    reminder3 = 1
        selectedVariable2DArrayAsListholder.append(selectedVariable2DArrayAsList)
        selectedcorrectionsholder.append(selectedcorrections)
        selectedmassfragsholder.append(selectedmassfrags)
        selectedrawsigsholder.append(selectedrawsigs)
    if reminder3 == 0:#if no solutions have been found
        solutions = []
        solvedmolecules = []
    else: #sets solutions based on the used molecules array that was made
        solutionsholder = numpy.zeros(len(solvedmolecules))
        place_holder = 0 
        
        for solvedmoleculesIndex in range(len(solvedmolecules)):#array-indexed for loop

            if solvedmolecules[solvedmoleculesIndex] == 1:#when there is a molecule in that position
                solutionsholder[solvedmoleculesIndex] = solutions[place_holder]
                place_holder = place_holder + 1 #helps add the two arrays together
        solutions = solutionsholder
    return [molecules_unedited,remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,solutions,remaining_molecules_SLS,solvedmolecules, uncertainties_dict]


#this function simply calls the other functions to be used, based on the user input pathway, that means that this
#function can send the sls to unique or common fragments, to inverse or brute method after, and sends the data back 
#and forth between the unique and common fragments for the common fragments method
def SLSMethod(molecules,monitored_reference_intensities,reciprocal_matching_correction_values,rawsignalsarrayline,timeIndex,conversionfactor,datafromcsv,molecules_copy,DataRangeSpecifierlist,SLSChoices,mass_fragment_numbers,permutationNum,scaledConcentrationsarray,objectiveFunctionType,time,maxPermutations=100001, bestMassFragChooser=False, uncertainties_dict={}):
    #FIXME: uncertainties_dict is supposed to be passed in as an argument, but for now I am taking it from othe globals because I had some strange empty list error for the solutions concentrations when I tried passing it in.
    uncertainties_dict = G.uncertainties_dict

    # This is creating a local copy of the monitored_reference_intensities which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reference_intensities_SLS = copy.deepcopy(monitored_reference_intensities)

    # This is creating a local copy of the reciprocal_matching_correction_values which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reciprocal_correction_factors_SLS = copy.deepcopy(reciprocal_matching_correction_values)

    # This is creating a local copy of the rawsignalsarrayline which will become
    # truncated as the molecules are solved and masses are removed
    remaining_rawsignals_SLS = copy.deepcopy(rawsignalsarrayline)

    # This is creating a local copy of 'molecules' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_molecules_SLS = copy.deepcopy(molecules)
    
    # Read sls options from input files
    [uniqueOrCommon,slsFinish,distinguished] = SLSChoices
                   
    if uniqueOrCommon == 'unique': #user input
        [remaining_molecules_SLS,remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,solutions,molecules_unedited,solvedmolecules, uncertainties_dict] = SLSUniqueFragments(remaining_molecules_SLS,remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS, timeIndex, time, uncertainties_dict)
    #should the user choose the common fragments method, it must be realized that once you solve a two by two or three
    #by three array using the common fragments method you can then possibly solve using the unique fragments method once 
    #again. This while loop makes sure that as long as there are possibilities after either method is used, it will keep
    #running until it has to use the inverse/brute method to finish
    elif uniqueOrCommon == 'common':
        place_holder = 1
        unique = []
        common = []
        order = []
        while place_holder != 0: # until a whole loop happens with no changes made, the while loop continues to run
            place_holder = 0
            for desirednum in range(1,4):#checks all sizes 1,2,3
                zerosrows = 0
                if remaining_reciprocal_correction_factors_SLS.size != 0: #if the array exists
                    for rowcounter in range(len(remaining_reciprocal_correction_factors_SLS[:,0])):#array-indexed for loop
                        zeros = 0
                        for columncounter in range(len(remaining_reciprocal_correction_factors_SLS[0,:])):#array-indexed for loop
                            if remaining_reciprocal_correction_factors_SLS[rowcounter,columncounter] == 0:#if there is a zero in the array
                                zeros = zeros + 1
                            if columncounter == len(remaining_reciprocal_correction_factors_SLS[0,:])-1:#if the whole row has been checked
                                if zeros == len(remaining_reciprocal_correction_factors_SLS[0,:])-(desirednum):#if the number of ones if the row is equal to the size being checked for
                                    zerosrows = zerosrows + 1
                    if zerosrows >= desirednum: #if the number of matching rows is equal to or greater than the size being looked for
                        if desirednum == 1: #unique method used if there are rows with only one correction value
                            [remaining_molecules_SLS,remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,solutions,molecules_unedited,solvedmolecules, uncertainties_dict] = SLSUniqueFragments (remaining_molecules_SLS,remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS, timeIndex, time, uncertainties_dict)
                            if len(solvedmolecules) != 0:
                                unique.append([solutions,molecules_unedited,solvedmolecules])
                                order.append('unique')
                            place_holder = place_holder +1 #if this happens, than place_holder is made equal to 1, so that the loop must run all the way through again
                            if len(solutions) == 0: #in case no solutions were found
                                place_holder = place_holder - 1 #while loop can stop- so it does not get stuck
                        else: #common fragments method 
                            
                            [molecules_unedited,remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,solutions,remaining_molecules_SLS,solvedmolecules, uncertainties_dict] = SLSCommonFragments (remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_SLS,remaining_molecules_SLS,scaledConcentrationsarray,molecules_copy,conversionfactor,datafromcsv,DataRangeSpecifierlist,objectiveFunctionType,timeIndex, maxPermutations, uncertainties_dict)
                            if len(solvedmolecules) != 0:
                                common.append([solutions,molecules_unedited,solvedmolecules])
                                order.append('common')
                            place_holder = place_holder + 1#if this happens, than place_holder is made equal to 1, so that the loop must run all the way through again
                            if len(solutions) == 0:#in case no solutions found 
                                place_holder = place_holder - 1#while loop can stop- so it does not get stuck
        #this section of if statements nested inside a for loop puts all the data that was just gathered together into one set of
        #solutions and molecules names, that way, this can be added to the solutions in the threshold function, and have values that
        #are obtained from the brute/inverse method stored here as well. 
        unique_place_holder = 0
        common_place_holder = 0    
        ulength = len(unique)
        clength = len(common)
        thenext = []
        for ordercounter in range(len(order)-1,-1,-1):#array-indexed for loop, goes backwards because we want to start at the end of the unique and common lists- whatever was added last
            if ordercounter == len(order)-1:#The first index
                if order[ordercounter] == 'unique': #if the last function was unique then the arrays from that are put in the 'next' list 
                    unique_place_holder = unique_place_holder + 1#unique has been used
                    thenext.append(unique[ulength-1])
                if order[ordercounter] == 'common':#if the last function was common then the arrays from that are put in the 'next' list 
                    common_place_holder = common_place_holder + 1 #common has been used 
                    thenext.append(common[clength-1])
            else:#array-indexed for loop
                if order[ordercounter] == 'unique':#This is how the information from after this (in next) is appended to the unique soltuions
                    adderhelp = 0
                    place_holder = 0
                    #unique_place_holder is the number of times unique has been called from, 1 is added to it because the last index of unique is ulength -1
                    for adder in range(len(unique[ulength-(1+unique_place_holder)][2])):#array-indexed for loop, each list item is a list, holding solutions in [0],stable molecules in [1] and used molecules in [2], only [0] and [2] are concatenated each time 
                        if unique[ulength-(1+unique_place_holder)][2][adder] == 0:#adds the one below where there are zeros in the used molecules
                            unique[ulength-(1+unique_place_holder)][2][adder] = unique[ulength-(1+unique_place_holder)][2][adder] +thenext[0][2][adderhelp]
                            if thenext[0][2][adderhelp] == 1:#if there is a solution in the one below, it will add the solutions together
                                unique[ulength-(1+unique_place_holder)][0][adder] = unique[ulength-(1+unique_place_holder)][0][adder] +thenext[0][0][adderhelp-place_holder]
                            else:
                                place_holder = place_holder + 1
                            adderhelp = adderhelp + 1
                    thenext.append(unique[ulength-(1+unique_place_holder)]) #this is saved for the next loop
                    unique_place_holder = unique_place_holder + 1 #this helps know how many times the unique method has had one of its solutions called (of the list)
                if order[ordercounter] == 'common': #all of the common additions are done the same way as the unique b/c they have the same names and formats
                    adderhelp = 0
                    place_holder = 0
                    for adder in range(len(common[clength-(1+common_place_holder)][2])):
                        if common[clength-(1+common_place_holder)][2][adder] == 0:
                            common[clength-(1+common_place_holder)][2][adder] = common[clength-(1+common_place_holder)][2][adder] +thenext[0][2][adderhelp]
                            if thenext[0][2][adderhelp] == 1:
                                common[clength-(1+common_place_holder)][0][adder] = common[clength-(1+common_place_holder)][0][adder] +thenext[0][0][adderhelp-place_holder]
                            else:
                                place_holder = place_holder + 1
                            adderhelp = adderhelp + 1
                    thenext.append(common[clength-(1+common_place_holder)])
                    common_place_holder = common_place_holder + 1
                thenext.pop(0) #since there was another item added to the list, this deleted the item just used
            if ordercounter == 0: #at the end of the loop
                if order[ordercounter] == 'unique':#if the last one of the order was unique
                    solutions = unique[0][0] #assigns the solutions
                    molecules_unedited = unique[0][1] #assigns the molecules from the last order (these are the stable molecules for this whole method)
                    solvedmolecules = unique[0][2] #solvedmolecules is set based on the first one used as well
                if order[ordercounter] == 'common':#if the last item of the order list was common
                    solutions = common[0][0]
                    molecules_unedited = common[0][1]
                    solvedmolecules = common[0][2]
    else: # throw error because uniqueOrCommon is not set properly
        raise ValueError("The value of 'uniqueOrCommon' is {}, it should be either 'common' or 'unique'".format(uniqueOrCommon)) 
    
    if bestMassFragChooser: return remaining_molecules_SLS
    
    #if the sls method does not solve for all the molecules, then the rest are sent to the inverse method 
    #where the remaining matrix is solved
    if remaining_reciprocal_correction_factors_SLS.size != 0:#if everything hasn't already been solved
        #If there are uncertainties, we need to make a copy of the existing uncertainties and then add in the remaining solved molecules at the end.
        if (G.calculateUncertaintiesInConcentrations == True):
                uncertainties_dict['concentrations_relative_uncertainties_one_time_before_finisher'] = uncertainties_dict['concentrations_relative_uncertainties_one_time']*1.0
                uncertainties_dict['concentrations_absolute_uncertainties_one_time_before_finisher'] = uncertainties_dict['concentrations_absolute_uncertainties_one_time']*1.0
        if slsFinish == 'inverse':#if inverse finish is chosen
             if distinguished == 'yes':#user input, choosing between distinguished inverse method or combinations method
                 concentrationsFromFinisher = InverseMethodDistinguished (remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS, uncertainties_dict)
             else:
                 concentrationsFromFinisher = InverseMethod (remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_SLS,mass_fragment_numbers,remaining_molecules_SLS,'composition')
           
        if slsFinish == 'brute':#if brute method is chosen
            if timeIndex == 0:#the first time is always run through the inverse method, where the ranges can use this information the loops afterwards
                if distinguished == 'yes':#user input, choosing between distinguished inverse method or combinations method
                    concentrationsFromFinisher = InverseMethodDistinguished (remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS, uncertainties_dict)
                else:
                     concentrationsFromFinisher = InverseMethod (remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_SLS,mass_fragment_numbers,remaining_molecules_SLS,'composition')
                if math.log(permutationNum,5) >= len(remaining_molecules_SLS):
                    #Ashi believes the above comparision is to ensure each molecule has at least 5 concentrations checked
                    specifications = DataRangeSpecifier(remaining_molecules_SLS,timeIndex,molecules_copy,conversionfactor,datafromcsv,DataRangeSpecifierlist,concentrationsFromFinisher)
                    concentrationsFromFinisher = OptimizingFinisher(remaining_molecules_SLS,specifications,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,objectiveFunctionType,maxPermutations)
                else:
                    print("Warning: The number of permutations requested is too small to allow for 5 possibilities per molecule. "
                          + "Switching to use Inverse instead of Brute for slsFinish for this datapont (at timeIndex of " +str(timeIndex)+ " where 0 is the first analyzed datapoint)."
                          + "Additional Info: There are " + str(len(remaining_molecules_SLS)) + " molecules that SLS was unable to solve! "
                          + "There are " + str(permutationNum) + " permutations allowed, and " + str(5**len(remaining_molecules_SLS)) + " would be needed."
                          + "If you wish to use Brute, increase the size of permutations in the user input file. ")
                    if distinguished == 'yes':#user input, choosing between distinguished inverse method or combinations method
                        concentrationsFromFinisher = InverseMethodDistinguished (remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS, uncertainties_dict)
                    else:
                        concentrationsFromFinisher = InverseMethod(remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_SLS,mass_fragment_numbers,remaining_molecules_SLS,'composition')               

            else: #after the first time these functions are called
                if math.log(permutationNum,5) >= len(remaining_molecules_SLS):
                    #Ashi believes the above comparision is to ensure each molecule has at least 5 concentrations checked
                    specifications = DataRangeSpecifier(remaining_molecules_SLS,timeIndex,molecules_copy,conversionfactor,datafromcsv,DataRangeSpecifierlist,
                                                        scaledConcentrationsarray)
                    concentrationsFromFinisher = OptimizingFinisher(remaining_molecules_SLS,specifications,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,objectiveFunctionType, maxPermutations)
                else:
                    print("Warning: The number of permutations requested is too small to allow for 5 possibilities per molecule. "
                          + "Switching to use Inverse instead of Brute for slsFinish for this datapont (at timeIndex of " +str(timeIndex)+ " where 0 is the first analyzed datapoint)."
                          + "Additional Info: There are " + str(len(remaining_molecules_SLS)) + " molecules that SLS was unable to solve! "
                          + "There are " + str(permutationNum) + " permutations allowed, and " + str(5**len(remaining_molecules_SLS)) + " would be needed."
                          + "If you wish to use Brute, increase the size of permutations in the user input file. ")
                    if distinguished == 'yes':#user input, choosing between distinguished inverse method or combinations method
                        concentrationsFromFinisher = InverseMethodDistinguished (remaining_reference_intensities_SLS,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS, uncertainties_dict)
                    else:
                        concentrationsFromFinisher = InverseMethod (remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_SLS,mass_fragment_numbers,remaining_molecules_SLS,'composition')               
        # the concentrations that were solved for by the Finisher are stored as a list
        # to make them easier to use and then discard in the for loop as they are added to solutions
        remainingMolecules = list(concentrationsFromFinisher.copy())
        if len(uncertainties_dict) > 0:
            remainingMolecules_concentrations_relative_uncertainties = list(uncertainties_dict['concentrations_relative_uncertainties_one_time'].transpose())
            remainingMolecules_concentrations_absolute_uncertainties = list(uncertainties_dict['concentrations_relative_uncertainties_one_time'].transpose())
            uncertainties_dict['concentrations_relative_uncertainties_one_time_after_finisher'] = uncertainties_dict['concentrations_relative_uncertainties_one_time_before_finisher'].transpose()*1.0 #initializing
            uncertainties_dict['concentrations_absolute_uncertainties_one_time_after_finisher'] = uncertainties_dict['concentrations_absolute_uncertainties_one_time_before_finisher'].transpose()*1.0  #initializing
        #adds the finisher solutions with the inital analysis solutions 
        for molecule_iiii in range(len(molecules_unedited)):
            if len(solvedmolecules) == 0:
                print("Warning: If you have chosen to use unique fragment SLS and your data has no unique fragments (not unique to any molecule), "\
                "then the program may be about to crash. If autosolver has been turned on, the program will first attempt to use SLS common and then inverse.")
            try: 
                #if the molecule wasn't solved for in the inital analysis, then it will have 0 for its solved molecules counter.
                if solvedmolecules[molecule_iiii] == 0:
                    # then add the appropriate Finisher concentration for that molecule  
                    solutions[molecule_iiii] = remainingMolecules.pop(0) 
                    try:
                        if len(uncertainties_dict) > 0: #Note: the indexing in the solutions must match the original indexing, or errors will occur.                       
                            uncertainties_dict['concentrations_relative_uncertainties_one_time_after_finisher'][molecule_iiii] = remainingMolecules_concentrations_relative_uncertainties.pop(0)
                            uncertainties_dict['concentrations_absolute_uncertainties_one_time_after_finisher'][molecule_iiii] = remainingMolecules_concentrations_absolute_uncertainties.pop(0)
                    except:
                        pass
            except IndexError:
                print("Warning: SLS could not solve this problem. If you are already using SLS Common, you can try raising the referenceValueThreshold within the feature Reference Mass Fragmentation Threshold. Alternatively, you can try using inverse.")
                solutions = numpy.array([None]) #This is just creating a numpy array with an element that has a None object, so that the main function can know that SLSMethod failed.
        if len(uncertainties_dict) > 0: #Note: the indexing in the solutions must match the original indexing, or errors will occur.  Now we take the "filled up" solutions and put everything together.                      
                uncertainties_dict['concentrations_relative_uncertainties_one_time'] = uncertainties_dict['concentrations_relative_uncertainties_one_time_after_finisher'].transpose()*1.0
                uncertainties_dict['concentrations_absolute_uncertainties_one_time'] = uncertainties_dict['concentrations_absolute_uncertainties_one_time_after_finisher'].transpose()*1.0
    return solutions
    
#this function actually calls the SLS function inside of it, because the SLS function is given a smaller array
#of correction values and a smaller array of raw signals, and the function here has a list of the numbers that 
#were eliminated from having a significant percentage. This is called the sensitivity threshold, where if the 
#value of the relative intensity of the molecule for that mass fragment that gives no raw signal is greater than 
#five, then the molecule is not present for that particular time. This whole row is then deleted instead of 
#altering the numbers in order that data is not skewed, of course, this can result in less mass fragments than 
#present molecules, which will give an error later on in the function (the inverse method section of the sls most likely)
def RawSignalThresholdFilter (distinguished,reciprocal_matching_correction_values,rawsignalsarrayline,monitored_reference_intensities,molecules,timeIndex,
                              mass_fragment_numbers,ThresholdList,answer,time,
                              conversionfactor = [],datafromcsv = [],DataRangeSpecifierlist = [],
                              SLSChoices = [],permutationNum = [],scaledConcentrationsarray = [],objectiveFunctionType = [],
                              maxPermutations = 100001, uncertainties_dict={}):

    # This is creating a local copy of the monitored_reference_intensities which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reference_intensities_filter = copy.deepcopy(monitored_reference_intensities)

    # This is creating a local copy of the reciprocal_matching_correction_values which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reciprocal_correction_factors_SLS = copy.deepcopy(reciprocal_matching_correction_values)

    # This is creating a local copy of the rawsignalsarrayline which will become
    # truncated as the molecules are solved and masses are removed
    remaining_rawsignals_SLS = copy.deepcopy(rawsignalsarrayline)

    # This is creating a local copy of 'molecules' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_molecules_SLS = copy.deepcopy(molecules)
    
    absentmolecules1 = numpy.ones([1,len(remaining_reciprocal_correction_factors_SLS[0,:])])
    absentmolecules = absentmolecules1[0]
    molecules_unedited = copy.deepcopy(molecules)
    place_holder = 0
    place_holder2 = 0
    monitored_reference_intensities_copy = monitored_reference_intensities
    rawsignalsarray_copy = remaining_rawsignals_SLS
    summation = sum(rawsignalsarray_copy)
    [rawSignalThresholdMethod,rawSignalThresholdValue,sensitivityThresholdValue,rawSignalThresholdDivider,rawSignalThresholdLimit,rawSignalThresholdLimitPercent] = ThresholdList
    #this section of the code enables the function to eliminate from the raw signal array (for the purpose of this function only)
    #the highest value in the array, if it makes up over 90 percent of the raw signals present. This is useful because if one of 
    #the signals does become this great then it will eliminate all the other signals present when it becomes very high
    if len(rawSignalThresholdLimitPercent) == 0:#if no user input
        rawSignalThresholdLimitPercent = 0.9
    if rawSignalThresholdLimit == 'yes':#from data edit file
        for rawsignalNumber in range(len(remaining_rawsignals_SLS)):    #array-indexed for loop
            if summation*rawSignalThresholdLimitPercent < max(rawsignalsarray_copy):#if the greatest value in the raw signal array is greater than 90 percent of the total sum it is eliminated
                for counternested in range(len(rawsignalsarray_copy)):#array-indexed for loop
                    if remaining_rawsignals_SLS[counternested] == max(rawsignalsarray_copy):# always checks against max of loop
                        rawsignalsarray_copy = numpy.delete(rawsignalsarray_copy,(counternested))
                        summation = sum(rawsignalsarray_copy)
    #these if statements define the thresholds if the user wants to use this feature yet does not input any values in the data
    #edit file - this sets the defaults
    if len(rawSignalThresholdValue) == 0:#user input
        if len(rawSignalThresholdDivider) == 0:#user input
            rawSignalThresholdValue = summation/float(100)
        else:#obtains raw signal threshold value
            rawSignalThresholdValue = summation/float(rawSignalThresholdDivider)
    if len(sensitivityThresholdValue) == 0:#user input
        #TODO FIXME: this needs to be a 0 or a required variable
        sensitivityThresholdValue = 5
    #TODO FIXME: we think the next two lines need to be deleted
    elif len(sensitivityThresholdValue) > 1:#user input
        sensitivityThresholdValue = sensitivityThresholdValue[0]
    #This for loop goes through all of the rows of the remaining_rawsignals_SLS and finds the values that are lower than the given 
    #threshold as defined above, and these rows are deleted from the remaining_rawsignals_SLS, the remaining_reciprocal_correction_factors_SLS and the matching
    #mass fragments array, all using a place hold so that even when the rows are deleted, the counting will still work for 
    #the loop. The a nested for loop looks in the matching mass fragments array (a copy, that is not changed in the previous
    #loop), and finds which molecules from each of those deleted rows have a relative intensity above the given threshold 
    #value in their respective mass fragmentation pattern. These molecules are saved as zeros in the absentmolecules array
    for rawcounter in range(len(remaining_rawsignals_SLS)):#array-indexed for loop
        if remaining_rawsignals_SLS[rawcounter-place_holder] < rawSignalThresholdValue:#any values below threshold are deleted
            remaining_rawsignals_SLS = numpy.delete(remaining_rawsignals_SLS,rawcounter-place_holder,axis = 0)
            remaining_reciprocal_correction_factors_SLS = numpy.delete(remaining_reciprocal_correction_factors_SLS,rawcounter-place_holder,axis = 0)
            remaining_reference_intensities_filter = numpy.delete(remaining_reference_intensities_filter,rawcounter-place_holder,axis = 0)
            if len(uncertainties_dict) > 0: #will mimic the above lines for uncertainties.          
                if type(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']) != type(None):
                    uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'] = numpy.delete(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'],rawcounter-place_holder,axis = 0)
                if type(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS']) != type(None):
                    uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'] = numpy.delete(uncertainties_dict['remaining_rawsignals_absolute_uncertainties_SLS'],rawcounter-place_holder,axis = 0)

            place_holder = place_holder + 1
            for monitored_reference_intensitiescounter in range(len(monitored_reference_intensities_copy[0,:])):#array-indexed for loop
                if monitored_reference_intensities_copy[rawcounter,monitored_reference_intensitiescounter] > sensitivityThresholdValue:#if the sensitivity is above the default/input then that mass fragment is not deleted.
                    absentmolecules[monitored_reference_intensitiescounter] = 0
    #the absent molecules array then, along with the help of a place_holder deletes all of these columns from the respective 
    #arrays and the values needed remain in this function, as the solvable leftover arrays are sent to the next function - sls
    for absentmoleculescounter in range(len(absentmolecules)):#array-indexed for loop
        if absentmolecules[absentmoleculescounter] == 0:#finds indexes where molecules can be deleted
            remaining_molecules_SLS = numpy.delete(remaining_molecules_SLS,absentmoleculescounter-place_holder2)
            remaining_reference_intensities_filter = numpy.delete(remaining_reference_intensities_filter,absentmoleculescounter-place_holder2,axis = 1)
            remaining_reciprocal_correction_factors_SLS = numpy.delete(remaining_reciprocal_correction_factors_SLS,absentmoleculescounter-place_holder2,axis = 1)
            if len(uncertainties_dict) > 0: #will mimic the above lines for uncertainties.          
                if type(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS']) != type(None):
                    uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'] = numpy.delete(uncertainties_dict['remaining_correction_values_relative_uncertainties_SLS'],absentmoleculescounter-place_holder2,axis = 1)
            place_holder2 = place_holder2 + 1    
    #this last section of code just adds the absent molecules together with the solutions to get all of the molecules needed, and
    #prints off the results using a copy of the molecules array made earlier; this part of the function is very similar to what the
    #sls function has to add its array to that of what it gets back from the inverse function that it uses.
    if any(absentmolecules) == 0:#everything is below threshold-not probable
        print('No Significant Peaks')
    else:
        if answer == 'inverse':#if the inverse method is wanted
            if distinguished == 'yes':#distinguished method
                solutions = InverseMethodDistinguished(remaining_reference_intensities_filter,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS)
            else:#combinations method
                solutions = InverseMethod(remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_filter,mass_fragment_numbers,remaining_molecules_SLS,'composition')
        if answer == 'sls' or answer == 'autosolver':#sls method
            solutions = SLSMethod(remaining_molecules_SLS,remaining_reference_intensities_filter,remaining_reciprocal_correction_factors_SLS,remaining_rawsignals_SLS,timeIndex,conversionfactor,datafromcsv,molecules_unedited,DataRangeSpecifierlist,SLSChoices,mass_fragment_numbers,permutationNum,scaledConcentrationsarray,objectiveFunctionType, time, maxPermutations)
        timeIndex = 0
        for moleculecounter in range(len(molecules_unedited)):#array-indexed for loop
            if absentmolecules[moleculecounter] == 1:#gets index for present molecule
                absentmolecules[moleculecounter] = solutions[moleculecounter-timeIndex]
            else:#keeps index in check
                timeIndex = timeIndex+1
    solutions = absentmolecules
    return solutions
    
#This little bit of code enables the user to input values for a raw signal intensity and the corresponding concentration of that molecule 
#such that the conversion between that mass fragment's signal_relative_to_CO and 
#concentration can be found, and as such the rest of the concentrations can be found as well. If no data is input here, then
#there are no concentrations printed out, only signals. (if they are both printed, it is done on separate excel sheets)
#TODO Make RatioFinder capable of using both numerous reference patterns and separate molecules (right now can only use one of these features at a time) and then remove the TODO comment with concentrationFinder in user input
#TODO continued: As seen below we have conversionFactorsAtEachTime initialized as an array of the same length as either the number of data points or number of molecules depending on which feature is being used.
#TODO continued: To make compatible with both features at the same time the array will need to be initialized as numpy.zeros((len(ExperimentData.times),len(ReferenceData[0].molecules))
#TODO continued: One solution would be to check if moleculesTSC_List is a list of lists (make this the first if statement after the line if concentrationFinder == 'yes')
#TODO continued: Initialize conversionFactorsAtEachTime (maybe change to conversionFactorArray) as an array of zeros that has the same shape as the experimental data (without the headers or abscissa headers) (maybe make a copy of experiment data or a sliced copy to remove headers) and populate it with the proper conversion factor
#TODO continued: To do so, find the first reference file's conversion factors based on the values in the first list of moleculeConcentrationTSC, moleculeSignalTSC, and moleculeTSC from the user input
#TODO continued: In conversionFactorAtEachTime, populate the rows associating the reference pattern's time ranges with the calculated conversion factors
#TODO continued: Repeat for each list in moleculeTSC_List
#TODO continued: In the case of gaps in time ranges, just interpolate the conversion factors (of like molecules) between the two time ranges
def RatioFinder (AllMoleculesReferenceDataList, AllMassFragmentsExperimentData, ReferenceData, ExperimentData, concentrationFinder,TSC_List_Type,moleculesTSC_List,moleculeConcentrationTSC_List,massNumberTSC_List,moleculeSignalTSC_List,referencePatternTimeRanges): 
    if concentrationFinder == 'yes':#user input
        if TSC_List_Type == 'MultipleReferencePatterns' or TSC_List_Type == 'UniformMolecularFactors':
            #TODO: change the variable name of conversionFactorsAtEachTime since it can now refer to numerous molecules
            conversionFactorsAtEachTime = numpy.zeros((len(ExperimentData.times),1)) #Initialize conversionFactorsAtEachTime to be an array that has 1 column and as many rows as there are data points
        elif TSC_List_Type == 'SeparateMolecularFactors':
            #TODO: change the variable name of conversionFactorsAtEachTime since it can now refer to numerous molecules
            conversionFactorsAtEachTime = numpy.zeros((1,len(ReferenceData[0].molecules))) #Initialize conversionFactorsAtEachTime to be an array that has 1 row and as many columns as there are molecules in the trimmed reference data
            
        if TSC_List_Type == 'MultipleReferencePatterns' or TSC_List_Type == 'UniformMolecularFactors':
            #initialize conversionFactorForEachReferenceFile as a numpy array with the same length as the number of reference files given
            conversionFactorForEachReferenceFile = numpy.zeros(len(ReferenceData))
            #Get the conversion factor for each reference pattern
            for referencePatternIndex in range(len(AllMoleculesReferenceDataList)): #loop through all reference patterns containing all molecules
                for moleculecounter in range(len(AllMoleculesReferenceDataList[referencePatternIndex].molecules)):#array-indexed for loop
                    for masscounter in range(len(AllMoleculesReferenceDataList[referencePatternIndex].matching_abscissa)):#array-indexed for loop
                        if moleculesTSC_List[referencePatternIndex] == AllMoleculesReferenceDataList[referencePatternIndex].molecules[moleculecounter]:#gets molecule index
                            if massNumberTSC_List[referencePatternIndex] == AllMoleculesReferenceDataList[referencePatternIndex].matching_abscissa[masscounter]:#gets index
                                #solve for the conversion factor of this reference file
                                conversionFactorForEachReferenceFile[referencePatternIndex] = (moleculeConcentrationTSC_List[referencePatternIndex]*AllMoleculesReferenceDataList[referencePatternIndex].reciprocal_matching_correction_values[masscounter,moleculecounter])/float(moleculeSignalTSC_List[referencePatternIndex]) #Use the matching correction value determined by using all molecules and mass fragments from the imported files
            #Now we need to populate conversionFactorsAtEachTime with the proper conversion factors
            if len(referencePatternTimeRanges) > 0: #If using reference pattern time chooser, loop through ExpData.times to determine which conversion factor goes where
                for timeIndex in range(len(ExperimentData.times)): #Looping through all times
                    for referencePatternTimeRangeIndex in range(len(referencePatternTimeRanges)): #Looping through all referencePatternTimeRanges to see which range the current time falls in
                        if ExperimentData.times[timeIndex] >= referencePatternTimeRanges[referencePatternTimeRangeIndex][0] and ExperimentData.times[timeIndex] <= referencePatternTimeRanges[referencePatternTimeRangeIndex][1]:
                            conversionFactorsAtEachTime[timeIndex] = conversionFactorForEachReferenceFile[referencePatternTimeRangeIndex]
                            break #Exit the for loop so the value does not get overwritten
                        #This elif needs to come before the last elif since the last elif checks the start time of the next time range and if in the last time range there is no 'next' time range to check so an index error occurs
                        elif referencePatternTimeRangeIndex == len(referencePatternTimeRanges)-1: #At the last referancePatternTimeRangeIndex if the time does not fall in any refPatternTimeRanges then the time either comes before the first time range begins or after the last time range ends
                            pass #Leave the value as a 0 #TODO Ask Ashi what to do in this scenario
                            break #Exit the for loop so the next elif statement is not read
                        elif ExperimentData.times[timeIndex] >= referencePatternTimeRanges[referencePatternTimeRangeIndex][1] and ExperimentData.times[timeIndex] <= referencePatternTimeRanges[referencePatternTimeRangeIndex+1][0]: #Check if in a gap
                            #If in a gap, linearly interpolate the conversion factor
                            conversionFactorsAtEachTime[timeIndex] = DataFunctions.analyticalLinearInterpolator(conversionFactorForEachReferenceFile[referencePatternTimeRangeIndex],conversionFactorForEachReferenceFile[referencePatternTimeRangeIndex+1],ExperimentData.times[timeIndex],referencePatternTimeRanges[referencePatternTimeRangeIndex][1],referencePatternTimeRanges[referencePatternTimeRangeIndex+1][0])
                            break #Exit the for loop so the value does not get overwritten
                    
            elif len(referencePatternTimeRanges) == 0: #User is not using RPTC so just make each value the single conversion factor that was calculated
                for index in range(len(ExperimentData.times)): 
                    conversionFactorsAtEachTime[index] = conversionFactorForEachReferenceFile[0]
            
        elif TSC_List_Type == 'SeparateMolecularFactors':
            #Default concentration factors at each molecule to match the first moleculeTSC input
            for moleculecounter in range(len(AllMoleculesReferenceDataList[0].molecules)): #loop over ALL molecules
                for masscounter in range(len(AllMoleculesReferenceDataList[0].matching_abscissa)): #loop over ALL mass fragments
                    if moleculesTSC_List[0] == AllMoleculesReferenceDataList[0].molecules[moleculecounter]: #gets index of first moleculeTSC in the reference data
                        if massNumberTSC_List[0] == AllMoleculesReferenceDataList[0].matching_abscissa[masscounter]: #Gets index of first massNumberTSC in the collected data
                            #Get the concentration factor for the first molecule listed
                            correction_values_masscounter = AllMoleculesReferenceDataList[0].matching_abscissa.index(AllMoleculesReferenceDataList[0].matching_abscissa[masscounter])
                            conversionFactorForFirstMoleculeTSC = (moleculeConcentrationTSC_List[0]*AllMoleculesReferenceDataList[0].reciprocal_matching_correction_values[correction_values_masscounter,moleculecounter])/float(moleculeSignalTSC_List[0]) #Use the matching correction value determined by using all molecules and mass fragments from the imported files
            #Overwrite all values in conversion factor at each time with the conversion factor of the first moleculeTSC
            #This for loop is looping conversionFactorsAtEachTime array which is the same length as the trimmed molecules
            for conversionIndex in range(len(conversionFactorsAtEachTime[0])): #index of 0 is needed because array is 2-D with 1 row and rows are indexed first and we want to loop over each column
                conversionFactorsAtEachTime[0][conversionIndex] = conversionFactorForFirstMoleculeTSC #index of 0 is needed because array is 2-D with 1 row and rows are indexed first
                
            #Now populate the conversion factors for the molecules that were listed
            for moleculeTSC_Index in range(len(moleculesTSC_List)): #Loop through the moleculesTSC_List
                for moleculecounter in range(len(AllMoleculesReferenceDataList[0].molecules)): #Looping through ALL molecules
                    for masscounter in range(len(AllMoleculesReferenceDataList[0].matching_abscissa)): #Looping through ALL mass fragments
                        if moleculesTSC_List[moleculeTSC_Index] == AllMoleculesReferenceDataList[0].molecules[moleculecounter]: #Gets the molecule index from all molecules
                            if massNumberTSC_List[moleculeTSC_Index] == AllMoleculesReferenceDataList[0].matching_abscissa[masscounter]: #Gets the mass fragment index from all mass fragments
                                if moleculesTSC_List[moleculeTSC_Index] in ReferenceData[0].molecules: #If the molecule is in the trimmed reference data find the index of where it appears
                                    ReferenceDataMoleculeIndex = numpy.where(numpy.array(ReferenceData[0].molecules) == numpy.array(moleculesTSC_List[moleculeTSC_Index]))[0][0] #np.where returns an array with the first element being a list of the indicies.  So using [0][0] as syntax we can pull the index out as an int assuming there are no repeats in molecule names
                                    #Solve for the new conversion factor and place it at the index of the molecule's appearance in the trimmed reference data
                                    #index of 0 is needed because array is 2-D with 1 row and rows are indexed first
                                    correction_values_masscounter = AllMoleculesReferenceDataList[0].matching_abscissa.index(AllMoleculesReferenceDataList[0].matching_abscissa[masscounter])
                                    conversionFactorsAtEachTime[0][ReferenceDataMoleculeIndex] = (moleculeConcentrationTSC_List[moleculeTSC_Index]*AllMoleculesReferenceDataList[0].reciprocal_matching_correction_values[correction_values_masscounter,moleculecounter])/float(moleculeSignalTSC_List[moleculeTSC_Index]) #Use the matching correction value determined by using all molecules and mass fragments from the imported files
                                else: #if the molecule is not in the trimmed data then just use the conversion factor of the first molecule listed which is what already populates conversionFactorAtEachTime
                                    pass
        else:
            raise ValueError('Invalid option for TSC_List_Type')
    elif concentrationFinder == 'no': #user input
        conversionFactorsAtEachTime = 0 #Originally defaulted to 0
    return conversionFactorsAtEachTime
    
    
#this function is going to be rather simple, but it will be the forward function, that simulates raw signals from the calculated
#signals that we acquired, which will later be printed out in the main() function. The inputs for this function are the signals 
#array as well as the correction values, in order to simulate the raw signals
def RawSignalsSimulation(scaledConcentrationsarray,reciprocal_matching_correction_values):
    simulateddata = numpy.zeros([len(scaledConcentrationsarray[:,0]),len(reciprocal_matching_correction_values[0][:,0])+1])#a simulated data array made the height of the signals array and the width equal to the correction arrays height
    times = scaledConcentrationsarray[:,0]#the first row of the signals array is the times
    scaledConcentrationsarray = scaledConcentrationsarray[:,1:] #because one of the lines here is the times, and it need to be removed
    for scaledtimeIndex in range(len(scaledConcentrationsarray[:,0])):#array-indexed for loop
        simulateddata[scaledtimeIndex:scaledtimeIndex+1,1:] = numpy.transpose(numpy.array(reciprocal_matching_correction_values[scaledtimeIndex]) @ numpy.array(numpy.vstack(scaledConcentrationsarray[scaledtimeIndex,:])))#the data is simulated by multiplying the matrix of correction values by the raw signals for each row. the "@" allows 'matrix multiplication' with numpy.
    simulateddata[:,0] = times #the times are added back in so they can be printed more easily
    return simulateddata

#This was written as a support function for NegativeAnalyzer in Jan 2020, but was put outside since it could be useful for other applications.
def littleSimulationFunction(concentrations_without_time, reciprocal_matching_correction_values_to_simulate):
    simulatedRawSignal = numpy.array(reciprocal_matching_correction_values_to_simulate)@numpy.array(numpy.vstack(concentrations_without_time))  #the "@" allows 'matrix multiplication' with numpy.
    return simulatedRawSignal

#a function to subtract the signals of the certain molecules, so that we can isolate the raw signals of  the non-subtracted molecules.
#This was written as a support function for NegativeAnalyzer in Jan 2020, but was put outside since it could be useful for other applications.
def subtract_simulated_signals_of_specific_molecules(moleculeIndicesToSubtract, allConcentrationsVector, all_reciprocal_matching_correction_values, rawsignalsarrayline):
    #allConcentrationsVector is a concetration vector for all molecules, and we will set to zeros for any molecule that should not be subtracted.
    #reciprocal_matching_correction_values_to_subtract is a concetration vector for all molecules and has zeros for any molecule that should not be subtracted.
    #First make a version of allConcentrationsVector and of reciprocal_matching_correction_values that has only the molecules to subtract, so we can simulate only those signals.
    concentrationsToSubtract = allConcentrationsVector*1.0
    reciprocal_matching_correction_values_to_subtract = all_reciprocal_matching_correction_values*1.0
    for moleculeIndex in range(len(allConcentrationsVector)):
        if moleculeIndex not in moleculeIndicesToSubtract:
            concentrationsToSubtract[moleculeIndex] = 0
            reciprocal_matching_correction_values_to_subtract[:, moleculeIndex] = 0 #matching correction values has indices of [massFragment, moleculeIndex]
        # else: #This else is implied, given how the vectors were created.
            # concentrationsToSubtract[moleculeIndex] = allConcentrationsVector[moleculeIndex]
            # reciprocal_matching_correction_values_to_subtract[:, moleculeIndex] = all_reciprocal_matching_correction_values[moleculeIndex]                    
    rawsignalstosubtract = littleSimulationFunction(concentrationsToSubtract, reciprocal_matching_correction_values_to_subtract)#The raw signals are simulated from the correction values and raw signals containing all molecules we don't want to consider.
    rawsignals_remaining = rawsignalsarrayline - rawsignalstosubtract
    return rawsignals_remaining



#This function is intended to add back any concentration / intensity that was removed or distorted by reference fragmentation pattern high pass filtering (which is UserChoices['minimalReferenceValue']['referenceValueThreshold']). It also adds half of the difference to the uncertainties.  As of 5/17/20, only a single recursion iteration is performed.
#The corrections are performed from largest to smallest (in percentage), with only one molecule's correction per sls mass. That is, the effect of molecule A's filtering on molecule B's concentration, for example. A molecule can be its largest self correction (A to A). Each molecule is corrected in kind, so serial (stacked) correction is among the possibilities that can occur.
#TODO: May 17 2020. Right now, we only take the LARGEST correction for each sls when deciding which sls mass to correct, then we apply that correction. But maybe we should take a single recursion of all molecules affecting? Then apply all molecules to that sls mass before moving to the next one? This would still be a single recursion, but would matter if (for example) a single molecule's sls was able to occur due to filtering out 20 other molecule's contributions at that mass.
class referenceThresholdFilterCorrectingSandbox():
    def __init__(self, solvedConcentrationsAtThisTime,referenceBeforeFiltering,referenceAfterFiltering, slsSolvedMasses):
        #The below few lines are to get the existing simulated signals for this time point.
        #We put list brackets around the arguments below because the RawSignalSimulation function actually expects an array of arrays.
        simulateddataAtThisTime = RawSignalsSimulation(numpy.array([solvedConcentrationsAtThisTime]), numpy.array([referenceAfterFiltering.reciprocal_matching_correction_values]))
        simulateddataAtThisTime1D = simulateddataAtThisTime[0] #The RawSignalsSimulation function returns a nested object, so we un-nest it.                                
        self.slsSolvedMasses = slsSolvedMasses
        self.referenceBeforeFiltering = referenceBeforeFiltering
        self.referenceAfterFiltering = referenceAfterFiltering
        self.solvedConcentrationsAtThisTimeUncorrected = solvedConcentrationsAtThisTime*1.0
        self.solvedConcentrationsCorrected1D = solvedConcentrationsAtThisTime*1.0 #just initializing. This *does* include the time value.
        self.correctedAlreadyCounter = self.solvedConcentrationsCorrected1D[1:]*0.0 #Just initializing with same length as molecules. Will put a 1 inside any index once it's corrected.  We use the 1: because we need to exclude the time column.
        
        for moleculesCounter in range(len(self.solvedConcentrationsCorrected1D)):#We loop across the number of molecules, not knowing which one we'll pick each time. Also, we may not pick them all. So this counter is not used for anything other than the number of times looping.
            self.populateVariablesFunction(self.solvedConcentrationsCorrected1D) #Get things ready to figure out which concentration needs to be decreased.
            self.calculateMaxContribution()  #This figures out which concentration needs to be decreased.
            #Now decrease the concentration that needs to be decreased.     
            if type(self.slsMassWithGreatestContributionToCorrectMoleculeIndex) != type(None):
                self.solvedConcentrationsCorrected1D[self.slsMassWithGreatestContributionToCorrectMoleculeIndex+1] = (1+self.slsMassWithGreatestContributionToCorrectContributionRatio)*self.solvedConcentrationsCorrected1D[self.slsMassWithGreatestContributionToCorrectMoleculeIndex+1]  #The +1 indexing is needed because the first column is time.  It's "1+...ratio" rather than "1-...ratio" because a negative ratio is defined as decreasing and a positive ratio as increasing.
                #Keep track of that having already been solved for.
                self.correctedAlreadyCounter[self.slsMassWithGreatestContributionToCorrectMoleculeIndex] =  1
            #Now the loop will occur again.
        
    def populateVariablesFunction(self, solvedConcentrationsCorrected):
        self.allMoleculesFilteredContributions = []
        self.allMoleculesSlsSimulatedSignals = []
        self.allMoleculesSlsMassIndices = []
        self.separatedSimulatedSignalsWithFiltering = []
        self.separatedSimulatedSignalsWithoutFiltering = []
        for moleculeIndex in range(len(self.slsSolvedMasses)):
            if self.correctedAlreadyCounter[moleculeIndex] == 0: #We only do this if the Counter's flag is still at 0 for this molecule. If it's at 1 that means any other contributions to this molecule have already occurred.
                slsSolvedMass = self.slsSolvedMasses[moleculeIndex]
                #Grab all of the masses and find the mass index for this mass.
                slsMassIndex = list(ExperimentData.mass_fragment_numbers).index(slsSolvedMass) #The reference object *does not* keep within it the mass fragment numbers of the observed signals. The mass fragments of the observed signals are an experimental quantity so are only in ExperimentData.
                #TODO: Consider that we are only going to do one recursive loop here. Maybe that should change.
                solvedConcentrationsWithoutOtherMolecules = solvedConcentrationsCorrected*0.0  #Need to keep the same shape but have zeros for other molecules.
                solvedConcentrationsWithoutOtherMolecules[moleculeIndex+1] = solvedConcentrationsCorrected[moleculeIndex+1] #the plus one is because this has "time" also.
                thisMoleculeSimulatedSignalWithReferenceWithoutFiltering = RawSignalsSimulation(numpy.array([solvedConcentrationsWithoutOtherMolecules]), numpy.array([self.referenceBeforeFiltering.reciprocal_matching_correction_values]))
                #There is a subtlety because the with filtering case must use the ORIGINAL concentration, because in the iterative way of doing things the solvedConcentrationsCorrected will no longer reflect he originally ascribed signal. We need to compare to that.
                solvedConcentrationsAtThisTimeUncorrectedWithoutOtherMolecules = self.solvedConcentrationsAtThisTimeUncorrected*0.0
                solvedConcentrationsAtThisTimeUncorrectedWithoutOtherMolecules[moleculeIndex+1] = solvedConcentrationsCorrected[moleculeIndex+1]
                thisMoleculeSimulatedSignalWithReferenceWithFilteringOriginal = RawSignalsSimulation(numpy.array([solvedConcentrationsAtThisTimeUncorrectedWithoutOtherMolecules]), numpy.array([self.referenceAfterFiltering.reciprocal_matching_correction_values]))
                thisMoleculefilteredContributions = (thisMoleculeSimulatedSignalWithReferenceWithFilteringOriginal - thisMoleculeSimulatedSignalWithReferenceWithoutFiltering)[0][1:]  #The 0 is necessary because this is actually a nested object. Then, the 1: is to remove the index that is related to time. We could have sliced before the subtraction, but we just do it after.
                thisMoleculeSlsSimulatedSignal = thisMoleculeSimulatedSignalWithReferenceWithFilteringOriginal[0][slsMassIndex+1]  #Need+1 because the simulated signal array includes the time.
                self.allMoleculesSlsSimulatedSignals.append(thisMoleculeSlsSimulatedSignal)
                self.allMoleculesFilteredContributions.append(thisMoleculefilteredContributions)
                self.allMoleculesSlsMassIndices.append(slsMassIndex)
                self.separatedSimulatedSignalsWithFiltering.append(thisMoleculeSimulatedSignalWithReferenceWithFilteringOriginal)
                self.separatedSimulatedSignalsWithoutFiltering.append(thisMoleculeSimulatedSignalWithReferenceWithoutFiltering)
            else: #The else implies correctedAlreadyCounter[moleculeIndex] == 1. if it's already been done, we'll append "none". 
                self.allMoleculesSlsSimulatedSignals.append(None)
                self.allMoleculesFilteredContributions.append(None)
                self.allMoleculesSlsMassIndices.append(None)
                self.separatedSimulatedSignalsWithFiltering.append(None)
                self.separatedSimulatedSignalsWithoutFiltering.append(None)
    def calculateMaxContribution(self):
        self.slsMassWithGreatestContributionToCorrectContributionRatio = 0
        self.slsMassWithGreatestContributionToCorrectMoleculeIndex = None
        self.slsMassWithGreatestContributionToCorrectSlsMassIndex = None
        self.slsMassWithGreatestContributionToCorrectSlsMassValue = None
        self.slsMassWithGreatestContributingMoleculeIndex = None
        #We're going to loop across each sls mass and find the maximum contribution ratio for that mass. We'll stick the ratio into maxContributionRatioThisSlsMass.
        for moleculeIndex, slsMassToCorrectFor in enumerate(self.slsSolvedMasses):#This is same length as molecules.
            if self.correctedAlreadyCounter[moleculeIndex] == 0: #We only check if the molecule hasn't been corrected yet.   
                slsSimulatedSignalInQuestion = self.allMoleculesSlsSimulatedSignals[moleculeIndex] #This pulls out a particular molecule's contribution to is 'sls signal'. We can't take the signal directly because that signal would (potentially) include other contributions.
                slsMassIndex=self.allMoleculesSlsMassIndices[moleculeIndex] #This is the mass index within the slsSimulatedSignalInQuestion
                individualMoleculeContributionsToThisMass = [] #Initializing
                maxContributionRatioThisSlsMass = 0  #Initializing
                maxContributingMoleculeIndex = None #Initializing
                #We're going to loop across all contributing molecules to find which one contributes most to this slsMass.\
                for MoleculeIndex in range(len(self.slsSolvedMasses)):#This is same length as molecules. We're going to get the contribution to this slsMass from each molecule.
                    if type(self.allMoleculesFilteredContributions[MoleculeIndex]) != type(None):
                        filteredContributionAllMasses = self.allMoleculesFilteredContributions[MoleculeIndex]
                        filteredContributionThisMass = filteredContributionAllMasses[slsMassIndex]
                        contributionRatio = filteredContributionThisMass / slsSimulatedSignalInQuestion 
                        #Here we keep track of whether this one was the greatest contributor or not.
                        if abs(contributionRatio) > abs(maxContributionRatioThisSlsMass) and abs(contributionRatio) > 0:  #need to take the absolute value because the biggest change could be a negative change, even though that's not we are normally searching for.
                            maxContributingMoleculeIndex = MoleculeIndex
                            maxContributionRatioThisSlsMass = contributionRatio
                #Now we check if this is the biggest ratio so far for *all* slsMasses or not.
                if abs(maxContributionRatioThisSlsMass) > abs(self.slsMassWithGreatestContributionToCorrectContributionRatio): #need to take the absolute value because the biggest change could be a negative change, even though that's not we are normally searching for.
                    self.slsMassWithGreatestContributionToCorrectContributionRatio= maxContributionRatioThisSlsMass
                    self.slsMassWithGreatestContributionToCorrectMoleculeIndex = moleculeIndex
                    self.slsMassWithGreatestContributionToCorrectSlsMassIndex = slsMassIndex
                    self.slsMassWithGreatestContributionToCorrectSlsMassValue = slsMassToCorrectFor
                    self.slsMassWithGreatestContributingMoleculeIndex = maxContributingMoleculeIndex
    
#this function can be useful for preventing negative concentrations. If it finds negatives it uses an optimizing finisher (grid based brute method) to search in the observed concentration range.
#The function finds other molecules that contribute to the signals that were relevant for that molecule, and then solves only those molecules.
#The baseNumberOfGridIntervals is per direction (so 5 actually means 2*5 = 10).
def NegativeAnalyzer(solutionsline,reciprocal_matching_correction_values,rawsignalsarrayline,molecules,objectiveFunctionType,maxPermutations=100001, topNcontributors=5, baseNumberOfGridIntervals = 5):
    NGstart = timeit.default_timer()
    solutionslinedata = solutionsline[1:]*1.0# gets rid of the times for our data array
    negatives = []
    negative_molecule_indexes = []
    num_MassFragments = len(reciprocal_matching_correction_values[:,0])  #The indices in reciprocal_matching_correction_values are [massfragment, moleculenumber]
    num_Molecules = len(solutionslinedata) #This is the same length as len(reciprocal_matching_correction_values[0,:])
    for solutionsIndex in range(len(solutionslinedata)): #looks through the line
        if solutionslinedata[solutionsIndex] < 0: #if there is a value below zero it keeps the value and the index
            negatives.append(solutionslinedata[solutionsIndex])
            negative_molecule_indexes.append(solutionsIndex)
    
    if len(negatives) > 0:#if there are negatives then the function runs
        for negativeXofN in range(len(negatives)):#does this for each negative value concentration. 
            #NOTE: the solutions line data *changes* within this loop. So it is actually a kind of iterative process. The first negative molecule solution will affect the later ones. The order of solution is not presently prioritized in any way.
            
            #We're going to pick a single mass fragment of the negative molecule to figure out which other molecules may have caused it to become negative.
            #First check if this molecule has been solved by a unique mass fragment, which means it was solved by SLS and tells us we should use that fragment.
            if G.massesUsedInSolvingMoleculesForThisPoint[ negative_molecule_indexes[negativeXofN] ] != 0: #the array being checked is not "0s and 1s", it has a mass fragment if one was used and the value of 0 otherwise.
                relevantMFofNegMolecule  = G.massesUsedInSolvingMoleculesForThisPoint[ negative_molecule_indexes[negativeXofN] ] 
                indexOfRelevantMFofNegMolecule = list(G.massFragmentsArrayForThisPoint).index(relevantMFofNegMolecule)
            else: #When the molecule does not have a unique mass fragment that it was solved by, then we get the biggest mass fragment. #TODO: probably should find the biggest significant mass fragment, since the biggest one may not have been the one used in solving.
                for massFragmentIndex in range(num_MassFragments):#looks through the correction values
                    if reciprocal_matching_correction_values[massFragmentIndex,negative_molecule_indexes[negativeXofN]] == max(reciprocal_matching_correction_values[:,negative_molecule_indexes[negativeXofN]]):#finds the index of the negative molecule's largest correction value 
                        indexOfBiggestMFofNegMolecule = massFragmentIndex
                        indexOfRelevantMFofNegMolecule = indexOfRelevantMFofNegMolecule
            
            #Now, we're going to find out which molecule(s) could have caused that "false" negative concentration. We'll focus on the biggest possible contributer(s)
            #FIXME: Below is using biggest concentration value. That's not the actual biggest contributor. Should be by simulated signal. So should simulate each molecule's contribution separately to this mass, and then find the one which has the maximum contribution. That will give the right choice for correction2index
            moleculesWithMassFragmentList = []
            for matchCorrIndexRow in range(num_Molecules):#goes through the correction values.  #num_Molecules is the same length as len(reciprocal_matching_correction_values[0,:])
                if reciprocal_matching_correction_values[indexOfRelevantMFofNegMolecule,matchCorrIndexRow] != 0:#if the molecule has a relative intensity (other than zero) at the mass fragment chosen (by the last loop)
                    moleculesWithMassFragmentList.append(1)
                else: #if there is no molecule a zero is appended to the list
                    moleculesWithMassFragmentList.append(0)
            moleculesWithMassFragmentArrayKey = numpy.array(moleculesWithMassFragmentList)
            solutionsline_MassFragmentPresentKey = solutionslinedata*moleculesWithMassFragmentArrayKey #the ones and zeros list is multiplied by the solutions, allowing only molecules with the mass fragment desired to be selected later
             
            #THIS CODE BLOCK IS NOW DEPRECATED AS OF JAN 2020, SKIP TO BELOW.
            # dominant_concentration_way = True
            # if dominant_concentration_way == True:
                # #HERE IS THE EXISTING CODE FOR DOMINANT CONCENTRATION
                # for moleculeConcentrationIndex in range(len(solutionslinedata)):#goes through the solution line, which means it's looking at each molecule/concentration.
                    # if max(solutionsline_MassFragmentPresentKey) != 0:#if there are any concentrations that are above zero
                        # if solutionsline_MassFragmentPresentKey[moleculeConcentrationIndex] == max(solutionsline_MassFragmentPresentKey):#the molecule with the highest concentration at the same mass fragment being investigated is kept for consideration.
                            # dominant_concentration_index = moleculeConcentrationIndex
                    # else:
                        # if solutionslinedata[moleculeConcentrationIndex] == max(solutionslinedata):#if there are none with that mass fragment, the highest solution (concentration) is chosen
                            # dominant_concentration_index = moleculeConcentrationIndex
                # dominant_concentration = solutionslinedata[dominant_concentration_index] #the higher (positive) concentration molecule is chosen to set a dominant concentration that could overwhelm the negative one, and will be used to set bounds as well.
                # #now we will make an array of these correction values of the two molecules chosen.
                # reciprocal_matching_correction_values_dominant_concentration_only = reciprocal_matching_correction_values[:,negative_molecule_indexes[negativeXofN]],reciprocal_matching_correction_values[:,dominant_concentration_index]#an array amalgam is made with  two columns for the two chosen molecules
                # reciprocal_matching_correction_values_dominant_concentration_only = numpy.array(reciprocal_matching_correction_values_dominant_concentration_only)
                # reciprocal_matching_correction_values_dominant_concentration_only = numpy.transpose(reciprocal_matching_correction_values_dominant_concentration_only) #the array is transposed so it can be used in matrix multiplication

                # #Now, a 'clever' algorithm is used in which we assume that the positive concentration molecule has caused the negative concentration molecule's negative concentration.
                # #In order to find a 'nice' solution, we will subtract the contributions of all other molecules, so we can analyze the signals that would be associated with these two molecules.
                # #To simulate the contribution of the other molecules (for subtraction) we'll set the concentrations of 'these' molecules to zero, and their correction values to zero also.
                # allMoleculeIndices = numpy.arange(num_Molecules)
                # indicesOfMoleculesToKeep = [negative_molecule_indexes[negativeXofN], dominant_concentration_index]
                # indicesOfMoleculesToKeep = list( dict.fromkeys(indicesOfMoleculesToKeep) ) #This removes duplicates, in case the negative_molecule is among the top matching ones.                
                
                # moleculeIndicesToSubtractList = list(allMoleculeIndices) #we start with a list of all molecules and then start removing some from the list.
                # for moleculeIndex in indicesOfMoleculesToKeep:
                    # moleculeIndicesToSubtractList.remove(moleculeIndex)
                # #Now we call a function to subtract the signals of the other molecules, so that we can isolate the raw signals of these molecules.
                # rawsignals_of_considered_molecules = subtract_simulated_signals_of_specific_molecules(moleculeIndicesToSubtractList, solutionslinedata, reciprocal_matching_correction_values,rawsignalsarrayline)
                
                # #For brute optimization, the 'specifications; are the min and max along an axis followed by the size of increments/spacing along the axis. The negative molecule is checked for between the higher molecule's signal/10, in 10ths of the range.
                # max_concentration_for_negative_molecule = dominant_concentration/float(10)
                # specifications = [(0,max_concentration_for_negative_molecule,max_concentration_for_negative_molecule/float(10)),(dominant_concentration/float(2),dominant_concentration*2,dominant_concentration*0.10)]#the specifications array is made here, with the higher molecule being checked in ten places within a factor of 2 of itself
                # if sum(specifications[0]) == 0 and sum(specifications[1]) == 0: #This can only happen if the higher of the two concentrations being considered is already 0.
                    # solutionslinedata[negativeXofN] = 0 #So we just set the concentration of our molecule in question to 0 and we skip the brute.
                # else: #The normal case is a brute optimization. Below, we pass the argument of molecules but that argument is not actually used within the optimizer.
                    # answers = OptimizingFinisher(molecules,specifications,reciprocal_matching_correction_values_dominant_concentration_only,rawsignals_of_considered_molecules,objectiveFunctionType,maxPermutations)#brute method used- 200 permutations- 20*10 from the increments above
                    # solutionslinedata[negative_molecule_indexes[negativeXofN]] = answers[0]#sets the first solution
                    # solutionslinedata[dominant_concentration_index] = answers[1]#sets the second

            dominant_contributor_way = True
            if dominant_contributor_way == True:
                #HERE IS THE NEW CODE FOR DOMINANT SIGNAL CONTRIBUTOR
                contributors_with_mass_fragment_rawsignals_list = []
                dominant_contributors_concentrations_list = []
                reciprocal_matching_correction_values_dominant_contributors_only = []
                solutionslinedata_dominant_contributors_only = []
                molecules_having_mass_fragment_indices = []
                key_for_negativeXofN_in_dominant_contributors_only = []
                for moleculeIndex in range(num_Molecules): #We will check each molecule.
                    if solutionsline_MassFragmentPresentKey[moleculeIndex] == 0:  #A zero in this array means that the molecule does not have that mass fragment, so cannot be a contributor.
                        contributors_with_mass_fragment_rawsignals_list.append(0)
                        dominant_contributors_concentrations_list.append(0)
                    if solutionsline_MassFragmentPresentKey[moleculeIndex] != 0: # #A non zero in this array means that the molecule does  have that mass fragment, so can be a contributor.
                        molecules_having_mass_fragment_indices.append(moleculeIndex) #There is probably a better way to do this than putting it in the loop, but for now it's convenient.
                        #we're going to make a solutionslinedata version that has *only* this molecule present. First make a version with all zeros, then put this molecule's concentration in.
                        solutionslinedata_this_molecule_only = solutionslinedata*0.0
                        solutionslinedata_this_molecule_only[moleculeIndex] = solutionslinedata[moleculeIndex]
                        all_simulated_raw_signals_this_molecule_only = littleSimulationFunction(solutionslinedata_this_molecule_only, reciprocal_matching_correction_values)
                        relevant_simulated_raw_signal_this_molecule_only = float(all_simulated_raw_signals_this_molecule_only[indexOfRelevantMFofNegMolecule]) #The "float" casting is to go from 1 value numpy matrices to floats.
                        contributors_with_mass_fragment_rawsignals_list.append(relevant_simulated_raw_signal_this_molecule_only) #note that this does already include the molecule chosen. 
                        dominant_contributors_concentrations_list.append(solutionsline_MassFragmentPresentKey[moleculeIndex])
                        reciprocal_matching_correction_values_dominant_contributors_only.append(reciprocal_matching_correction_values[:,moleculeIndex])
                        solutionslinedata_dominant_contributors_only.append(solutionslinedata[moleculeIndex])
                        if moleculeIndex != negativeXofN: #We append to find the useful information of where the "negative" concentration we're looking is within the list of dominant contributors (really all contributors).
                            key_for_negativeXofN_in_dominant_contributors_only.append(0)
                        if moleculeIndex == negative_molecule_indexes[negativeXofN]:
                            key_for_negativeXofN_in_dominant_contributors_only.append(1)
                #now we are going to find the indicdes that we want, in a descending manner (that is, by descending level of contribution)
                sorted_contributors_with_mass_fragment_molecular_indices_ascending = numpy.argsort(contributors_with_mass_fragment_rawsignals_list) #unfortunately, this has no reverse option, so we will flip it.
                sorted_contributors_with_mass_fragment_molecular_indices_descending = numpy.flipud(sorted_contributors_with_mass_fragment_molecular_indices_ascending) #Ideally we could use fliplr, but numpy treats 1D arrays weirdly, this time we don't want to use the trick of making it 2D.
                #Now, we will only keep the concentrations of the topNcontributors or less.
                molecularIndicesOfTopNorLessContributors = []
                concentrations_topNorLess_contributors_only = []
                for majorContributorNumber in range(0,topNcontributors): #this loop works correctly for array indexing.
                    currentContributorMolecularIndex = sorted_contributors_with_mass_fragment_molecular_indices_descending[majorContributorNumber] #Get the molecule index
                    if currentContributorMolecularIndex in molecules_having_mass_fragment_indices: #Molecules can get into the top N even with a value of 0 raw signal contribution, so we will only even consider it for optimizing if it has the mass fragment index. 
                        molecularIndicesOfTopNorLessContributors.append(currentContributorMolecularIndex)
                        concentrations_topNorLess_contributors_only.append(solutionslinedata[currentContributorMolecularIndex]) #Keep the data from that particular index.

                #Now, a 'clever' algorithm is used in which we assume that dominant contributing molecules have caused the negative concentration molecule's negative concentration.
                #In order to find a 'nice' solution, we will subtract the contributions of all other molecules from the raw signals, so we can analyze the signals that would be associated with these relevant molecules.
                #we are going to subtract the contributions from everything except the molecularIndicesOfTopNorLessContributors and the negative molecule.
                allMoleculeIndices = numpy.arange(num_Molecules)
                indicesOfMoleculesToKeep = copy.deepcopy(molecularIndicesOfTopNorLessContributors)
                indicesOfMoleculesToKeep.append(negative_molecule_indexes[negativeXofN]) #we take the topNorLess and add the negative molecules.
                indicesOfMoleculesToKeep = list( dict.fromkeys(indicesOfMoleculesToKeep) ) #This removes duplicates, in case the negative_molecule is among the top matching ones.                

                moleculeIndicesToSubtractList = list(allMoleculeIndices) #we start with a list of all molecules and then start removing some from the list.
                for moleculeIndex in indicesOfMoleculesToKeep:
                    moleculeIndicesToSubtractList.remove(moleculeIndex)
                #Now we call a function to subtract the signals of the other molecules, so that we can isolate the raw signals of these molecules.
                rawsignals_of_considered_molecules_topNorLess_contributors = subtract_simulated_signals_of_specific_molecules(moleculeIndicesToSubtractList, solutionslinedata, reciprocal_matching_correction_values,rawsignalsarrayline)
                #Now we have the raw signals that we need to go into brute optimization.
                
                #Now we need the concentrations of *only* the dominant contributors and the negative molecule. To make sure we don't double count, we'll check if the negative molecule ended up in the topNorLess and work with it separately if it ended up there.
                if negative_molecule_indexes[negativeXofN] not in molecularIndicesOfTopNorLessContributors:
                    concentrations_topNorLess_contributors_only_without_negativeXofN = concentrations_topNorLess_contributors_only #No removal needed.
                    molecularIndicesOfTopNorLessContributors_without_negativeXofN = molecularIndicesOfTopNorLessContributors
                if negative_molecule_indexes[negativeXofN] in molecularIndicesOfTopNorLessContributors:
                    index_of_negativeXofN_in_molecularIndicesOfTopNorLessContributors = molecularIndicesOfTopNorLessContributors.index(negative_molecule_indexes[negativeXofN])
                    concentrations_topNorLess_contributors_only_without_negativeXofN = list(copy.deepcopy(concentrations_topNorLess_contributors_only)) # First make a copy, then remove on next line
                    concentrations_topNorLess_contributors_only_without_negativeXofN.pop(index_of_negativeXofN_in_molecularIndicesOfTopNorLessContributors)
                    concentrations_topNorLess_contributors_only_without_negativeXofN = numpy.array(concentrations_topNorLess_contributors_only_without_negativeXofN)
                    molecularIndicesOfTopNorLessContributors_without_negativeXofN = copy.deepcopy(molecularIndicesOfTopNorLessContributors)
                    molecularIndicesOfTopNorLessContributors_without_negativeXofN.pop(index_of_negativeXofN_in_molecularIndicesOfTopNorLessContributors)

                if sum(concentrations_topNorLess_contributors_only_without_negativeXofN) == 0: #This is for exceptions and can only happen if all of the contributing concentrations being considered are already 0.
                     solutionslinedata[negativeXofN] = 0 #In this case we just set the concentration of our molecule in question to 0 and we skip the brute.
                else: #This is for the normal cases, including when some fo the contributors are negative already.                    
                    #Now to set the grid limits for brute optimization. We'll take the maximum of the contributing concentrations divided by 10 for the possible positive size increments of the negative molecule.
                    #Will make the  specifications into a list, because I need to separate them. Each is supposed to be a tuple of min, max, spacing/increment.
                    specifications_topNorLess_contributors = [] #This will include the negative molecule.
                    max_concentration_for_negative_molecule = max(solutionslinedata_dominant_contributors_only)/float(baseNumberOfGridIntervals)
                    specifications_topNorLess_contributors.append(  (0, max_concentration_for_negative_molecule,max_concentration_for_negative_molecule/float(10)) )
                    #Now a loop with decreasing numbers of intervals as we go down in "rank".
                    currentNumberOfGridIntervals = baseNumberOfGridIntervals
                    for topNorLessConcentration in concentrations_topNorLess_contributors_only_without_negativeXofN: #now loop across the topNorLessContributors indices and make 5 points per grid side.
                        currentNumberOfGridIntervals = currentNumberOfGridIntervals - 1 #decreases by one for each loop.
                        if currentNumberOfGridIntervals > 10: #can't go above 10, because that can cause negative values with current algorithm of checking by 10% intervals.
                            currentNumberOfGridIntervals = 10
                        if currentNumberOfGridIntervals < 1: #can't go below 1.
                            currentNumberOfGridIntervals = 1
                        currentIntervalSize = topNorLessConcentration*0.10
                        specifications_topNorLess_contributors.append(  (topNorLessConcentration-(currentNumberOfGridIntervals*currentIntervalSize), topNorLessConcentration+(currentNumberOfGridIntervals*currentIntervalSize), currentIntervalSize)  ) #(min, max, spacing/increment) in tuple.

                    #now need to make a new correction values array. 
                    reciprocal_matching_correction_values_topNorLess_contributors_and_negative_molecule = []
                    reciprocal_matching_correction_values_topNorLess_contributors_and_negative_molecule.append(reciprocal_matching_correction_values[:,negative_molecule_indexes[negativeXofN]])#first add the negative molecule because it's simple.
                    for topNorLessIndex in molecularIndicesOfTopNorLessContributors_without_negativeXofN: #now loop across the topNorLessContributors indices and add their correction values too.
                        reciprocal_matching_correction_values_topNorLess_contributors_and_negative_molecule.append(reciprocal_matching_correction_values[:,topNorLessIndex])
                    reciprocal_matching_correction_values_topNorLess_contributors_and_negative_molecule = numpy.array(reciprocal_matching_correction_values_topNorLess_contributors_and_negative_molecule)
                    reciprocal_matching_correction_values_topNorLess_contributors_and_negative_molecule = numpy.transpose(reciprocal_matching_correction_values_topNorLess_contributors_and_negative_molecule) #the array is transposed to from list of 1D rows to 1D columns (or vice versa)

                    answers_topNorLessContributors = OptimizingFinisher(molecules,specifications_topNorLess_contributors,reciprocal_matching_correction_values_topNorLess_contributors_and_negative_molecule,rawsignals_of_considered_molecules_topNorLess_contributors,objectiveFunctionType,maxPermutations)#brute method used- 200 permutations- 20*10 from the increments above
                    #now put things into the solutions line data.
                    solutionslinedata[negative_molecule_indexes[negativeXofN]] = answers_topNorLessContributors[0] #The negative molecule's information is here.
                    for topNorLessIndex, topNorLessConcentration in enumerate(answers_topNorLessContributors[1:]):  #we start from 1 because we've already got the negative molecule's information. Now topNorLessIndex of 0 is the first of these concentrations.
                        moleculeIndex = molecularIndicesOfTopNorLessContributors_without_negativeXofN[topNorLessIndex]
                        solutionslinedata[moleculeIndex] = topNorLessConcentration

    solutionsline[1:] = solutionslinedata #sets the new data, with the times
    NGtimeSpent = (timeit.default_timer() - NGstart)
    if NGtimeSpent > 10:
        print("Time Spent for Negative Analyzer on this data point", NGtimeSpent, "original data was:", rawsignalsarrayline)
    return solutionsline

'''This function exports all XYYY Data to the User specified document (usually a CSV)'''
#TODO Future Development: This function could benefit from creating folders to store the different
#runs of the program if it is run in the same directory as before
## NOTE: This function leaves a tailing comma on each line of the
## created csv file. When pandas is used to read the csv this
## will result in the creation of a column of 'nan'
## ImportWorkingData() has been modified to remove this column
## If the data header and the data are the same length, the abscissa header is disregarded.
def ExportXYYYData(outputFileName, data, dataHeader, abscissaHeader = 'Mass', fileSuffix = '', dataType = None, rowIndex = [], units = None): 
    formatedDataHeader = dataHeader
    if dataType == 'preProcessed' or dataType == 'simulated' or dataType == 'Experiment':
        formatedDataHeader = ['m%s' % MFNumber for MFNumber in dataHeader]
    if dataType == 'scaled':
        formatedDataHeader = ['%s Concentration Relative to CO' % molecule for molecule in dataHeader]
    if dataType == 'concentration':
        label = 'Concentration(%s)' % units
        formatedDataHeader = [molecule + ' ' + label for molecule in dataHeader]
    if dataType == 'percent_concentration': #deprecated, used to be used for relative_uncertainty
        label = 'Percent' 
        formatedDataHeader = [molecule + label for molecule in dataHeader]
    if dataType == 'relative_uncertainty':
        label = 'relative uncertainty' 
        formatedDataHeader = [molecule + ' ' + label for molecule in dataHeader]
    #extraLine is used to create CSV files that conform to MSRESOLVE's import requirements i.e. having a row for comments at the top
    extraLine = False
    if dataType == 'Experiment':
        extraLine = len(data[0,1:])
        
    #If future applications of Export XYYY are desired, the new formats can be 
    #specified by additional keywords and if statements.

    #if iterative analysis is being used and the suffix is wanted
    if not fileSuffix =='':
        #then the filename will have a suffix attached
        outputFileName = outputFileName[:-4] + fileSuffix + outputFileName[-4:]

    #testing if file is open, and rename if it is
    #create list of name options
    nameOptions = [''] + list(range(100))
    for x in nameOptions:
        # create new name
        filename = (outputFileName[:-4] + '%s' + outputFileName[-4:]) %x
        #Test if it can be opened
        try:
            open(filename, 'w')  
            break
        except(IOError):
            pass
    #combine the column headers and data into one array
    try:
        fullArrayToExport = numpy.vstack((formatedDataHeader,data))
    #occasionally, abscissaHeader needs to be inserted without rowIndex being used
    except ValueError: 
        formatedDataHeader = numpy.hstack((abscissaHeader,formatedDataHeader))
        fullArrayToExport = numpy.vstack((formatedDataHeader,data))
        
    #if the row index is available, then add it 
    if len(rowIndex) > 0:    
        abscissaHeader = numpy.transpose((numpy.array([[abscissaHeader]])))
        rowIndex = numpy.transpose([rowIndex])
        abscissaArrayToExport =  numpy.vstack((abscissaHeader,rowIndex))
        fullArrayToExport = numpy.hstack((abscissaArrayToExport,fullArrayToExport))
    #insert an extra line with a header of the data type. Included to allow exported files to be uploaded during iterative analysis.
    if not extraLine == False:
        lineToInsert = "%s,%s" %(dataType, ',' * (extraLine))
        lineToInsert = numpy.array(lineToInsert.split(','))
        fullArrayToExport = numpy.vstack((lineToInsert, fullArrayToExport))
    #save the file to the correct name
    
    delimiter = getDelimiterFromExtension(filename)
    numpy.savetxt(filename, fullArrayToExport, delimiter = delimiter, fmt ="%s")
  

'''This function inserts rows of percentages into arrays of data'''
def GeneratePercentages(scaledConcentrationsarray):
    #FYI GeneratePercentages function can currently only deal with data sets that 
    #contain times 
    
    #find shape of data array
    cols = len(scaledConcentrationsarray[0,:]) 
    rows = len(scaledConcentrationsarray[:,0])
    #initilize size of output array
    newsignalsarray = numpy.zeros(cols)
    #for all data points
    for row in range(rows):
        #reset sum
        sum = 0
        #insert time in [0] position
        percentagesArray = [scaledConcentrationsarray[row,0]]
        #for loop calculates sum of a row
        for col in range(1, cols): #starts at 1 to not include time
            sum = sum + scaledConcentrationsarray[row,col]
        #for loop constructs an array of the percentages
        for col in range(1, cols): #starts at 1 to not include time
            percentagesArray.append((scaledConcentrationsarray[row,col]/sum)*100)
        #add row of percentages
        newsignalsarray = numpy.vstack((newsignalsarray, percentagesArray))    
    #return complied array, excluding size-setting row
    return newsignalsarray[1:,:]  

'''
This function is the standard function used to graph 
molecule concentrations or mass fragments signals. 
#TODO: Update Draw to use the figure objects and set functions like in DrawUncertainties (that is the newer way of doing things).
'''
def Draw(times, data, molecules, concentrationFinder, units, graphFileName = '', fileSuffix = '', label="000", stopAtGraphs = True, figureNumber=1):
    import matplotlib.pyplot as plt
    plt.figure(figureNumber) #should number figure before doing more (answer with wording "when you call") https://stackoverflow.com/questions/6916978/how-do-i-tell-matplotlib-to-create-a-second-new-plot-then-later-plot-on-the-o
    if concentrationFinder == 'yes':
        plt.ylabel('Concentration (%s)'%(units))
    else:
        plt.ylabel('Concentration Relative to CO')
    colormap = plt.cm.gist_ncar
    colorListNumbers = numpy.linspace(0.0,0.9,len(data[0,:])) #these choices were made to be aesthetically pleasing. The len part is the number of data series. The 0.9 cuts off part of the color spectrum  (can go from 0 to 1)
    colorList = []
    for color in colorListNumbers:
        colorList.append(colormap(color))
    ax = plt.subplot(111, label=label)
    for x in range(len(data[0,:])):
        ax.plot(times,data[:,x], color=colorList[x]) #In 2019, needed to change how colors were implemented due to plt change.
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(molecules,loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel(ExperimentData.abscissaHeader) #Get the x-axis title from the collected data csv file

    # Now save the plot as a png
    # assuming a file name was passed for it
    if graphFileName != '':
        
        #if a file suffix has been provided, append it to the file name
        if not fileSuffix == '':
            graphFileName =  graphFileName + fileSuffix 
            
        # need to save to plot directory
        # directory containing MSRESOLVE.py
        currentDirectory = os.path.dirname(os.path.realpath(__file__))
        
        # subdirectory for plots
        graphDirectory = os.path.join(currentDirectory, 'Graphs')
        
        # make sure it exists and create it if not
        if not os.path.exists(graphDirectory):
            os.makedirs(graphDirectory)
            
        plt.savefig(os.path.join(graphDirectory, graphFileName))
    #Need to check for ipython because it affects how graphs will be displayed.
    #based below function on https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    def checkForIpython():
        try:
            cfg = get_ipython().config 
            return True
        except NameError:
            return False        
    #Figured out "block" argument works. Got it from: https://community.esri.com/thread/185110-matplotlib-show-prevents-script-from-completing
    if checkForIpython()==True:
        plt.show(block=False)
    if checkForIpython()==False:
        if stopAtGraphs == False:
            plt.show(block=False)
        if stopAtGraphs == True:
            plt.show(block=True)                                       
    return figureNumber


'''
This function is the standard function used to graph 
molecule concentrations or mass fragments signals with uncertanties.
'''
def DrawUncertaintiesConcentrations(times, data, data_uncertainties, molecules, concentrationFinder, units, graphFileName = '', fileSuffix = '', label="000", stopAtGraphs = True, figureNumber=1):
    import matplotlib.pyplot as plt
    numConcentrations = len(molecules)
    colormap = plt.cm.gist_ncar
    colorListNumbers = numpy.linspace(0.0,0.9,len(data[0,:])) #these choices were made to be aesthetically pleasing. The len part is the number of data series. The 0.9 cuts off part of the color spectrum  (can go from 0 to 1)
    colorList = []
    for color in colorListNumbers:
        colorList.append(colormap(color))
    
    figureNumber = figureNumber
    for concentrationIndex in range(numConcentrations):
        thisFigure = plt.figure(figureNumber) #should number figure before doing more (answer with wording "when you call") https://stackoverflow.com/questions/6916978/how-do-i-tell-matplotlib-to-create-a-second-new-plot-then-later-plot-on-the-o
        thisFigure.ax = thisFigure.add_subplot(111, label=label)
        moleculeName = str(molecules[concentrationIndex])
        
        if concentrationFinder == 'yes':
            thisFigure.ax.set(  ylabel = str(molecules[concentrationIndex])+'\n Concentration (%s)'%(units))
        else:
            thisFigure.ax.set(  ylabel = str(molecules[concentrationIndex])+'\n Concentration Relative to CO')

        thisFigure.ax.set(xlabel = ExperimentData.abscissaHeader) #Get the x-axis title from the collected data csv file
        thisFigure.ax.errorbar(times, data[:,concentrationIndex], yerr=data_uncertainties[:,concentrationIndex], color=colorList[concentrationIndex], ecolor="lightgray")
        box = thisFigure.ax.get_position()
        thisFigure.ax.set_position([box.x0 +0.1, box.y0, box.width * 0.8, box.height])
        thisFigure.ax.legend(moleculeName)# ,loc='center left', bbox_to_anchor=(1, 0.5)) # Could not get text in legend to work correctly.

        # Now save the plot as a png
        # assuming a file name was passed for it
        if graphFileName != '':
            #if a file suffix has been provided, append it to the file name
        
            graphFileNameFinal =  graphFileName + "_" + str(moleculeName) + fileSuffix 
                
            # need to save to plot directory
            # directory containing MSRESOLVE.py
            currentDirectory = os.path.dirname(os.path.realpath(__file__))
            
            # subdirectory for plots
            graphDirectory = os.path.join(currentDirectory, 'Graphs')
            
            # make sure it exists and create it if not
            if not os.path.exists(graphDirectory):
                os.makedirs(graphDirectory)
                
            plt.savefig(os.path.join(graphDirectory, graphFileNameFinal))
        #Need to check for ipython because it affects how graphs will be displayed.
        #based below function on https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
        def checkForIpython():
            try:
                cfg = get_ipython().config 
                return True
            except NameError:
                return False        
        #Figured out "block" argument works. Got it from: https://community.esri.com/thread/185110-matplotlib-show-prevents-script-from-completing
        if checkForIpython()==True:
            plt.show(block=False)
        if checkForIpython()==False:
            if stopAtGraphs == False:
                plt.show(block=False)
            if stopAtGraphs == True:
                plt.show(block=True)                                       
        #In the end, need to increment the figure number after each cycle.
        figureNumber = figureNumber+1
    
    return figureNumber



'''This function is called to create the log file and to record the time'''
def CreateLogFile():
    filename6 = 'LogFile.txt'
    f6 = open(filename6,'w')
    f6.write(time.ctime(time.time()))
    f6.write('\n')
    f6.close()


# Since SlSUniqueFragments is potentially used in a number of the analysis options
# set up the 'SLSUniqueOrder.csv' file headers here
# note that this open overwrites previous file contents
# while subsequent opens to this file will append
## Parameters:
# abscissaHeader - string name of the experiment data type (e.g. 'Temp' or 'time')
# molecules - list of strings of molecule names from referenceData.molecules
def createSLSUniqueExportFiles(abscissaHeader, molecules):
    outputMoleculeOrderFileName = 'ExportedSLSUniqueMoleculesOrder.csv'
    if G.iterativeAnalysis:
        #then the filename will have a suffix attached
        outputMoleculeOrderFileName = outputMoleculeOrderFileName[:-4] + '_iter_%s' %G.iterationNumber + outputMoleculeOrderFileName[-4:]
    
    with open(outputMoleculeOrderFileName,'w') as fp:
        # Headers
        fp.write('{},{}'.format(
            'Data Point',abscissaHeader))
        # in writing molecule names in the header here
        # I assume the actual numbers are written in the same
        # order later 
        for molecule in molecules:
            fp.write(",{}".format(molecule))
        fp.write('\n')
    
    outputMassFragmentOrderFileName = 'ExportedSLSUniqueMassFragmentsUsed.csv'
    if G.iterativeAnalysis:
        #then the filename will have a suffix attached
        outputMassFragmentOrderFileName = outputMassFragmentOrderFileName[:-4] + '_iter_%s' %G.iterationNumber + outputMassFragmentOrderFileName[-4:]    
    
    with open(outputMassFragmentOrderFileName, 'w') as f:
        # Headers
        f.write('{},{},'.format(
            'Data Point',abscissaHeader))        
        f.write(str(list(ExperimentData.mass_fragment_numbers))[1:-1] + "\n") #The [1:-1] is to get rid of the brackets from the list during export.

'''
This function takes in the end result of everything and exports it to the 
log file. 
'''
def PopulateLogFile():
    filename6 = 'LogFile.txt' #the log file is printed here
    f6 = open(filename6,'a')
    f6.write('\n')
    f6.write('referenceFileName = %s \n'%(G.referenceFileNamesList))
    f6.write('form  = %s \n'%(G.referenceFormsList))
    f6.write('collectedFileName = %s \n'%(G.collectedFileName ))
    if G.timeRangeLimit == 'yes':#some of the lines in the backgroundinput file don't need to be printed unless a selection is made, so the if statements here make that happen
        f6.write('timeRangeLimit = %s \n'%(G.timeRangeLimit))
        f6.write('timeRangeStart = %s \n'%(G.timeRangeStart))
        f6.write('timeRangeFinish = %s \n'%(G.timeRangeFinish))
    if G.backgroundMassFragment != []:#some of the lines in the backgroundinput file don't need to be printed unless a selection is made, so the if statements here make that happen
        f6.write('backgroundMassFragment = %s \n'%(G.backgroundMassFragment))
        f6.write('backgroundSlopes = %s \n'%(G.backgroundSlopes))
        f6.write('backgroundIntercepts = %s \n'%(G.backgroundIntercepts))
    if G.measuredReferenceYorN == 'yes':
        f6.write('measuredReferenceYorN = %s \n'%G.measuredReferenceYorN)
        f6.write('referenceCorrectionCoefficientA = %s \n'%(G.referenceCorrectionCoefficients[0]))
        f6.write('referenceCorrectionCoefficientB = %s \n'%(G.referenceCorrectionCoefficients[1]))
        f6.write('referenceCorrectionCoefficientC = %s \n'%(G.referenceCorrectionCoefficients[2]))
        f6.write('referenceCorrection_abcCoefficients = %s \n'%(list(G.referenceCorrectionCoefficients)))
        if type(G.referenceCorrectionCoefficients_cov) == type(None):
            strVersionOf_abcCoefficients_cov = str(G.referenceCorrectionCoefficients_cov)
        if type(G.referenceCorrectionCoefficients_cov) != type(None):
            if numpy.sum(numpy.array(G.referenceCorrectionCoefficients_cov)) == 0:
                listVersionOf_abcCoefficients_cov = [[0,0,0],[0,0,0],[0,0,0]]
            else: listVersionOf_abcCoefficients_cov = list(G.referenceCorrectionCoefficients_cov) #make the outer structure a list.
            for elementIndex in range(len(listVersionOf_abcCoefficients_cov)): #make the inner structures lists, also.
                listVersionOf_abcCoefficients_cov[elementIndex] = list(listVersionOf_abcCoefficients_cov[elementIndex])
            strVersionOf_abcCoefficients_cov = str(listVersionOf_abcCoefficients_cov)
        f6.write('referenceCorrection_abcCoefficients_cov = %s \n'%(strVersionOf_abcCoefficients_cov))
    if G.specificMolecules == 'yes':
        f6.write('specificMolecules = %s \n'%(G.specificMolecules))
        f6.write('chosenMoleculesNames = %s \n'%(G.chosenMoleculesNames))
    if G.specificMassFragments == 'yes':
        f6.write('specificMassFragments = %s \n'%(G.specificMassFragments))
        f6.write('chosenMassFragments = %s \n'%(G.chosenMassFragments))
    if G.linearBaselineCorrectionSemiAutomatic == 'yes':
        f6.write('linearBaselineCorrectionSemiAutomatic = %s \n'%(G.linearBaselineCorrectionSemiAutomatic))
        f6.write('selection = %s \n'%(G.baselineType))
        f6.write('massesToBackgroundCorrect = %s \n'%(G.massesToBackgroundCorrect))
        f6.write('earlyBaselineTimes = %s \n'%(G.earlyBaselineTimes))
        f6.write('lateBaselineTimes = %s \n'%(G.lateBaselineTimes))
    if G.dataRangeSpecifierYorN == 'yes':
        f6.write('dataRangeSpecifierYorN = %s \n'%(G.dataRangeSpecifierYorN ))
        f6.write('signalOrConcentrationRange = %s \n'%(G.signalOrConcentrationRange))
        f6.write('csvFile = %s \n'%(G.csvFile))
        f6.write('moleculesRange = %s \n'%(G.moleculesToRestrict))
        if G.csvFile == 'yes':
            f6.write('csvFileName = %s \n'%(G.csvFileName))
        else:
            f6.write('lowerBound = %s \n'%(G.dataLowerBound))
            f6.write('higherBound = %s \n'%(G.dataUpperBound))
            f6.write('increments = %s \n'%(G.bruteIncrements))
        f6.write('permutationNum = %s \n'%(G.permutationNum))
    if G.extractReferencePatternFromDataOption == 'yes':
        f6.write('extractReferencePatternFromDataOption = %s \n'%(G.extractReferencePatternFromDataOption))
        f6.write('rpcMoleculesToChange = %s \n'%(G.rpcMoleculesToChange))
        f6.write('rpcMoleculesToChangeMF = %s \n'%(G.rpcMoleculesToChangeMF))
        f6.write('rpcTimeRanges = %s \n'%(G.rpcTimeRanges))
    if G.minimalReferenceValue == 'yes':
        f6.write('minimalReferenceValue = %s \n'%(G.minimalReferenceValue))
        f6.write('referenceValueThreshold = %s \n'%(G.referenceValueThreshold))
        f6.write('referenceSignificantFragmentThresholds = %s \n'%(G.referenceSignificantFragmentThresholds))
    if G.lowerBoundThresholdChooser == 'yes':
        f6.write('lowerBoundThresholdChooser = %s \n'%(G.lowerBoundThresholdChooser))
        f6.write('massesToLowerBoundThresholdFilter  = %s \n'%(G.massesToLowerBoundThresholdFilter ))
        f6.write('lowerBoundThresholdPercentage = %s \n'%(G.lowerBoundThresholdPercentage))
        f6.write('lowerBoundThresholdAbsolute = %s \n'%(G.lowerBoundThresholdAbsolute))
    if G.dataSmootherYorN == 'yes':
        f6.write('dataSmootherYorN = %s \n'%(G.dataSmootherYorN))
        f6.write('dataSmootherChoice = %s \n'%(G.dataSmootherChoice))
        f6.write('dataSmootherTimeRadius = %s \n'%(G.dataSmootherTimeRadius))
        f6.write('dataSmootherPointRadius = %s \n'%(G.dataSmootherPointRadius))
        f6.write('dataSmootherHeadersToConfineTo = %s \n'%(G.dataSmootherHeadersToConfineTo))
    if G.rawSignalThresholdMethod == 'yes':
        f6.write('rawSignalThresholdMethod = %s \n'%(G.rawSignalThresholdMethod))
        f6.write('rawSignalThresholdValue = %s \n'%(G.rawSignalThresholdValue))
        f6.write('sensitivityThresholdValue = %s \n'%(G.sensitivityThresholdValue))
        f6.write('rawSignalThresholdDivider = %s \n'%(G.rawSignalThresholdDivider))
    if G.rawSignalThresholdLimit  == 'yes':
        f6.write('rawSignalThresholdLimit = %s \n'%(G.rawSignalThresholdLimit))
        f6.write('rawSignalThresholdLimitPercent  = %s \n'%(G.rawSignalThresholdLimitPercent))
    if G.negativeAnalyzerYorN == 'yes':
        f6.write('negativeAnalyzerYorN = %s \n'%(G.negativeAnalyzerYorN))
    if G.dataAnalysis == 'yes':
        f6.write('answer = %s \n'%(G.answer))
        if G.answer == 'sls':
            f6.write('uniqueOrCommon = %s \n'%(G.uniqueOrCommon))
            f6.write('slsFinish = %s \n'%(G.slsFinish))
            if G.slsFinish == 'brute':
                f6.write('objectiveFunctionType = %s \n'%(G.objectiveFunctionType))
            if G.slsFinish == 'inverse':
                f6.write('distinguished = %s \n'%(G.distinguished))
        if G.answer == 'inverse':
            f6.write('distinguished = %s \n'%(G.distinguished))
    if G.concentrationFinder == 'yes':
        f6.write('concentrationFinder = %s \n'%(G.concentrationFinder))
        f6.write('molecule = %s \n'%(G.moleculesTSC_List))
        f6.write('moleculeSignal = %s \n'%(G.moleculeSignalTSC_List))
        f6.write('massNumber = %s \n'%(G.massNumberTSC_List))
        f6.write('moleculeConcentration = %s \n'%(G.moleculeConcentrationTSC_List))
        f6.write('units = %s \n'%(G.unitsTSC))
    f6.write('scaleRawDataFactor  = %s \n'%(G.scaleRawDataFactor ))
    f6.write('resolvedScaledConcentrationsOutputName  = %s \n'%(G.resolvedScaledConcentrationsOutputName ))
    f6.write('concentrationsOutputName = %s \n'%(G.concentrationsOutputName))
    f6.write('simulatedSignalsOutputName = %s \n'%(G.simulatedSignalsOutputName))
    f6.write('Run Time %.5f seconds \n'%(G.checkpoint-G.start))
    f6.write('###################################################################### \n')
    f6.close()#once the log file is printed the program is finished
    return None
    
    
def debuggingExportIndividualItem (prefix, objectToExport, objectName = ""):

    prefix = str(prefix)
    objectName = str(objectName)
    print(numpy.shape(objectToExport), type([])); 
    if type(objectToExport) == type([]):
        fileName = "DebugExport" + prefix + objectName + ".csv"
        #use numpy.savetxt here 
        #numpy.savetxt(fileName, prefix, objectName)
        numpy.savetxt(fileName, objectToExport, delimiter = ',', fmt ="%s")
    if type(objectToExport) == type(numpy.array([0])):

        fileName = "DebugExport" + prefix + objectName + ".csv"
        #use numpy.savetxt here 
        #numpy.savetxt(fileName, prefix, objectName)
        numpy.savetxt(fileName, objectToExport, delimiter = ',', fmt ="%s")
    if type(objectToExport) == type(str):

        fileName = "DebugExport" + prefix + objectName + ".txt"
            #use numpy
        f6 = open(fileName,'a')
        f6.write('\n')
        f6.write('prefix = %s \n'%(prefix))
        f6.write('Object  = %s \n'%(objectToExport))
        f6.write('objectName = %s \n'%(objectName))
        f6.close()
    else:

        fileName = "DebugExportWARNING" + prefix + objectName + ".txt"
        print("warning " + prefix + objectName + " is not supported")
        #write warning text to file name 
        f6 = open(fileName,'a')
        f6.write('\n')
        f6.write('prefix = %s \n'%(prefix))
        f6.write('Object  = %s \n'%(objectToExport))
        f6.write('objectName = %s \n'%(objectName))
        f6.close()

def debuggingExportListofItem (prefix, listOfobject = [], listOfObjectNames = ""):

    prefix = str(prefix)
    listOfObjectNames = str(listOfObjectNames)
    
    if type(listOfobject) == type([]):
        fileName = "DebugExport" + prefix + listOfObjectNames + ".csv"
        #use numpy.savetxt here 
        numpy.savetxt(fileName, listOfobject, delimiter = ',', fmt ="%s")
    if type(listOfobject) == type(numpy):

        fileName = "DebugExport" + prefix + listOfObjectNames + ".csv"
        #use numpy.savetxt here 
        numpy.savetxt(fileName, listOfobject, delimiter = ',', fmt ="%s")
    if type(listOfobject) == type(str):

        fileName = "DebugExport" + prefix + listOfObjectNames + ".txt"
            #use numpy
        f6 = open(fileName,'a')
        f6.write('\n')
        f6.write('prefix = %s \n'%(prefix))
        f6.write('Object  = %s \n'%(listOfObjectNames))
        f6.write('listOfObjectNames = %s \n'%(listOfObjectNames))
        f6.close()
    else:

        fileName = "DebugExportWARNING" + prefix + listOfObjectNames + ".txt"
        print("warning" + prefix + listOfObjectNames + "is not supported")
        #write warning text to file name 
        f6 = open(fileName,'a')
        f6.write('\n')
        f6.write('prefix = %s \n'%(prefix))
        f6.write('Object  = %s \n'%(listOfObjectNames))
        f6.write('listOfObjectNames = %s \n'%(listOfObjectNames))
        f6.close()
   
##################################################################################################################
###############################################Algorithm Part 3: Main Control Function ###################################
##################################################################################################################
def main():
    global G #This connects the local variable G to the global variable G, so we can assign the variable G below as needed.    
    G.lastFigureNumber = 0
    
    global debuggingExportIndividualItem; debuggingExportIndividualItem = False #This line turns on and off the "debuggingExportIndividualItem" funciton.
    global debuggingExportListofItem; debuggingExportListofItem = False #This line turns on and off the "debuggingExportListofItem" funciton.
  
    filesAndDirectories = os.listdir()
    for name in filesAndDirectories:
        if name.startswith("Exported") and (name.endswith(".csv") or name.endswith(".tsv") or name.endswith(".skv")):
            print("Previous run Exported file detected. Deleting file", name)
            os.remove(name)
    # #The below try statement is to check the user input dictionary's existence. Older MSRESOLVE did not use a dictionary.
    # for now, these types of lines are at the bottom of the UserInput and DefaultUserInput files. I'm considering keeping there permanently and then deleting these commented out lines.
    #from userInputValidityFunctions import parseUserInput, userInputValidityCheck, settingsCompatibilityCheck, #settingsDependenciesCheck,populateModuleVariablesFromDictionary,populateModuleVariablesFromNestedDictionary
    # try:
        # len(G.UserChoices) #this will fail if the dictionary has not yet been defined.
        # UserChoicesExists = True
    # except:
        # UserChoicesExists = False
    # if UserChoicesExists == True:
        # #if dictionary exists, then code will continue with userinput validation. It is only valid when a UserChoices dictionary exists, which does not exist in 'old' versions of MSRESOLVE (also not used in some unit tests).
        # G.SettingsVDictionary = userInputValidityCheck(G.UserChoices)
        # #NOTE: SettingsVDictionary is created inside userInputValidityCheck(G)
        # populateModuleVariablesFromDictionary(G, G.SettingsVDictionary)
    try:
        len(G.AllMID_ObjectsDict)
    except:
    	#Create the dictionary storing the Molecular Ionization Data objects if the ionization data file exists in the main directory
        #During normal use, this except statement will always be entered except during the iterative analysis unit test.
        #It is actually okay if this except statement runs every time, but it saves a bit of time in the iterative analysis unit test if it does not get run each time and is retained through each iteration.
        #while it would be desirable to put his below the iterative analysis  if statement, with the currentcode flow, it must happen before the directory switch and therefore must be before the iterative code block.
        G.AllMID_ObjectsDict = populateAllMID_ObjectsDict(G.ionizationDataFileName)
    if G.iterativeAnalysis:
        #This section is to overwrite the UI if iterative analysis is in the process of being run. 
        highestIteration = int(FindHighestDirNumber("_iter_"))
        iterationDirectorySuffix = '_iter_%s' %str(highestIteration)
        for directoryName in filesAndDirectories:
            if iterationDirectorySuffix in directoryName:
                userInputName = 'UserInput{}'.format(iterationDirectorySuffix)
                userInputPath = '{}.{}'.format(directoryName, userInputName)
                UserInputCurrentIteration = importlib.import_module(str(userInputPath))
                AllMID_ObjectsDict = G.AllMID_ObjectsDict #Temporarily store the global variable as a local variable
                G = UserInputCurrentIteration
                G.AllMID_ObjectsDict = AllMID_ObjectsDict #Put the local variable back into the namespace after G points to the new namespace.  This retains the MID dictionary through different iterations
                G.lastFigureNumber = 0
                break
        if G.iterativeAnalysis:
            G.iterationNumber = highestIteration
            G.iterationSuffix = iterationDirectorySuffix
        elif not G.iterativeAnalysis:
            G.iterationSuffix = ''
    
    try:
        type(G.calculateUncertaintiesInConcentrations) #Need to check if this variable exists. If it does not, then fill it with None.
    except:
        G.calculateUncertaintiesInConcentrations = None

    #if this is not the first iterative run, then the required files are all stored in the highest iteration directory
    if G.iterativeAnalysis and G.iterationNumber != 1:
        #implied arguments for this function are G.referenceFileNamesList and G.collectedFileName
        IterationDirectoryPreparation(G.iterativeAnalysis, G.iterationNumber) #This function also changes the working directory

    #Read in the molecules used before parsing the user input file    
    G.referenceFileNamesList = parse.listCast(G.referenceFileNamesList)
    G.referenceFormsList = parse.listCast(G.referenceFormsList)
    
    [provided_reference_patterns, electronnumbers, molecules, molecularWeights, 
        SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType,
        mass_fragment_numbers_monitored, referenceFileName, form]=readReferenceFile(G.referenceFileNamesList[0], G.referenceFormsList[0])
    G.moleculesNames = molecules

    
    #If a tuning correction is going to be done, and if requested, we'll make a mixed reference pattern.
    #G.moleculesNamesExtended needs to be populated before the first tuning correction and before parsing of userinput.
    if str(G.measuredReferenceYorN).lower() == 'yes': 
        if G.createMixedTuningPattern == True:
            if len(G.referenceFileExistingTuningAndForm) == 0:
                G.referenceFileExistingTuningAndForm = G.referenceFileStandardTuningAndForm #Use the standard tuning file if blank.
            [provided_reference_patterns, electronnumbers, molecules, molecularWeights, 
            SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType,
            mass_fragment_numbers_monitored, referenceFileName, form] = readReferenceFile(G.referenceFileExistingTuningAndForm[0], G.referenceFileExistingTuningAndForm[1])
            G.moleculesNamesExistingTuning = molecules
            moleculesToAddToReferencePattern = []
            for moleculeName in list(G.moleculesNamesExistingTuning):
                if moleculeName in list(G.moleculesNames):
                    pass
                else:
                    moleculesToAddToReferencePattern.append(moleculeName)
            G.moleculesNamesExtended = list(G.moleculesNames) + list(moleculesToAddToReferencePattern)
        #We are reading the experimental data in and this must be before user input processing so we have the mass fragments
    
    
    G.exp_mass_fragment_numbers = getMassFragmentsFromCollectedData(G.collectedFileName)

    #Save an MSReference object containing all molecules and an MSData object containing all mass fragments
    if G.iterativeAnalysis and G.iterationNumber != 1: #If using iterative and not on the first iteration we will need to remove _iter_x from the file names
        AllMoleculesReferenceFileNamesList = [] #Initialize AllMoleculesReferenceDataList as an empty list
        for referenceFileNameIndex in range(len(G.referenceFileNamesList)): #Loop through the reference file names list
            AllMoleculesReferenceFileName = remove_iter_fromFileName(G.referenceFileNamesList[referenceFileNameIndex]) #Remove the _iter_ from the name so the program has the original filename to access from the parent directory
            AllMoleculesReferenceDataFilePath = os.path.normpath(os.path.join(os.curdir, os.pardir,AllMoleculesReferenceFileName)) #This function will get the path of the reference file from the parent directory 
            AllMoleculesReferenceFileNamesList.append(AllMoleculesReferenceDataFilePath) #Append the path to the list and the program will read the reference file from the path name
        AllMassFragmentsExperimentDataFileName = remove_iter_fromFileName(G.collectedFileName) #Remove _iter_ from the data filename so the program has the original filename to access from the parent directory
        AllMassFragmentsExperimentDataFileNamePath = os.path.normpath(os.path.join(os.curdir, os.pardir, AllMassFragmentsExperimentDataFileName)) #This function will get the path of the data file from the parent directory
    else: #Otherwise not running iterative or in the first iteration, just copy the filename
        AllMoleculesReferenceFileNamesList = copy.copy(G.referenceFileNamesList)
        AllMassFragmentsExperimentDataFileNamePath = copy.copy(G.collectedFileName)
    #Create the MSReference and MSData objects containing all molecules and all mass fragments, respectively
    [exp_mass_fragment_numbers, exp_abscissaHeader, exp_times, exp_rawCollectedData, exp_collectedFileName]=readDataFile(AllMassFragmentsExperimentDataFileNamePath)
    AllMassFragmentsExperimentData = MSData(exp_mass_fragment_numbers, exp_abscissaHeader, exp_times, exp_rawCollectedData, collectedFileName=exp_collectedFileName)        
    AllMoleculesReferenceDataList = GenerateReferenceDataList(AllMoleculesReferenceFileNamesList,G.referenceFormsList,G.AllMID_ObjectsDict)
    #Then prepare AllMoleculesReferenceDataList to get reciprocal_matching_correction_values, this value is fed into RatioFinder
    global PreparingAllMoleculesReferenceDataList; PreparingAllMoleculesReferenceDataList = False #initializing this flag          
    for referenceObjectIndex in range(len(AllMoleculesReferenceDataList)):
        PreparingAllMoleculesReferenceDataList = True #Set flag to true before doing each preparation.  
        AllMoleculesReferenceDataList[referenceObjectIndex].ExportAtEachStep = 'no'
        PrepareReferenceObjectsAndCorrectionValues(AllMoleculesReferenceDataList[referenceObjectIndex],AllMassFragmentsExperimentData.mass_fragment_numbers, AllMassFragmentsExperimentData)
        PreparingAllMoleculesReferenceDataList = False #Set flag to false after doing each preparation.                                 
    #beforeParsedGDict this will be needed for iterative. This actually contains "UserChoices" when that's available, but we won't use that.
    beforeParsedGDict = {}
    if len(G.__dict__.items())>0:       
        for (key,value) in G.__dict__.items():  #This gets a tuple of all of the items in G.
            if type(value) != type(copy): #don't want to do this with use modules.
                if not key.startswith("__"): #We don't want the ones that start with "__" since those are not normal variables.
                    try: #We're going to take deep copies, but some of the items are still not normal (like modules) so can be non-copyable and thus we have an except statement below.
                        beforeParsedGDict[key] = copy.deepcopy(value)  
                    except TypeError:
                        pass
    G.beforeParsedGDict = beforeParsedGDict # now we save it into the globals module for using later.
    from userInputValidityFunctions import parseUserInput
    parseUserInput(G) #This parses the variables in the user input file
            
    #it is useful to trim whitespace from each chosenMolecules string. The same thing is done to the molecule names of each reference pattern when an MSReference object is created.
    for moleculeIndex, moleculeName in enumerate(G.chosenMoleculesNames):
        G.chosenMoleculesNames[moleculeIndex] = moleculeName.strip()
    
    #Record the time
    G.start = timeit.default_timer()
    G.checkpoint = timeit.default_timer()

    #initalize the data classes with the data from given Excel files
    #These are being made into globals primarily for unit testing and that functions are expected to receive the data as arguments rather than accessing them as globals
    global ReferenceDataList
    global ExperimentData
    global prototypicalReferenceData
    global currentReferenceData
    global resultsObjects
    resultsObjects = {}
    [exp_mass_fragment_numbers, exp_abscissaHeader, exp_times, exp_rawCollectedData, exp_collectedFileName]=readDataFile(G.collectedFileName)
    ExperimentData = MSData(exp_mass_fragment_numbers, exp_abscissaHeader, exp_times, exp_rawCollectedData, collectedFileName=exp_collectedFileName)
    ReferenceDataList = GenerateReferenceDataList(G.referenceFileNamesList,G.referenceFormsList,G.AllMID_ObjectsDict)
    ExperimentData.provided_mass_fragment_numbers = ExperimentData.mass_fragment_numbers
    #This is where the experimental uncertainties object first gets populated, but it does get modified later as masses are removed and time-points are removed.
    if type(G.collectedFileUncertainties) != type(None):
        if type(G.collectedFileUncertainties) != type("str") and type(G.collectedFileUncertainties) != type(['list']) : # if it's not a string or a list, then it's assumed to be an integer for which local region to use during/like datasmoother.
            G.collectedFileUncertainties = int(G.collectedFileUncertainties)

    #Prints a warning if the user has more reference files than specified time ranges
    if len(G.referencePatternTimeRanges) > 0 and (len(G.referenceFileNamesList) > len(G.referencePatternTimeRanges)):
        print("WARNING: There are more reference files given than time ranges")
    #save global variable into the class objects 
    ExperimentData.ExportAtEachStep = G.ExportAtEachStep
   
    #if this is the first iterative run, then the reference and experimental files need to have been imported before the iteration can begin
    if G.iterativeAnalysis and G.iterationNumber == 1 :
        #implied arguments for the following function are G.referenceFileNamesList and G.collectedFileName
        IterationFirstDirectoryPreparation(G.iterativeAnalysis, G.iterationNumber)

    # Skip preProcessing all together if we are loading analyzed data
    if(G.dataAnalysis == 'load'):
        print("DataAnalysis set to 'load': skipping preprocessing")

        # Skip preprocessing
        G.preProcessing = 'skip'
        
    if(G.preProcessing == 'yes'):
        
        #TODO Make a new global to remove mass fragments from the experimental data in preprocessing
        #using trimDataMassesToMatchChosenMassFragments
        
        # Perform the actual data preprocessing on ExperimentData
        ExperimentData = DataInputPreProcessing(ExperimentData)
        print("Data PreProcessing completed")
        #This graph call is graphing fully preprocessed data.
        if G.grapher == 'yes':
            print('PreProcessed Data Graph')
            Draw(ExperimentData.times, ExperimentData.workingData, ExperimentData.mass_fragment_numbers, 'no', 'Amp', graphFileName = 'PreprocessingAfterSmoothing', fileSuffix = G.iterationSuffix, label="PreProcessed Data Graph", stopAtGraphs=G.stopAtGraphs, figureNumber=G.lastFigureNumber+1 )
            G.lastFigureNumber = G.lastFigureNumber+1

        #Exports the Preprocessed Data
        ExportXYYYData(G.preProcessedDataOutputName, ExperimentData.workingData, ExperimentData.mass_fragment_numbers, 
                       abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = 'preProcessed', rowIndex = ExperimentData.times)
        print("Preprocessed data exported")
        
        #Export collected export data
        ExperimentData.ExportMSData()
        
        #show net time for preProcessing
        #Record time in case data analysis isn't used
        G.timeSinceLastCheckPoint = timeit.default_timer() - G.checkpoint
        G.checkpoint = timeit.default_timer()
        print('PreProcessing Time: ', (G.timeSinceLastCheckPoint))
    
    elif(G.preProcessing == 'skip'):
        
        # Output to make sure user knows we are skipping Preprocessing
        G.timeSinceLastCheckPoint = timeit.default_timer() - G.checkpoint
        G.checkpoint = timeit.default_timer()
        print("Preprocessing set to 'skip'. Performed only mandatory PreProcessing with Time: {:.4f}".format(G.timeSinceLastCheckPoint))
            
    elif(G.preProcessing == 'load' and G.dataAnalysis == 'yes'):
            
        #This function call loads the preprocessed data
        print("Loading the preprocessed data from file '{}'".format(G.preProcessedDataOutputName))
        ExperimentData.workingData, ExperimentData.mass_fragment_numbers, ExperimentData.times = ImportWorkingData(G.preProcessedDataOutputName)
             
        # Output to make sure user knows we are loading Preprocessing
        G.timeSinceLastCheckPoint = timeit.default_timer() - G.checkpoint
        G.checkpoint = timeit.default_timer()
        print("Preprocessed data loaded. Time: {:.4f}".format(G.timeSinceLastCheckPoint))

    else:
        # if we are here then 'G.preProcessing' != ('yes' or 'skip' or 'load')
        raise ValueError("The value of preProcessing is not set appropriately, it should be 'yes', 'skip' or 'load'." +
                         "Or you are attempting to load pre-processed data without running data analysis")

    # Reference Pattern Changer
    if G.extractReferencePatternFromDataOption == 'yes':
        for i in range(len(ReferenceDataList)):#perform for each reference pattenr.
            ReferenceDataList[i].provided_reference_patterns, ReferenceDataList[i].provided_reference_patterns_absolute_uncertainties = ExtractReferencePatternFromData(ExperimentData, ReferenceDataList[i], G.rpcMoleculesToChange,G.rpcMoleculesToChangeMF,G.rpcTimeRanges)
            ReferenceDataList[i].ExportCollector('ExtractReferencePatternFromData',use_provided_reference_patterns = True)
            ReferenceDataList[i].ExportCollector('ExtractReferencePatternFromData_absolute_uncertainties', export_uncertainties= True)
        print('ReferencePatternChanger complete')

    #This codeblock is for the TuningCorrectorGasMixture feature. It should be before the prototypicalReferenceData is created.
    #A measured gas mixture spectrum is compared to a simulated gas mixture spectrum, and the tuning correction is then made accordingly.
    if G.measuredReferenceYorN == 'yes':
        if len(G.UserChoices['measuredReferenceYorN']['tuningCorrectorGasMixtureMoleculeNames']) > 0:
            ReferenceDataList = tuningCorrectorGasMixture(ReferenceDataList, G, ExperimentData)
    #Creating prototypicalReferenceData which will be interrogated later for which molecules and masses to expect in the ReferenceDataObjects.
    prototypicalReferenceData = copy.deepcopy(ReferenceDataList[0])
    #TODO make a variable allMoleculesAnalyzed that is a list containing all the molecules analyzed so far
    #Steps required if preprocessing has also been run
    if (G.dataAnalysis == 'yes' and G.preProcessing == 'yes'):
        
        if G.iterativeAnalysis:
            ReferenceDataListFullCopy = []
            for RefObjectIndex, RefObject in enumerate(ReferenceDataList): #a list
                #create a copy of the Reference Data
                ReferenceDataListFullCopy.append(copy.deepcopy(RefObject)) 
        ##Start: Preparing data for data analysis based on user input choices
        # Trim the reference data according to the selected molecules list
        #TODO: Ashi thinks this can be moved into PrepareReferenceObjectsAndCorrectionValues. But care of not messing up iterativeAnalysis would be needed.
        #TODO continued: In that case, we could change the trimming function to work on the standardized reference patterns.
        #TODO continued: This would leave provided_reference_intensities unchanged, which is an additional bonus.
        #TODO continued: If we parse chosenMolecules appropriately ahead of time, we can just pass that as an argument to prepare function.
        #TODO continued: If it matches (as a set, not an equal comparison of lists), then no trimming occurs. Else, trimming occurs.
        #TODO continued: but the above changes would probably also require trimDataMassesToMatchReference to occur on Experimental data, again, after prepare reference patterns, including on ExperimentDataCopy.
        if G.specificMolecules == 'yes' or G.iterativeAnalysis:
            if G.specificMolecules == 'no': #This should really be in user input validation.
                print("WARNING: You have chosenMolecules / specificMolecules set to no, but iterativeAnalysis is set to True and requires chosenMolecules. ChosenMolecules will be used as part of iterative. If you did not fill out the chosenMolecules list correctly, you must do so and run again (you may also need to delete the directory for the next iteration).")
            for RefObjectIndex, RefObject in enumerate(ReferenceDataList): #a list
                ReferenceDataList[RefObjectIndex] = trimDataMoleculesToMatchChosenMolecules(RefObject, G.chosenMoleculesNames)
                #ReferenceDataList[RefObjectIndex].ClearZeroRowsFromProvidedReferenceIntensities()
                if hasattr(ReferenceDataList[RefObjectIndex], 'standardized_reference_patterns'):
                    pass#ReferenceDataList[RefObjectIndex].ClearZeroRowsFromStandardizedReferenceIntensities()
            prototypicalReferenceData = trimDataMoleculesToMatchChosenMolecules(prototypicalReferenceData, G.chosenMoleculesNames)
            #prototypicalReferenceData.ClearZeroRowsFromProvidedReferenceIntensities()
            if hasattr(prototypicalReferenceData, 'standardized_reference_patterns'):
                pass#prototypicalReferenceData.ClearZeroRowsFromStandardizedReferenceIntensities()
	
        if G.iterativeAnalysis:
            #make a copy of the experimental data for later use in iterative processing
            ExperimentDataFullCopy = copy.deepcopy(ExperimentData)
            #make a copy of Experimental data specifically to be used in signal simulation. i.e. will have mass fragments trimmed if they aren't referenced by the current molecules. 
            ExperimentDataCopy = copy.deepcopy(ExperimentData)
            #remove any unreference masses from the signal simulation copy of experimental data
            ExperimentDataCopy = trimDataMassesToMatchReference(ExperimentDataCopy, prototypicalReferenceData)       
        
        
        # If we are only interested in a subset of the MS data
        # remove the irrelevant mass data series from ExperimentData.mass_fragment_numbers
        # and the corresponding colums from ExperimentData.workingData
        if G.specificMassFragments == 'yes':
            print("MassFragChooser")
             #Trim the experimental data according to the mass fragments in G.chosenMassFragments 
            ExperimentData = trimDataMassesToMatchChosenMassFragments(ExperimentData, G.chosenMassFragments) 
            ExperimentData.ExportCollector("MassFragChooser")
        #Trim the experimental data according to the mass fragments in referenceData
        ExperimentData = trimDataMassesToMatchReference(ExperimentData, prototypicalReferenceData) 

    ## Here perform the ReferenceData preprocessing that is required regardless of the selection for 'G.preProcessing'
    # and needed if G.dataAnalysis == 'load' or 'yes'  
    if (G.dataAnalysis == 'yes' or G.dataAnalysis =='load'):
        #Prepare prototypicalReferenceData which is currently the first reference object in the list
        extractReferencePatternFromDataOptionHere = 'no' #should be extracted already.
        prototypicalReferenceData = PrepareReferenceObjectsAndCorrectionValues(prototypicalReferenceData, ExperimentData.mass_fragment_numbers, ExperimentData,extractReferencePatternFromDataOptionHere, G.rpcMoleculesToChange,G.rpcMoleculesToChangeMF,G.rpcTimeRanges)
        #for loop to preprocess the remaining MSReference objects and match correction values
        for i in range(len(ReferenceDataList)):
            try:
                ReferenceDataList[i].populateIonizationEfficiencies(G.AllMID_ObjectsDict)
            except:
                ReferenceDataList[i].populateIonizationEfficiencies()
            extractReferencePatternFromDataOptionHere = 'no' #should be extracted already
            ReferenceDataList[i] = PrepareReferenceObjectsAndCorrectionValues(ReferenceDataList[i],ExperimentData.mass_fragment_numbers, ExperimentData, extractReferencePatternFromDataOptionHere, G.rpcMoleculesToChange,G.rpcMoleculesToChangeMF,G.rpcTimeRanges)
            
    if (G.dataAnalysis == 'yes'):
        print("Entering Data Analysis")
        #The iterative analysis preprocessing creates the proper export folder and exports the unused reference data
        if G.iterativeAnalysis:
            ReferenceDataSSreciprocal_matching_correction_valuesList, G.unusedMolecules = IADirandVarPopulation(G.iterativeAnalysis, G.chosenMassFragments, G.chosenMoleculesNames, ExperimentData, ExperimentDataFullCopy, ReferenceDataList, ReferenceDataListFullCopy)
                
        # Reset the checkpoint timer for the data analysis section
        G.checkpoint = timeit.default_timer()
        #check to make sure that there are enough mass fragments to solve for each variable. 
        if len(ExperimentData.mass_fragment_numbers) < len(ReferenceDataList[0].molecules):
            raise SystemError("There are too few mass fragments to solve for the number of molecules provided. \n Data Analysis has been ended.")
        
        # Since SlSUniqueFragments is potentially used in a number of the analysis options
        # set up the 'SLSUniqueOrder.csv' file headers here
        # note that this open overwrites previous file contents
        # while subsequent opens to this file will append
        if G.SLSUniqueExport == 'yes':
            createSLSUniqueExportFiles(ExperimentData.abscissaHeader,
                                     prototypicalReferenceData.molecules)
            
        #this numpy.zeros line is going to be the array that holds all of the answers before they are printed out, which
        #is done in order to save time and decrease expense
        concentrationsScaledToCOarray = numpy.zeros(len(prototypicalReferenceData.molecules)+1)
        concentrationsarray = numpy.zeros(len(prototypicalReferenceData.molecules)+1)
        correctionFactorArraysList = [] #initialize the correctionFactorArray as an empty list
        SS_reciprocal_matching_correction_values_TimesList = [] #initialze ReferenceDataSSreciprocal_matching_correction_valuesList as an empty list.  This list will store correction values used at each time point
        
        # Loading user choices for data analysis
        DataRangeSpecifierlist = [G.dataRangeSpecifierYorN, G.signalOrConcentrationRange,
                                  G.csvFile, G.moleculesToRestrict, G.csvFileName,G.dataUpperBound,
                                  G.dataLowerBound, G.bruteIncrements, G.permutationNum]
        SLSChoices = [G.uniqueOrCommon, G.slsFinish, G.distinguished]
        ThresholdList = [G.rawSignalThresholdMethod, G.rawSignalThresholdValue, G.sensitivityThresholdValue,
                         G.rawSignalThresholdDivider, G.rawSignalThresholdLimit, G.rawSignalThresholdLimitPercent]
        
        currentReferenceData = ReferenceDataList[0] #TODO this line is placeholder by charles to fix currentRefenceData issue until Alex has a better solution 
    
        # Calculate a coefficient for doing a unit conversion on concentrations #TODO resolve Ratio Finder issue, i.e. list of conversionValues
        conversionFactorsAtEachTime = RatioFinder(AllMoleculesReferenceDataList, AllMassFragmentsExperimentData, ReferenceDataList, ExperimentData, G.concentrationFinder, G.TSC_List_Type,
                                      G.moleculesTSC_List, G.moleculeConcentrationTSC_List, G.massNumberTSC_List, G.moleculeSignalTSC_List,G.referencePatternTimeRanges)
	##End: Preparing data for data analysis based on user input choices

        #Initialize a current reference pattern index
        currentReferencePatternIndex = 0
        numberOfTimePoints = len(ExperimentData.workingData[:,0])
        for timeIndex in range(numberOfTimePoints):#the loop that runs the program to get a set of signals/concentrations for each time  
            #These print statements are used to track the progress of the program during long analysis runs
            if (timeIndex == 0):
                print("TimeIndex", timeIndex, "\n Estimated runtime remaining:", "unknown", "seconds.") 
                lastPrintedIndexClockTime = timeit.default_timer() 
            if ((timeIndex % 100) == 0 and timeIndex != 0):
                averageRunTimePerPoint = (timeit.default_timer() - G.checkpoint)/timeIndex
                estimatedRemainingRuntime = (numberOfTimePoints-timeIndex)*averageRunTimePerPoint
                print("TimeIndex", timeIndex, "\n Estimated runtime remaining:", estimatedRemainingRuntime, "seconds.")
                lastPrintedIndexClockTime = timeit.default_timer() 
            if timeit.default_timer() - lastPrintedIndexClockTime > 10:
                averageRunTimePerPoint = (timeit.default_timer() - G.checkpoint)/timeIndex
                estimatedRemainingRuntime = (numberOfTimePoints-timeIndex)*averageRunTimePerPoint
                print("TimeIndex", timeIndex, ": Ten seconds or more checkpoint. Actually", timeit.default_timer() - lastPrintedIndexClockTime, "seconds. \n Estimated runtime remaining:", estimatedRemainingRuntime, "seconds.")
                lastPrintedIndexClockTime = timeit.default_timer() #reset this value.
            #If referencePatternTimeRanges has anything in it, then the user has opted to use the Reference Pattern Time Chooser feature
            #FIXME interpolation portion of reference pattern time chooser does not work yet
            if len(G.referencePatternTimeRanges) != 0:
                #If we are in the last time range, calling this function will result in an index error
                #If using this feature, (len(G.referencePatternTimeRanges)) will always be at least 2 time ranges so use len(G.referencePatternTimeRanges)-1
                if currentReferencePatternIndex < (len(G.referencePatternTimeRanges)-1):    
                    currentReferenceData, currentReferencePatternIndex = SelectReferencePattern(currentReferencePatternIndex, G.referencePatternTimeRanges, ExperimentData.times[timeIndex], ReferenceDataList[currentReferencePatternIndex], ReferenceDataList[currentReferencePatternIndex+1], ReferenceDataList)
            else: #referencePatternTimeRanges is empty so user is opting to not use reference pattern time chooser
                currentReferenceData = ReferenceDataList[0]
            
            #populate the mass fragments monitored subobject for the current reference pattern
            currentReferenceData.mass_fragment_numbers_monitored = ExperimentData.mass_fragment_numbers
            
            #populate the uncertainties_dict now that the reference pattern has been chosen. First make it empty, then add more if warranted.
            uncertainties_dict = {}
            if (G.calculateUncertaintiesInConcentrations == True):
                if type(G.referenceFileUncertainties) != type(None):  #This means we should have uncertainties in the reference file.
                    uncertainties_dict['correction_values_relative_uncertainties_one_time'] = currentReferenceData.correction_values_relative_uncertainties*1.0
                    uncertainties_dict['reciprocal_matching_correction_values_relative_uncertainties_one_time'] = currentReferenceData.reciprocal_matching_correction_values_relative_uncertainties*1.0
                elif type(G.referenceFileUncertainties) == type(None):
                    uncertainties_dict['correction_values_relative_uncertainties'] = None
                #G.collectedFileUncertainties = None
                if (type(G.collectedFileUncertainties) != type(None)) and (str(G.collectedFileUncertainties).lower() != "none"):
                    uncertainties_dict['rawsignals_absolute_uncertainties'] = ExperimentData.rawsignals_absolute_uncertainties #This is for *all* times.
                    uncertainties_dict['rawsignals_absolute_uncertainties_one_time']=uncertainties_dict['rawsignals_absolute_uncertainties'][timeIndex] #This is for one time.                                                                                                                                         
                elif type(G.collectedFileUncertainties) == type(None):
                    uncertainties_dict['rawsignals_absolute_uncertainties'] = None
                    uncertainties_dict['rawsignals_absolute_uncertainties_one_time']= ExperimentData.workingData[timeIndex]*0.0                                                                                                             
                uncertainties_dict['concentrations_relative_uncertainties_one_time'] = numpy.zeros(len(currentReferenceData.molecules))
                uncertainties_dict['concentrations_absolute_uncertainties_one_time'] = numpy.zeros(len(currentReferenceData.molecules))
                    
                # referenceCorrectionCoefficientsUncertainties does not matter because it is used to get correction_values_relative_uncertainties.
                #referenceCorrectionCoefficientsIonizationUncertainties does not matter because it is used to get correction_values_relative_uncertainties.
                
            G.massesUsedInSolvingMoleculesForThisPoint = numpy.zeros(len(currentReferenceData.molecules)) #initializing this for use inside SlSUnique and possibly elsewhere.
            G.massFragmentsArrayForThisPoint = currentReferenceData.mass_fragment_numbers_monitored*1.0
            
            ## TODO: Find out why RawSignalsArrayMaker() takes longer to run when preprocessed data is
            # computed directly rather than loaded. It doesn't seem to effect rawsignalsarrayline in
            # either case so not a priority. 
            #TODO : these are not actually raw signals, they are preprocessed to be relevant masses from those which mass.
            rawsignalsarrayline = RawSignalsArrayMaker(currentReferenceData.mass_fragment_numbers_monitored,
                                                       ExperimentData.mass_fragment_numbers,ExperimentData.workingData,
                                                       timeIndex,currentReferenceData.referenceabscissa)#gets the collected values that will be solved
            
            #TODO: rename to excludeMoleculesIfMajorFragmentNotObserved, and call the other valuables things like minimumSignalRequired, and minimumStandardizedReferenceHeightToBeSignificant
            #TODO continued: I am putting some variables here to make that process easier by getting some of it done already, then only the user input file needs to be changed.
            #This feature is intended to remove molecules that have major fragments not observed. previously, it was done in a more complicated way.
            # now, to simplify things, is being used as a filter that simply sets standardized intensities in the reference patterns to zero.
            G.excludeMoleculesIfSignificantFragmentNotObserved = G.rawSignalThresholdMethod
            G.minimumSignalRequired = G.rawSignalThresholdValue 
            G.minimumStandardizedReferenceHeightToBeSignificant = G.sensitivityThresholdValue
            if G.excludeMoleculesIfSignificantFragmentNotObserved == 'yes':
                #FIXME: Check if this is compatible with uncertainties. I don't think it is, as of Feb 2nd 2020.                                                                                                
                currentReferenceData = signalThresholdFilter(currentReferenceData, rawsignalsarrayline, G.minimumSignalRequired, G.minimumStandardizedReferenceHeightToBeSignificant)
                    
                    #solutions =RawSignalThresholdFilter(G.distinguished, currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline,
                                                         # currentReferenceData.monitored_reference_intensities,currentReferenceData.molecules,timeIndex,ExperimentData.mass_fragment_numbers,
                                                         # ThresholdList,G.answer,ExperimentData.times[timeIndex],ExperimentData.conversionfactor,ExperimentData.datafromcsv,
                                                         # DataRangeSpecifierlist,SLSChoices,G.permutationNum,concentrationsScaledToCOarray,G.objectiveFunctionType, G.maxPermutations)           
            if G.answer == 'inverse':#user input, the inverse method
                if G.distinguished == 'yes':#user input, choosing between distinguished inverse method or combinations method
                    solutions = InverseMethodDistinguished(currentReferenceData.monitored_reference_intensities,currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline, uncertainties_dict)
                else:
                    solutions = InverseMethod(currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline,currentReferenceData.monitored_reference_intensities,ExperimentData.mass_fragment_numbers,currentReferenceData.molecules,'composition', uncertainties_dict)

            elif G.answer == 'sls' or G.answer == 'autosolver':#user input, the SLS method is chosen)
                #FIXME: should be passing uncertainties_dict into SLSMethod, but when I try I get some kind of solutions error. Comes out as []. Since I am having trouble, I am just passing it through globals.
                G.uncertainties_dict = uncertainties_dict
                solutions = SLSMethod(currentReferenceData.molecules,currentReferenceData.monitored_reference_intensities,currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline, timeIndex, conversionFactorsAtEachTime, ExperimentData.datafromcsv,currentReferenceData.molecules,DataRangeSpecifierlist,SLSChoices,ExperimentData.mass_fragment_numbers,G.permutationNum,concentrationsScaledToCOarray,G.objectiveFunctionType,ExperimentData.times[timeIndex],G.maxPermutations)
                if G.answer == 'autosolver':
                    if solutions.any() == None:
                        if SLSChoices[0] == "unique":
                            print("SLS Unique has failed, trying SLS common.")
                            solutions = SLSMethod(currentReferenceData.molecules,currentReferenceData.monitored_reference_intensities,currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline, timeIndex, conversionFactorsAtEachTime, ExperimentData.datafromcsv,currentReferenceData.molecules,DataRangeSpecifierlist,["common", G.slsFinish, G.distinguished],ExperimentData.mass_fragment_numbers,G.permutationNum,concentrationsScaledToCOarray,G.objectiveFunctionType,ExperimentData.times[timeIndex],G.maxPermutations)
                            if solutions.any() == None:
                                print("The SLS method failed to solve this problem even with SLS common. Attempting an inverse method solution. To solve without inverse, consider raising the Reference Mass Fragmentation Threshold.")
                                if G.distinguished == 'yes':#user input, choosing between distinguished inverse method or combinations method
                                    solutions = InverseMethodDistinguished(currentReferenceData.monitored_reference_intensities,currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline)
                                else:
                                    solutions = InverseMethod(currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline,currentReferenceData.monitored_reference_intensities,ExperimentData.mass_fragment_numbers,currentReferenceData.molecules,'composition')
            solutions = solutions.flatten() #In the major numpy update before summer 2019, this line suddenly started becoming filled with nested arrays sometimes (not sure why) but flattening fixed it.
            arrayline = []
            for moleculecounter in range(len(currentReferenceData.molecules)):#array-indexed for loop, this is the same data structure as the inverse method above once the line of solutions is found, see above for comments
                if moleculecounter == 0:#only for the first loop will times be added to new collected data
                    arrayline.append(ExperimentData.times[timeIndex])
                arrayline.append(solutions[moleculecounter])
            arrayline = numpy.array(arrayline)
#            if G.fullBrute == 'yes':
#                specifications = 
#                arrayline = OptimizingFinisher(ReferenceData.molecules,  G.objectiveFunctionType)
            try:
                if G.UserChoices['dataAnalysisMethods']['finalOptimization'] != 'None':
                    optimizer = G.UserChoices['dataAnalysisMethods']['finalOptimization']
                    specifications=[solutions]
                    concentrationsFromFinisher = OptimizingFinisher(currentReferenceData.molecules,specifications,currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline,G.UserChoices['dataAnalysisMethods']['objectiveFunctionType'],G.maxPermutations, optimizer=optimizer)
                    solutions = concentrationsFromFinisher
            except:
                pass

            if G.SLSUniqueExport == 'yes':
                outputMassesUsedInSolvingMoleculesFilename = 'ExportedSLSUniqueMassesUsedInSolvingMolecules.csv'
                if G.iterativeAnalysis:
                    #then the filename will have a suffix attached
                    outputMoleculesOrderFileName = outputMassesUsedInSolvingMoleculesFilename[:-4] + '_iter_%s' %G.iterationNumber + outputMassesUsedInSolvingMoleculesFilename[-4:] 
                with open(outputMassesUsedInSolvingMoleculesFilename, 'a') as f:
                    f.write('%s,' %timeIndex)
                    f.write('%s,' %ExperimentData.times[timeIndex]) 
                    f.write( str(list(G.massesUsedInSolvingMoleculesForThisPoint))[1:-1] + "\n" )                    
            
            if G.negativeAnalyzerYorN == 'yes':
                arrayline = NegativeAnalyzer(arrayline,currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline,currentReferenceData.molecules,G.objectiveFunctionType, G.maxPermutations, G.NegativeAnalyzerTopNContributors, G.NegativeAnalyzerBaseNumberOfGridIntervals)

            if G.implicitSLScorrection == True: #Note: Not compatibile with iterative. Needs to be executed before the concentrations (arrayline) are stacked.
            #TODO: May 17 2020. Right now, we only take the LARGEST correction for each sls when deciding which sls mass to correct, then we apply that correction. But maybe we should take a single recursion of all molecules affecting? Then apply all molecules to that sls mass before moving to the next one?
                G.implicitSLScorrectionOccurred = False #This is a flag that will get set to true after the correction.
                if G.iterativeAnalysis == False:
                    if G.uniqueOrCommon=='unique' and (G.answer == 'sls' or G.answer =='autosolver'):
                        solvedConcentrationsAtThisTime = arrayline
                        referenceBeforeFiltering = G.currentReferenceDataUnfiltered
                        referenceAfterFiltering = currentReferenceData
                        slsSolvedMasses = G.massesUsedInSolvingMoleculesForThisPoint
                        if 0 in slsSolvedMasses:
                            print("The concentrations are currently not completely solvable by SLS at this time so implicitSLScorrection is being skipped and changed to False. One option to consider is increasing UserChoices['minimalReferenceValue']['referenceValueThreshold'] to a highernumber." )
                            G.implicitSLScorrection = False
                        else:
                            #Now need to move the class out of this area. Currently test 2 is taking 0.08 to 0.20 anlaysis time versus 0.014 for test_1.
                            filterCorrectingObject = referenceThresholdFilterCorrectingSandbox(solvedConcentrationsAtThisTime,referenceBeforeFiltering,referenceAfterFiltering, slsSolvedMasses)
                            implicitSLScorrectionAmounts = filterCorrectingObject.solvedConcentrationsCorrected1D[1:] - arrayline[1:]  #The 1: is to remove the "time" value.
                            arrayline = filterCorrectingObject.solvedConcentrationsCorrected1D
                            G.implicitSLScorrectionOccurred = True
            if timeIndex == 0: #Can't vstack with an array of zeros or else the first time is 0 with all data points at 0 so make first row the first arrayline provided
                concentrationsScaledToCOarray = arrayline*1.0
            elif timeIndex > 0: #Everything else is appended via numpy.vstack
                concentrationsScaledToCOarray = numpy.vstack((concentrationsScaledToCOarray,arrayline))
            correctionFactorArraysList.append(currentReferenceData.reciprocal_matching_correction_values) #populate the list with the proper correction values
            if (G.calculateUncertaintiesInConcentrations == True):  #Keep track of uncertainties same way as concentrationsScaledToCOarray 
                #We have an 'optional'/'added' correction here related to G.implicitSLScorrection
                if G.implicitSLScorrection == True: 
                    if G.implicitSLScorrectionOccurred == True: #if implicitSLScorrectionOccurred then need to add half of implicitSLScorrectionAmounts to the absolute uncertainties (assume orthogonal). We use half because that would make the full correction approximated as a 95% confidence interval and half the correction as 1 sigma, which is how our errors are defined.
                        implicitSLScorrectionAmounts1sigma = implicitSLScorrectionAmounts*0.5
                        concentrations_absolute_uncertainties_one_time_revised = (implicitSLScorrectionAmounts1sigma**2 + uncertainties_dict['concentrations_absolute_uncertainties_one_time']**2)**0.5
                        #Now to use the same ratio on the relative uncertainties.
                        concentrations_relative_uncertainties_one_time_revised = uncertainties_dict['concentrations_relative_uncertainties_one_time']*concentrations_absolute_uncertainties_one_time_revised/uncertainties_dict['concentrations_absolute_uncertainties_one_time']
                        #Now to populate both variables in the uncertainties dictionary.
                        uncertainties_dict['concentrations_relative_uncertainties_one_time']=concentrations_relative_uncertainties_one_time_revised
                        uncertainties_dict['concentrations_absolute_uncertainties_one_time']=concentrations_absolute_uncertainties_one_time_revised
                if timeIndex == 0: #NOTE: In the below lines we use a variable outside the dictionary because the dictionary gets reset with each iteration.
                    concentrations_relative_uncertainties_all_times = uncertainties_dict['concentrations_relative_uncertainties_one_time'] *1.0
                    concentrations_absolute_uncertainties_all_times = uncertainties_dict['concentrations_absolute_uncertainties_one_time']*1.0
                elif timeIndex > 0: #Everything else is appended via numpy.vstack
                    concentrations_absolute_uncertainties_all_times = numpy.vstack((concentrations_absolute_uncertainties_all_times,uncertainties_dict['concentrations_absolute_uncertainties_one_time']))
                    concentrations_relative_uncertainties_all_times = numpy.vstack((concentrations_relative_uncertainties_all_times,uncertainties_dict['concentrations_relative_uncertainties_one_time']))
            if G.iterativeAnalysis: #If using iterative analysis, append the subtracted signals' matching correction values that were used to SS_reciprocal_matching_correction_values_TimesList
                SS_reciprocal_matching_correction_values_TimesList.append(currentReferenceData.SSreciprocal_matching_correction_values)
                
            

                            
                

                        
        concentrationsScaledToCOarray[:,1:] = concentrationsScaledToCOarray[:,1:]/G.scaleRawDataFactor #correct for any scaling factor that was used during analysis.
        if (G.calculateUncertaintiesInConcentrations == True):#correct for any scaling factor that was used during analysis.
            concentrations_absolute_uncertainties_all_times = concentrations_absolute_uncertainties_all_times/G.scaleRawDataFactor
            #The relative absolute uncertainties should not be divided.
        
        #Now store things in resulsObjects.
        resultsObjects['concentrationsScaledToCOarray'] = concentrationsScaledToCOarray #Store in the global resultsObjects dictionary
        if (G.calculateUncertaintiesInConcentrations == True):
            resultsObjects['concentrations_absolute_uncertainties_all_times'] = concentrations_absolute_uncertainties_all_times*1.0
            resultsObjects['concentrations_relative_uncertainties_all_times'] = concentrations_relative_uncertainties_all_times*1.0
        if G.concentrationFinder == 'yes': #If using concentration finder
            concentrationsarray = copy.copy(concentrationsScaledToCOarray) #point concentrationsarray to a copy of concentrationsScaledToCOArray
            concentrationsarray[:,1:] = concentrationsarray[:,1:]*conversionFactorsAtEachTime #Multiply the data points by the appropriate conversion factor
            resultsObjects['concentrationsarray'] = concentrationsarray
            if (G.calculateUncertaintiesInConcentrations == True): #If there is aconversion factor that also gets applies to absolute_uncertainties
                resultsObjects['concentrations_absolute_uncertainties_all_times_scaled_to_unit'] = concentrations_absolute_uncertainties_all_times*conversionFactorsAtEachTime
        print('Data Analysis Finished.')
        #show net time for Data Analysis
        G.timeSinceLastCheckpoint = timeit.default_timer() - G.checkpoint
        G.checkPoint = timeit.default_timer()
        print('Data Analysis Time: ', (G.timeSinceLastCheckpoint))
        
        #this section exports and graphs the analyzed signals 
        if G.generatePercentages == 'yes':
            percentagesOutputArray = GeneratePercentages(concentrationsScaledToCOarray)
            resultsObjects['percentagesOutputArray'] = percentagesOutputArray
            ExportXYYYData(G.scaledConcentrationsPercentages, percentagesOutputArray, currentReferenceData.molecules, fileSuffix = G.iterationSuffix)
        ExportXYYYData(G.resolvedScaledConcentrationsOutputName, concentrationsScaledToCOarray, currentReferenceData.molecules, abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = str('scaled'))
        times = concentrationsScaledToCOarray[:,0]#the times are just the first column of the array
        data = concentrationsScaledToCOarray[:,1:]#the data is the whole array except the first column, which is the times 
        if (G.calculateUncertaintiesInConcentrations == True): #This is scaled to CO, not scaled to unit (which is concentration finder case).
            times_array_2D = numpy.transpose(numpy.atleast_2d(times)) #need to make it a 2D array with right direction to stack.
            concentrations_absolute_uncertainties_all_times_with_abscissa = numpy.hstack((times_array_2D,resultsObjects['concentrations_absolute_uncertainties_all_times']))
            concentrations_relative_uncertainties_all_times_with_abscissa = numpy.hstack((times_array_2D,resultsObjects['concentrations_relative_uncertainties_all_times']))
            if '.tsv' in G.resolvedScaledConcentrationsOutputName:                            
                ExportXYYYData(G.resolvedScaledConcentrationsOutputName[:-4]+"_absolute_uncertainties.tsv" , concentrations_absolute_uncertainties_all_times_with_abscissa, currentReferenceData.molecules, abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = 'scaled', units = None)
                ExportXYYYData(G.resolvedScaledConcentrationsOutputName[:-4]+"_relative_uncertainties.tsv", concentrations_relative_uncertainties_all_times_with_abscissa, currentReferenceData.molecules, abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = 'relative_uncertainty', units = None)
            if '.csv' in G.resolvedScaledConcentrationsOutputName:                            
                ExportXYYYData(G.resolvedScaledConcentrationsOutputName[:-4]+"_absolute_uncertainties.csv" , concentrations_absolute_uncertainties_all_times_with_abscissa, currentReferenceData.molecules, abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = 'scaled', units = None)
                ExportXYYYData(G.resolvedScaledConcentrationsOutputName[:-4]+"_relative_uncertainties.csv", concentrations_relative_uncertainties_all_times_with_abscissa, currentReferenceData.molecules, abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = 'relative_uncertainty', units = None)

        
        G.generateStatistics = True #FIXME: This variable needs to be added to user Input somehow. Also does not work with iterative.
        if G.generateStatistics == True and G.iterativeAnalysis == False:
            ScaledConcentrations_Statistics = []
            import numpy as np
            if '.csv' in G.resolvedScaledConcentrationsOutputName:
                ScaledConcentrations = numpy.genfromtxt(G.resolvedScaledConcentrationsOutputName, delimiter=',',skip_header=1)
            if '.tsv' in G.resolvedScaledConcentrationsOutputName:
                ScaledConcentrations = numpy.genfromtxt(G.resolvedScaledConcentrationsOutputName, delimiter='\t',skip_header=1)
            x_axis = ScaledConcentrations[:,0]
            ScaledConcentrations_data = ScaledConcentrations[:,1:]
            ScaledConcentrations_Means = numpy.atleast_2d(numpy.mean(ScaledConcentrations_data, axis=0)) #atleast_2D is just making it 2D
            ScaledConcentrations_StdDevs = numpy.atleast_2d(numpy.std(ScaledConcentrations_data, axis=0))
            ScaledConcentrations_StdErrs = ScaledConcentrations_StdDevs/(len(x_axis)**0.5)
            ScaledConcentrations_StdErrs_data = ScaledConcentrations_StdErrs #storing this for use below.
            # numpy.savetxt("ScaledConcentrations_TimeRange_Means.csv", ScaledConcentrations_Means, delimiter = ',', fmt ="%s")
            # numpy.savetxt("ScaledConcentrations_TimeRange_StdDevs.csv", ScaledConcentrations_StdDevs, delimiter = ',', fmt ="%s")
            # numpy.savetxt("ScaledConcentrations_TimeRange_StdErrs.csv", ScaledConcentrations_StdDevs/(len(x_axis)**0.5), delimiter = ',', fmt ="%s")
            ScaledConcentrations_Means = ["ScaledConcentrations_Means"] + list(ScaledConcentrations_Means.flatten())
            ScaledConcentrations_StdDevs = ["ScaledConcentrations_StdDevs"] + list(ScaledConcentrations_StdDevs.flatten())
            ScaledConcentrations_StdErrs = ["ScaledConcentrations_StdErrs"] + list(ScaledConcentrations_StdErrs.flatten())
            #We don't want things to crash if the uncertainties have a problem, so will use try and except.
            ScaledConcentrations_Statistics = ScaledConcentrations_Statistics + [ScaledConcentrations_Means,ScaledConcentrations_StdDevs,ScaledConcentrations_StdErrs]
            if (G.calculateUncertaintiesInConcentrations == True):
                if '.tsv' in G.resolvedScaledConcentrationsOutputName:
                    try:
                        ScaledConcentrations_absolute_uncertainties = numpy.genfromtxt( "ScaledConcentrations_absolute_uncertainties.tsv", delimiter='\t',skip_header=1, encoding='utf8')
                    except:
                        ScaledConcentrations_absolute_uncertainties = numpy.genfromtxt( "ScaledConcentrations_absolute_uncertainties.tsv", delimiter='\t',skip_header=1, encoding='utf16')
                if '.csv' in G.resolvedScaledConcentrationsOutputName:
                    ScaledConcentrations_absolute_uncertainties = numpy.genfromtxt( "ScaledConcentrations_absolute_uncertainties.csv", delimiter=',',skip_header=1)
                
                x_axis = ScaledConcentrations_absolute_uncertainties[:,0]
                #The "A" below is for absolute uncertainties, the "R" is for relative uncertainties.
                ScaledConcentrations_absolute_uncertainties_data = ScaledConcentrations_absolute_uncertainties[:,1:]
                ScaledConcentrations_absolute_uncertainties_Means = numpy.atleast_2d(numpy.mean(ScaledConcentrations_absolute_uncertainties_data, axis=0)) #atleast_2D is just making it 2D
                ScaledConcentrations_absolute_uncertainties_StdDevs = numpy.atleast_2d(numpy.std(ScaledConcentrations_absolute_uncertainties_data, axis=0))
                ScaledConcentrations_relative_uncertainties_Means =   numpy.atleast_2d  (ScaledConcentrations_absolute_uncertainties_Means.flatten()/np.array(ScaledConcentrations_Means[1:]).flatten()) #Need to skip first value since now that is a word.
                ScaledConcentrations_relative_uncertainties_StdDevs =  numpy.atleast_2d  (ScaledConcentrations_absolute_uncertainties_StdDevs.flatten()/np.array(ScaledConcentrations_Means[1:]).flatten()) #Need to skip first value since now that is a word.
                ScaledConcentrations_absolute_uncertainties_Total = numpy.atleast_2d ( (ScaledConcentrations_StdErrs_data.flatten()**2 + ScaledConcentrations_absolute_uncertainties_Means.flatten()**2)**0.5 )
                # numpy.savetxt("ScaledConcentrations_TimeRange_Absolute_Uncertainties_Means.csv", ScaledConcentrations_absolute_uncertainties_Means, delimiter = ',', fmt ="%s")
                # numpy.savetxt("ScaledConcentrations_TimeRange_Absolute_Uncertainties_StdDevs.csv", ScaledConcentrations_absolute_uncertainties_StdDevs, delimiter = ',', fmt ="%s")
                # numpy.savetxt("ScaledConcentrations_TimeRange_Relative_Uncertainties_Means.csv", ScaledConcentrations_relative_uncertainties_Means, delimiter = ',', fmt ="%s")
                ScaledConcentrations_absolute_uncertainties_Means = ["ScaledConcentrations_absolute_uncertainties_Means"] + list(ScaledConcentrations_absolute_uncertainties_Means.flatten())
                ScaledConcentrations_absolute_uncertainties_StdDevs = ["ScaledConcentrations_absolute_uncertainties_StdDevs"] + list(ScaledConcentrations_absolute_uncertainties_StdDevs.flatten())
                ScaledConcentrations_relative_uncertainties_Means = ["ScaledConcentrations_relative_uncertainties_Means"] + list(ScaledConcentrations_relative_uncertainties_Means.flatten())
                ScaledConcentrations_relative_uncertainties_StdDevs = ["ScaledConcentrations_relative_uncertainties_StdDevs"] + list(ScaledConcentrations_relative_uncertainties_StdDevs.flatten())
                ScaledConcentrations_relative_uncertainties_TotUnc = ["ScaledConcentrations_absolute_uncertainties_Total"] + list(ScaledConcentrations_absolute_uncertainties_Total.flatten())
                ScaledConcentrations_Statistics = ScaledConcentrations_Statistics + [ScaledConcentrations_absolute_uncertainties_Means, ScaledConcentrations_absolute_uncertainties_StdDevs,ScaledConcentrations_relative_uncertainties_Means, ScaledConcentrations_relative_uncertainties_StdDevs, ScaledConcentrations_relative_uncertainties_TotUnc]
            if '.tsv' in G.resolvedScaledConcentrationsOutputName:
                moleculesHeader = 'molecules \t' + '\t'.join(currentReferenceData.molecules)
                numpy.savetxt(G.resolvedScaledConcentrationsOutputName.replace(".tsv",'') + "_Statistics.tsv", np.array(ScaledConcentrations_Statistics), delimiter = '\t', fmt ="%s", header=moleculesHeader)
            if '.csv' in G.resolvedScaledConcentrationsOutputName:
                moleculesHeader = 'molecules ,' + ', '.join(currentReferenceData.molecules)
                numpy.savetxt(G.resolvedScaledConcentrationsOutputName.replace(".csv",'') + "_Statistics.csv", np.array(ScaledConcentrations_Statistics), delimiter = ',', fmt ="%s", header=moleculesHeader)

            
        #Graph the concentration/relative signal data
        if G.grapher == 'yes':
            Draw(times, data, currentReferenceData.molecules, G.concentrationFinder, G.unitsTSC, graphFileName='scaledConcentrationsAfterAnalysis', fileSuffix = G.iterationSuffix, label="Relative Concentrations Graph", stopAtGraphs=G.stopAtGraphs, figureNumber= G.lastFigureNumber+1)
            G.lastFigureNumber = G.lastFigureNumber+1
            if (G.calculateUncertaintiesInConcentrations == True):
                G.lastFigureNumber = DrawUncertaintiesConcentrations(times, data, resultsObjects['concentrations_absolute_uncertainties_all_times'], currentReferenceData.molecules, G.concentrationFinder, G.unitsTSC, graphFileName='scaledConcentrationsAfterAnalysis_uncertainties', fileSuffix = G.iterationSuffix, label="Concentrations With Uncertainties", stopAtGraphs=G.stopAtGraphs, figureNumber= G.lastFigureNumber+1)

            
        if G.concentrationFinder == 'yes':
            ExportXYYYData(G.concentrationsOutputName, concentrationsarray, currentReferenceData.molecules, abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = 'concentration', units = G.unitsTSC)
            times = concentrationsarray[:,0]
            data = concentrationsarray[:,1:]
            if (G.calculateUncertaintiesInConcentrations == True):
                concentrations_absolute_uncertainties_all_times_with_abscissa = numpy.hstack((times_array_2D,resultsObjects['concentrations_absolute_uncertainties_all_times_scaled_to_unit']))
                concentrations_relative_uncertainties_all_times_with_abscissa = numpy.hstack((times_array_2D,resultsObjects['concentrations_relative_uncertainties_all_times']))
                ExportXYYYData(G.concentrationsOutputName[:-4]+"_absolute_uncertainties.csv" , concentrations_absolute_uncertainties_all_times_with_abscissa, currentReferenceData.molecules, abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = 'concentration', units = G.unitsTSC)
                ExportXYYYData(G.concentrationsOutputName[:-4]+"_relative_uncertainties.csv", concentrations_relative_uncertainties_all_times_with_abscissa, currentReferenceData.molecules, abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = 'relative_uncertainty', units = None)

        
            #Graph the concentration/relative signal data
            if G.grapher == 'yes':
                Draw(times, data, currentReferenceData.molecules, G.concentrationFinder, G.unitsTSC, graphFileName='resolvedConcentrationsAfterAnalysis', fileSuffix = G.iterationSuffix, label="Concentrations", stopAtGraphs=G.stopAtGraphs, figureNumber= G.lastFigureNumber+1)
                G.lastFigureNumber = G.lastFigureNumber+1
                if (G.calculateUncertaintiesInConcentrations == True):
                    G.lastFigureNumber = DrawUncertaintiesConcentrations(times, data, resultsObjects['concentrations_absolute_uncertainties_all_times_scaled_to_unit'], currentReferenceData.molecules, G.concentrationFinder, G.unitsTSC, graphFileName='resolvedConcentrationsAfterAnalysis_uncertainties', fileSuffix = G.iterationSuffix, label="Concentrations With Uncertainties", stopAtGraphs=G.stopAtGraphs, figureNumber= G.lastFigureNumber+1)

            
    if G.dataSimulation =='yes':
        
        if G.dataAnalysis == 'skip':
            print("Warning: Simulation cannot be run without data")
            sys.exit()
        if G.dataAnalysis == 'load':
            #This function call loads the preprocessed data
            concentrationsScaledToCOarray = ImportAnalyzedData(G.resolvedScaledConcentrationsOutputName)
        #reset timer so that data simiulation can be timed
        G.checkpoint = timeit.default_timer()
        if G.iterativeAnalysis:
            #TODO when RawSignalSimulation is rewritten for multiple reference patterns, ReferenceDataSSreciprocal_matching_correction_valuesList[0] should be replaced with ReferenceDataSSreciprocal_matching_correction_valuesList
            SS_reciprocal_matching_correction_values_TimesArray = numpy.stack(SS_reciprocal_matching_correction_values_TimesList) #use numpy.stack to make an array of all the arrays
            #now the simulated data function uses the answer array and finds what the raw signals would be produced with their signals
            if G.scaleRawDataFactor == 1.0: 
                simulateddata = RawSignalsSimulation(concentrationsScaledToCOarray, SS_reciprocal_matching_correction_values_TimesArray)
            else:
                print("Warning: With iterative analysis, using scaleRawDataOption (also called scaleRawDataYorN) ends up scaling the simulated raw signals, also.")
                simulateddata = RawSignalsSimulation(concentrationsScaledToCOarray, SS_reciprocal_matching_correction_values_TimesArray)*G.scaleRawDataFactor
        if not G.iterativeAnalysis:
            correctionFactorsAtEachTime = numpy.stack(correctionFactorArraysList) #use numpy.stack to make an array of all the arrays
            #now the simulated data function uses the answer array and finds what the raw signals would be produced with their signals
            simulateddata = RawSignalsSimulation(concentrationsScaledToCOarray, correctionFactorsAtEachTime)
        #Exporting the simulated signals data
        if G.iterativeAnalysis:
            ExportXYYYData(G.simulatedSignalsOutputName, simulateddata, ExperimentDataCopy.mass_fragment_numbers, abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = 'simulated')
        if not G.iterativeAnalysis:
            ExportXYYYData(G.simulatedSignalsOutputName, simulateddata, ExperimentData.mass_fragment_numbers, abscissaHeader = ExperimentData.abscissaHeader, fileSuffix = G.iterationSuffix, dataType = 'simulated')
        #show net time for simulation
        G.timeSinceLastCheckPoint = timeit.default_timer() - G.checkpoint
        G.checkPoint = timeit.default_timer()
        print("Simulation Finished")
        print('Simulation Time: ', (G.timeSinceLastCheckPoint))
        
        resultsObjects['simulateddata'] = simulateddata
        
    CreateLogFile()    
    PopulateLogFile()
    print("LogFile complete")

    if G.iterativeAnalysis:
        IterativeAnalysisPostProcessing(ExperimentData, simulateddata, ExperimentDataCopy.mass_fragment_numbers, ExperimentDataFullCopy, times, data, currentReferenceData.molecules)

if __name__ == '__main__':
    main()
