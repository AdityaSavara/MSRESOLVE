import copy 
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
import UserInput 
G = UserInput
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
        
#this function works by drawing the input values from the input file which includes a baslineType, a time range, and the mass fragment
#that is being altered. the time range is taken, and using a for loop, all the values between those two times (including those two times)
#are used in finding the average or polyfit, the collected signals are also picked up here, in order to use them for finding the polyfit. 
#So then, depending of the baslineType, the mass fragments row will be found using a nested for loop with an if statement, which then
# subtracts what ever was found- the polyfit or the average from each data point in the collected data set. 
def LinearBaselineCorrectorSemiAutomatic(ExperimentData,baslineType,massesToBackgroundCorrect,earlyBaselineTimes,lateBaselineTimes):
    
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

        #these are the if statements that choose what happens based on user baslineType in the input data file
    if len(baslineType) == 1:
        baslineTypeholder = []
        for length in range(len(massesToBackgroundCorrect)):
            baslineTypeholder.append(baslineType[0])
        baslineType = baslineTypeholder
    for MassFragmentIndex in range(len(massesToBackgroundCorrect)):#array-indexed for loop
        if baslineType[MassFragmentIndex] == 'flat': #the different baslineTypes subtract different things
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
        if baslineType[MassFragmentIndex] == 'linear':#other option
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
THis function determines and returns the ABC correction coefficients based off 
of a literature (NIST) reference pattern and a measured (Experimental)reference pattern
'''
def ABCDetermination(ReferencePatternMeasured, ReferencePatternLiterature):
    '''
    Step 1: Read in all neccessary information from fragment patterns
    Step 2: populate all variables and standardize signals by 100 

    '''
    #first we import neccessary variables for Literature reference values: 
    MeasuredReferencePattern=ReferencePatternMeasured
    ReferenceFragmentationPattern=ReferencePatternLiterature

    #get reference numbers
    reference = genfromtxt( '%s' %ReferenceFragmentationPattern,delimiter=',',skip_header=1)

    #get molecule names: read in all as strings
    spamReader = csv.reader(open('%s' %ReferenceFragmentationPattern), delimiter=' ')
    list_holder=[]
    molecules_holder=[]
    #find row of strings with molecules
    for row in spamReader:#array-indexed for loop
        list_holder.append(row)
    molecules_holder = list_holder[1][0]
    #seperate molecules into iteratble entities
    for x in range(len(list_holder[1])-1):#array-indexed for loop
        molecules_holder = molecules_holder + ' ' +list_holder[1][x+1]
    molecules_holder = molecules_holder.split(',')
    MoleculesReference= molecules_holder[1:]

    #standardize to 100: Reference
    NewReference = StandardizeReferencePattern(reference[3:,:], len(MoleculesReference))
    MassFragmentsReferenceHolder=NewReference[0:,0]
    SignalsReferenceHolder=NewReference[:,1:]

    #IMPORT Experimental variables

    #get reference numbers
    referenceMeasured = genfromtxt( '%s' %MeasuredReferencePattern,delimiter=',',skip_header=1)
     
     #standardize to 100: Measured
    NewReferenceMeasured=StandardizeReferencePattern(referenceMeasured[3:,:], len(MoleculesReference))
    
    #set variables
    MassFragmentsMeasuredHolder=NewReferenceMeasured[0:,0]
    SignalsMeasuredHolder=NewReferenceMeasured[0:, 1:]
    matchingMassFragmentsHolder=[]

    #determine matching mass fragment numbers
    matchingMassFragmentsHolder=(MatchingNumbers(MassFragmentsMeasuredHolder, MassFragmentsReferenceHolder))
    matchingMassFragmentsHolder.sort(key = int)

    
    #set to numpy arrays
    MassFragmentsReference=numpy.array(MassFragmentsReferenceHolder)
    SignalsReference=numpy.array(SignalsReferenceHolder)
    MassFragmentsMeasured=numpy.array(MassFragmentsMeasuredHolder)
    SignalsMeasured=numpy.array(SignalsMeasuredHolder)

    

    '''
    Step 3: Determine Matching Mass fragments{future work will be done to do generate matching mass numbers}
     and generate Ratio array

    '''
    #later iterations will correct this
    MatchingMassFragments=numpy.array(matchingMassFragmentsHolder)


    #delete all indexes not at matching mass fragments: In Measured Data
    newMassFragmentsMeasured =[]
    newSignalsMeasured=[]
    MassIndex=0
    for  MatchingMassIndex in range(0, len(MatchingMassFragments)):

            if (MatchingMassFragments[MatchingMassIndex] != MassFragmentsMeasured[MassIndex]): 
                while(MatchingMassFragments[MatchingMassIndex] != MassFragmentsMeasured[MassIndex]):
                    MassIndex=MassIndex +1    
            newMassFragmentsMeasured.append(MassFragmentsMeasured[MassIndex])
            newSignalsMeasured.append(SignalsMeasured[MassIndex])

    

    MassFragmentsMeasured=numpy.array(newMassFragmentsMeasured)
    SignalsMeasured=numpy.array(newSignalsMeasured)

    #delete all indexes not at matching mass fragments: In Reference Data

    newMassFragmentsReference =[]
    newSignalsReference=[]
    MassIndex=0
    for  MatchingMassIndex in range(0, len(MatchingMassFragments)):

            if (MatchingMassFragments[MatchingMassIndex] != MassFragmentsReference[MassIndex]): 
                while(MatchingMassFragments[MatchingMassIndex] != MassFragmentsReference[MassIndex]):
                    
                    MassIndex=MassIndex +1    
            newMassFragmentsReference.append(MassFragmentsReference[MassIndex])
            newSignalsReference.append(SignalsReference[MassIndex])


    MassFragmentsReference=numpy.array(newMassFragmentsReference)
    SignalsReference=numpy.array(newSignalsReference)

    FactorListHolder=[]
    FactorList=[]
    MassFragmentsAdjustedHolder=[]

    

    #create and generate a factor list; also adjust mass fragment format to [x,x], [y,y]
    #instead of [x,y], [x,y]
    for  moleculeIndex in range(0, len(MoleculesReference)):
        for MassIndex in range (0, len(MatchingMassFragments)):
            
            if SignalsReference[MassIndex,moleculeIndex] !=0 and SignalsMeasured[MassIndex,moleculeIndex] != 0:
                Factor =(SignalsMeasured[MassIndex,moleculeIndex]/SignalsReference[MassIndex,moleculeIndex])
                
                FactorListHolder.append(Factor)
                MassFragmentsAdjustedHolder.append(MatchingMassFragments[MassIndex])




                
    FactorList=FactorListHolder
    MassFragmentsAdjusted=MassFragmentsAdjustedHolder

    
    '''
    Find a,b,c:
    '''


    (a,b,c)=numpy.polyfit(MassFragmentsAdjusted,FactorList,2)

    
    return a,b,c

 
#this function either creates or gets the three coefficients for the polynomial correction and calculates
#the correction factor for the relative intensities of each mass fragment, outputting a corrected set
#of relative intensities
def CorrectionValueCorrector(reference,referenceCorrectionCoefficients,referenceFileDesiredTuningAndForm,referenceFileExistingTuningAndForm,tuningCorrection):
    
    if tuningCorrection =='yes':
        (referenceCorrectionCoefficients['A'],referenceCorrectionCoefficients['B'],referenceCorrectionCoefficients['C'])=ABCDetermination(referenceFileExistingTuningAndForm,referenceFileDesiredTuningAndForm )
    
    referenceabscissa = reference[:,0] #gets arrays of just data and abscissa
    referencedata = reference[:,1:]
    for massfrag_counter in range(len(referenceabscissa)):#array-indexed for loop, only the data is altered, based on the abscissa (mass-dependent correction factors)
        factor = referenceCorrectionCoefficients['A']*(referenceabscissa[massfrag_counter]**2)  + referenceCorrectionCoefficients['B']*referenceabscissa[massfrag_counter]+referenceCorrectionCoefficients['C'] #obtains the factor from molecular weight of abscissa
        referencedata[massfrag_counter,:] = referencedata[massfrag_counter,:]*factor
    reference[:,0] = referenceabscissa
    return reference
    
        
#this function eliminates any fragments that are below a certain threshold, for solving 
#data that is giving negatives or over emphasizing small mass fragments, this will eliminate 
#those below a certain user-input value
def ReferenceThreshold(reference,referenceMassFragmentFilterThreshold):
    referencedata = reference[:,1:] #all the data except the line of abscissa- mass fragment numbers
    for rowcounter in range(len(referencedata[:,0])):#goes through all rows in references
        for columncounter in range(len(referencedata[0,:])):#goes through all columns in all rows in reference
            if referencedata[rowcounter,columncounter] < referenceMassFragmentFilterThreshold:#user input changes
                referencedata[rowcounter,columncounter] = 0 #made to be equal to zero
    reference[:,1:] = referencedata #this puts changed reference back with mass fragment numbers
    return reference
    
    
#The function just uses the input sheet to get values to change the reference sheet
#this is done by getting a list of all the values needed from the collected sheet at
#the right time and then gets each number following the first, and finds its ratio
#with the first, and multiplies that number by the number in the reference sheet in 
#order to change the second mass fragments number in the table
def ExtractReferencePatternFromData (ExperimentData, ReferenceData, rpcChosenMolecules,rpcChosenMoleculesMF,rpcTimeRanges):
    copyOfReferenceData = copy.deepcopy(ReferenceData)    
    for chosenmoleculescounter in range(len(rpcChosenMolecules)):#array-indexed for loop
        extractedIntensities = []
        allExtractedIntensities = []
        massfragindexer = []
        for moleculecounter in range(len(copyOfReferenceData.molecules)):#array-indexed for loop
            if copyOfReferenceData.molecules[moleculecounter] == rpcChosenMolecules[chosenmoleculescounter]:#finds index of molecule
                if len(rpcChosenMoleculesMF[chosenmoleculescounter]) == 1:#if only one number is given then the function changes all the other values according to this one
                    for x in range(len(ExperimentData.mass_fragment_numbers)):#array-indexed for loop
                        if ExperimentData.mass_fragment_numbers[x] != 0:#finds the other mass fragments and appends them
                            rpcChosenMoleculesMF[chosenmoleculescounter].append(ExperimentData.mass_fragment_numbers[x])
                            rpcTimeRanges[chosenmoleculescounter].append(rpcTimeRanges[chosenmoleculescounter][0])
                            if ExperimentData.mass_fragment_numbers[x] == rpcChosenMoleculesMF[chosenmoleculescounter][0]:#if the mass fragment is equal to the one being checked with, it is not need so it is deleted
                                rpcChosenMoleculesMF[chosenmoleculescounter].pop()
                                rpcTimeRanges[chosenmoleculescounter].pop()
                for eachChosenMoleculesMF in range(len(rpcChosenMoleculesMF[chosenmoleculescounter])):#array-indexed for loop
                    for refMFCounter in range(len(copyOfReferenceData.mass_fragment_numbers_monitored)): #checks whole input list (or list that was made by previous loop)
                        if rpcChosenMoleculesMF[chosenmoleculescounter][eachChosenMoleculesMF] == copyOfReferenceData.mass_fragment_numbers_monitored[refMFCounter]:#gets index of mass fragment number
                            massfragindexer.append(refMFCounter)
                for eachChosenMoleculesMF in range(len(rpcChosenMoleculesMF[chosenmoleculescounter])):
                    for refMFCounter in range(len(ExperimentData.mass_fragment_numbers)):
                        if rpcChosenMoleculesMF[chosenmoleculescounter][eachChosenMoleculesMF] == ExperimentData.mass_fragment_numbers[refMFCounter]:
                            for timecounter in range(len(ExperimentData.times)):#array-indexed for loop
                                if (rpcTimeRanges[chosenmoleculescounter][0] <= ExperimentData.times[timecounter]) and (rpcTimeRanges[chosenmoleculescounter][1] >= ExperimentData.times[timecounter]):#gets index of time
                                    extractedIntensities.append(ExperimentData.workingData[timecounter,refMFCounter])
                            #Place current extractedIntensities in a larger list using copy so it can be cleared without affecting allExtractedIntensities
                            allExtractedIntensities.append(copy.copy(extractedIntensities))
                            #clear extracted intensities for the next mass fragment
                            extractedIntensities.clear()
                #Convert list to an array
                allExtractedIntensitiesArray = numpy.array(allExtractedIntensities)
                #Initialize empty list to store average values
                allExtractedIntensitiesAverage = []
                #For loop to find the average intensity values
                #Then store the average value in the allExtractedIntensityAverage list
                for eachChosenMoleculesMF in range(len(allExtractedIntensitiesArray)):
                    intensitiesAverage = numpy.average(allExtractedIntensitiesArray[eachChosenMoleculesMF])
                    allExtractedIntensitiesAverage.append(intensitiesAverage)
                #For loop to overwrite the reference file with the value of the reference signal of the chosen mass fragment with the product of the signal and the ratio of the averages
                for eachChosenMoleculesMF in range(len(rpcChosenMoleculesMF[chosenmoleculescounter])):
                    copyOfReferenceData.provided_reference_intensities[massfragindexer[eachChosenMoleculesMF],moleculecounter+1] = (allExtractedIntensitiesAverage[eachChosenMoleculesMF]/allExtractedIntensitiesAverage[0])*copyOfReferenceData.provided_reference_intensities[massfragindexer[eachChosenMoleculesMF],moleculecounter+1]
    return copyOfReferenceData.provided_reference_intensities

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

    return (massFragmentNumbers, workingData)


'''
MassFragChooser() compares chosenMassFragments and ExperimentData.mass_fragment_numbers. If there is a 
mass in mass_fragment_nubmers that isn't in chosenMassFragments it removes that mass from
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

    return (massFragmentNumbers,workingData)


    
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
    
#This function operates in a parallel way to trimDataMasses, but it operates on the reference data and all of it's constituent variables  
def TrimDataMolecules(ReferenceData, chosenMolecules):
    
    print("MoleculeChooser")
    #the copy is required because the keep only selected columns function is called twice with the same rows to clear
    copy_moleculeselecNum = copy.deepcopy(ReferenceData.molecules)
    
    #shorten the reference fragmentation pattern to the required length
    (Temp_Reference_Data, ReferenceData.molecules) = DataFunctions.KeepOnlySelectedYYYYColumns(ReferenceData.provided_reference_intensities[:,1:],
                                                                                                                    ReferenceData.molecules, chosenMolecules)  
    #Shorten the electronnumbers to the correct values, using the copy of molecules 
    ArrayOneD = True
    (ReferenceData.electronnumbers, copy_molecules) = DataFunctions.KeepOnlySelectedYYYYColumns(ReferenceData.electronnumbers, copy_moleculeselecNum, chosenMolecules, ArrayOneD)
    #add a second dimension to the reference data
    newReferenceMF = numpy.reshape(ReferenceData.mass_fragment_numbers_monitored,(-1,1))
    
    #Add the abscissa back into the reference values
    ReferenceData.provided_reference_intensities = numpy.hstack((newReferenceMF,Temp_Reference_Data))
    
    #remove any zero rows that may have been created
    ReferenceData.ClearZeroRows()
    #update the mass fragment list from the posibly shortened reference spectrums
    ReferenceData.mass_fragment_numbers_monitored = ReferenceData.provided_reference_intensities[:,0]
    
    ReferenceData.ExportCollector("MoleculeChooser")
    
    return ReferenceData.provided_reference_intensities, ReferenceData.electronnumbers, ReferenceData.molecules, ReferenceData.mass_fragment_numbers_monitored
    
'''
trimDataMasses() is just a wrapper function for two calls to DataFunctions.KeepOnlySelectedYYYYColumns(). 
Both of the calls trim ExperimentData.workingData and ExperimentData.mass_fragment_numbers. 
The first call trims the data according to the mass fragment selections in G.chosenMassFragments.
The second call trims the data to remove any mass fragments for which there is no ReferenceData. 

Parameters:
ExperimentData - of type MSData, the one instantiated in main() named ExperimentData is a good example of one
    that will work here
ReferenceData - of type MSReference, ReferenceData from main() is a good example
'''
def trimDataMasses(ExperimentData, ReferenceData):

    # If we are only interested in a subset of the MS data
    # and that subset is a subset of the loaded data
    # remove the irrelevant mass data series from ExperimentData.mass_fragment_numbers
    # and the corresponding colums from ExperimentData.workingData
    if G.specificMassFragments == 'yes':
        print("MassFragChooser")
        if len(G.chosenMassFragments) < len(ReferenceData.molecules):
            print("Selected Mass Fragments are too few to solve for the number of molecules provided")
            print("Mass fragment selection has been canceled")
        else:
            (ExperimentData.workingData, ExperimentData.mass_fragment_numbers) = DataFunctions.KeepOnlySelectedYYYYColumns(ExperimentData.workingData,
                                                                                                            ExperimentData.mass_fragment_numbers,
                                                                                                            G.chosenMassFragments)
        ExperimentData.ExportCollector("MassFragChooser")

    # Remove elements of ExperimentData.mass_fragment_numbers for which there is no matching mass in the reference data.
    # Also remove the corresponding mass data column from Experiment.workingData.
    (ExperimentData.workingData, ExperimentData.mass_fragment_numbers) = DataFunctions.KeepOnlySelectedYYYYColumns(ExperimentData.workingData,
                                                                                                        ExperimentData.mass_fragment_numbers,
                                                                                                        ReferenceData.provided_reference_intensities[:,0])

    return (ExperimentData.workingData, ExperimentData.mass_fragment_numbers)

    
#The correction values are calculated based on each molecule's molecular weight, number of electrons and mass fragmentation pattern
#This obtains the ionization energy, transmission and electron multiplier gains from the mw, e- and the application of simple equations
#Then it gets the sum of each molecules frag pattern and uses all these values to get each moleule's correction value via the
#system created by Madix and Ko and puts this in an answers array
def CorrectionValuesObtain(ReferenceData):
    reference_width = len(ReferenceData.standardized_reference_intensities[0,:]) 
    reference_height = len(ReferenceData.standardized_reference_intensities[:,0]) 
    correction_values = numpy.zeros([1,reference_height])
    #the first for loop here gets all of the values for e- and mw and uses them to get the
    #respective values that can find the correction factor for each mass fragment of each molecule
    for column_counter in range(1,reference_width): #array-indexed for loop, skips column one b/c that is the mass fragment numbers, not relative intensities
        ionization_efficiency = (0.6*ReferenceData.electronnumbers[column_counter-1]/14)+0.4
        answer_array_row = []
        quotients = numpy.zeros(len(ReferenceData.standardized_reference_intensities[:,0]))
        #This first loop goes through the row and gets all of the mass fragment's Tms and Gms
        #and then uses those along with the relative intensity of the mass fragment itself to
        #calculate the sum of each column's Fm/(Gm*Tm)
        for row_counter in range(0,reference_height):  #array-indexed for loop
            fragmentmass = ReferenceData.standardized_reference_intensities[row_counter,column_counter] 
            if fragmentmass != 0: #only gets the Gm and Tm if relative intensity is not equal to zero
                electron_multiplier_gain = (28/ReferenceData.standardized_reference_intensities[row_counter,0])**0.5 
                if (ReferenceData.standardized_reference_intensities[row_counter,0] < 30): #transmission gain depends on the mass fragment mass
                    transmission_gain = 1
                else:
                    transmission_gain = 10**((30-ReferenceData.standardized_reference_intensities[row_counter,0])/155)
                quotient = fragmentmass/(transmission_gain*electron_multiplier_gain)
                quotients[row_counter] = quotient
                a = sum(quotients)
        # This second part of the loop must be separate only because the sums of each column's 
        #Fm/(Gm*Tm) has to be found first, but then the loop gets each mass fragment's correction
        #value by dividing the sum by that relative intensity and ionization efficiency for the 
        #molecule, then creating an correction_values of correction values
        for row_counter in range(0,reference_height):#array-indexed for loop, there is a second loop for each column, because the first only gets the sum of (Fm)/(Gm*Tm)
            fragmentmass = ReferenceData.standardized_reference_intensities[row_counter,column_counter] 
            if fragmentmass != 0: #once again, if the relative intensity is zero, then the correction value will be zero as well, this is done by the else statement
                fragmentmass = ReferenceData.standardized_reference_intensities[row_counter,column_counter]
                correction = a/(ionization_efficiency*fragmentmass)
                answer_array_row.append(correction)
            else: 
                correction = 0
                answer_array_row.append(correction)
            if row_counter == (reference_height-1):#if the loop is on the last index
                if column_counter == 1: #the first column can be the beginning of the array
                    answer_array_row = numpy.array(answer_array_row)
                    correction_values = correction_values + answer_array_row
                else: #afterwards, all the rows are stacked
                    correction_values = numpy.vstack([correction_values,answer_array_row]) 
    return correction_values


#Populate_reciprocal_matching_correction_values-> for loops that build the reciprocal_matching_correction_values and raw signals into multiple arrays, inside of lists
#that are ready to be solved in the next section. It does this by using a itertools combination function in order
#to create all the combinations of row and thereby creating all the possible answer arrays for the data given
#need these indexes for later
def Populate_reciprocal_matching_correction_values(mass_fragment_numbers, ReferenceData):
    ReferenceData.referenceabscissa = ReferenceData.standardized_reference_patterns[:,0]
    referenceDataArray = ReferenceData.standardized_reference_patterns[:,1:]
    correction_values = numpy.array(list(zip(*ReferenceData.correction_values)))
    if G.calculateUncertaintiesInConcentrations == True:
        if type(G.referencePatterns_uncertainties) != type(None): #Just mimicing the above lines. 
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
        if type(G.referencePatterns_uncertainties) != type(None): #Just mimicing the above lines.
            ReferenceData.reciprocal_matching_correction_values_relative_uncertainties, ReferenceData.matching_abscissa = ArrayRowReducer(mass_fragment_numbers,ReferenceData.referenceabscissa,correction_values_relative_uncertainties)
    ReferenceData.reciprocal_matching_correction_values = ArrayElementsInverser(ReferenceData.matching_correction_values)
    return ReferenceData
    
    
#This function will take the reference data and eliminate any molecules from the data that do not contain any mass fragment
#data this is done so that there will not be errors in the code later (as other molecules may also have 0 signal 
#relative to CO) It does this by looking at the matching mass fragments and deleting any columns which contain only
#zeros, and then also deletes that molecule form the molecules array and the correction values array.
def  UnnecessaryMoleculesDeleter(ReferenceData):
    width = len(ReferenceData.monitored_reference_intensities[0,:])
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
                    place_holder = place_holder + 1
    return ReferenceData

    
#this little function lets you choose your own times range, the inputs are the start of the range,
#the end of the range, a 'yes' or 'no' (timerangelimit), and the times and collected arrays
#the collected data will be shortened to the length of the new chosen times abscissa
def  TimesChooser (ExperimentData,timeRangeStart,timeRangeFinish):
    place_holder = 0 #enables indexing when parts of the array are being deleted
    for timescounter in range(len(ExperimentData.times)): #array indexed for loop
        if ExperimentData.times[timescounter-place_holder] < timeRangeStart: #all rows that are before the time range are deleted from the collected data and times abscissa
            ExperimentData.times = numpy.delete(ExperimentData.times,timescounter-place_holder) #place holder subtracts from the for loop so that the correct index is maintained
            ExperimentData.workingData = numpy.delete(ExperimentData.workingData,timescounter-place_holder,axis = 0)
            place_holder = place_holder + 1 #the place holder increased by one with every deleted row to maintain array indexing
        if ExperimentData.times[timescounter-place_holder] > timeRangeFinish: #once the time is greater than the time range finish, all values after are deleted
            ExperimentData.times = numpy.delete(ExperimentData.times,timescounter-place_holder)
            ExperimentData.workingData = numpy.delete(ExperimentData.workingData,timescounter-place_holder,axis = 0)
            place_holder = place_holder + 1
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
        minNumber = float(numpy.min(a1DArray[numpy.nonzero(a1DArray)]))
        # Confirm that the array needs scaling
        if minNumber > Base:
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
    if scaleRawDataOption == 'auto':
        #automatically scales the smallest point to 1
        data = ScaleRawDataAuto(data)
    elif scaleRawDataOption == 'manual':
        #because a multiplier is given, it doesn't matter if ScaleUp or ScaleDown is used
        data, multiplier = ScaleUp(data, multiplier = scaleRawDataFactor)
    return data

''' ScaleRawDataAuto operates in a similar way to ScaleUp or Down, except that 
    it will always scale the smallest value to 1, regardless of the scaling direction'''
def ScaleRawDataAuto(data): 
    a1Dholder=[]
    a1Dholder=data[0]
    minNum=numpy.min(a1Dholder[numpy.nonzero(a1Dholder)])

    for index in range(0, len(data)):
        a1Dholder =data[index]
        
        minNumHolder=numpy.min(a1Dholder[numpy.nonzero(a1Dholder)])
        if (minNumHolder < minNum):
            minNum = minNumHolder
            
    ScaleFactor= 1/minNum

    for dataPointIndex in range(0, len(data)):
        for referenceIndex in range(0, len(data[0])):
            data[dataPointIndex, referenceIndex]=data[dataPointIndex, referenceIndex]*ScaleFactor
    return data

'''
Standardize is a simple algorithim that reads in one numpy array
(oneArray) and scales the whole array so that the
greatest number is 100 within the array
'''

def StandardizeTo100(a1DArray,n):
    maxNumber=float(numpy.max(a1DArray))
    multiplier=100/maxNumber
    for index2 in range(0,len(a1DArray)):
        a1DArray[index2]=(a1DArray[index2]*multiplier)
    return a1DArray



'''
StandardizeReferencePattern uses StandardizeTo100 to modify reference values so that each column is scaled to
100. NOTE THAT THIS FUNCTION ASSUMES THAT THE FIRST COLUMN in reference contains the 
mass framgnet numbers. 
Parameters: 
standardizedReference -  a name chosen for the numpy array that contains reference values
num_of_molecues-  an integer describing the number of  molecues that contributed to the reference file
'''
def StandardizeReferencePattern(referenceUnstandardized,num_of_molecules):
    # preallocate new array for standardized values
    standardizedReference = copy.deepcopy(referenceUnstandardized)

    # standardize
    for moleculeIndex in range(1,num_of_molecules+1):
        standardizedReference[0:,moleculeIndex]=StandardizeTo100(referenceUnstandardized[0:,moleculeIndex],1)

    return standardizedReference

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
    masses = dfmass.as_matrix()
    #sort through the matrix and remove labels
    #masses = numpy.delete(masses, -1)
    for i in range(0,len(masses)):
        masses[i] = masses[i].replace('m','')
    #convert the matrix to floats if they aren't already 
    mass_fragment_numbers = masses.astype(numpy.float)
    
    '''generate time list'''
    #select column of times
    dftimes = dataFrame.iloc[1:][0]
    #convert to matrix
    times = dftimes.as_matrix()
    #save with type float
    fulltimes = times.astype(numpy.float)
    
    '''collect preprocessed data'''
    #select matrix of signals
    dfpreprocessed = dataFrame.iloc[1:,1:]
    #convert to matrix
    preprocessed = dfpreprocessed.as_matrix()
    #save  with type float
    preprocessedData = preprocessed.astype(numpy.float)

    
    '''create data set to work on'''
    workingData = preprocessedData
    
    return [workingData, mass_fragment_numbers, fulltimes]

def ImportAnalyzedData(concentrationsOutputName):

    dataFrame = pandas.read_csv('%s' %concentrationsOutputName, header = None)

    '''collect preprocessed data'''
    #select matrix of signals
    dfanalyzed = dataFrame.iloc[1:,0:]
    #convert to matrix
    analyzed = dfanalyzed.as_matrix()
    analyzed = numpy.delete(analyzed, -1, 1)
    #save  with type float
    analyzedData = analyzed.astype(numpy.float)
    
    return analyzedData


'''
Performs some manipulations related to the reference pattern
'''
def ReferenceInputPreProcessing(ReferenceData):

    # standardize the reference data columns such that the maximum value is 100 and everything else is
    # linearly scaled according that the maximum value scaling
    ReferenceData.standardized_reference_intensities=StandardizeReferencePattern(ReferenceData.provided_reference_intensities,len(ReferenceData.molecules))
    ReferenceData.ExportCollector('StandardizeReferencePattern')
    
    print('beginning CorrectionValueCorrector')
    ReferenceData.standardized_reference_intensities = CorrectionValueCorrector(ReferenceData.standardized_reference_intensities, G.referenceCorrectionCoefficients,
                                                       G.referenceFileDesiredTuningAndForm, G.referenceFileExistingTuningAndForm,
                                                       G.tuningCorrection)
    ReferenceData.ExportCollector('CorrectionValueCorrector')
    
    if G.applyReferenceMassFragmentsThresholds == 'yes':
        ReferenceData.standardized_reference_intensities = ReferenceThreshold(ReferenceData.standardized_reference_intensities,G.referenceMassFragmentFilterThreshold)
        ReferenceData.ExportCollector('ReferenceThreshold')
    
    ReferenceData.correction_values = CorrectionValuesObtain(ReferenceData)
    print('CorrectionValuesObtain')

    return ReferenceData

'''
GenerateReferenceDataAndFormsList takes in the list of referenceFileNames and the
list of forms.  A list is generated containing MSReference objects created based
on the referenceFileName and the corresponding form
It allows MSRESOLVE to be backwards compatible with previous user input files
'''
def GenerateReferenceDataList(referenceFileNames,referencePatternForm):
    #referencePatternForm can take values of 'xyyy' or 'xyxy' and must be a string
    ##If referenceFileName is a string or if form is a string then make them lists
    if isinstance(referenceFileNames,str):
        referenceFileNames = [referenceFileNames]
    if isinstance(referencePatternForm,str):
        referencePatternForm = [referencePatternForm]
    #If referenceFileNames and forms are lists of 1 then create a list of the single MSReference object
    #This allows MSRESOLVE to be backwards compatible with previous user input files while still incorporating the reference pattern time chooser feature
    if len(referencePatternForm) == 1 and len(referenceFileNames) == 1:
        [provided_reference_intensities, electronnumbers, molecules, molecularWeights, sourceInfo, mass_fragment_numbers_monitored, referenceFileName, form]=readReferenceFile(referenceFileNames[0],referencePatternForm[0])
        ReferenceDataList = [MSReference(provided_reference_intensities, electronnumbers, molecules, molecularWeights, sourceInfo, mass_fragment_numbers_monitored, referenceFileName=referenceFileName, form=form)]
        return ReferenceDataList
    #Otherwise we have multiple reference files and forms
    #If just one form is used, make a list of forms that is the same length as referenceFileNames
    if len(referencePatternForm) == 1:
        #Generate a copy of referenceFileNames to be overwritten with forms
        listOfForms = copy.copy(referenceFileNames)
        #replace each value with the given form
        for i in range(len(referenceFileNames)):
            listOfForms[i] = referencePatternForm[0]
    #If list of forms is the same length of referenceFileNames then each form should correspond to the referenceFile of the same index
    elif len(referencePatternForm) == len(referenceFileNames):
        #So just set listOfForms equal to forms
        listOfForms = referencePatternForm
    #Initialize ReferenceDataAndFormsList so it can be appended to
    ReferenceDataAndFormsList = []
    #For loop to generate each MSReferenceObject and append it to a list
    for i in range(len(referenceFileNames)):
        [provided_reference_intensities, electronnumbers, molecules, molecularWeights, sourceInfo, mass_fragment_numbers_monitored, referenceFileName, form]=readReferenceFile(referenceFileNames[i],referencePatternForm[i])
        ReferenceDataAndFormsList.append(MSReference(provided_reference_intensities, electronnumbers, molecules, molecularWeights, sourceInfo, mass_fragment_numbers_monitored, referenceFileName=referenceFileName, form=form))
    return ReferenceDataAndFormsList


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
    
    #Scaling Raw Data
    #This if statement only applies ScaleRawData if scaling is set to automatic
    # or if there is a manual factor that is not 1.
    #We skip this function in manual factor was 1 because that would not change the data.
    if (G.scaleRawDataOption == 'auto' or (G.scaleRawDataOption == 'manual' and G.scaleRawDataFactor != 1)):
        ExperimentData.workingData  = ScaleRawData(ExperimentData.workingData, G.scaleRawDataOption, G.scaleRawDataFactor)
        ExperimentData.ExportCollector("ScaleRawData")

    #displays graph of raw data
    if G.grapher == 'yes':
        print("Raw Signal Graph")
        Draw(ExperimentData.times, ExperimentData.workingData, ExperimentData.mass_fragment_numbers, 'no', 'Amp', graphFileName = 'rawData', fileSuffix = G.iterationSuffix)
        
    if len(G.backgroundMassFragment) != 0:
        SlopeEliminator (ExperimentData,G.backgroundMassFragment,G.backgroundSlopes,G.backgroundIntercepts)
        print('Linear Baseline Correction, Manual, Complete')
        ExperimentData.ExportCollector("SlopeEliminator")

    if G.linearBaselineCorrectionSemiAutomatic   == 'yes': #the data edit sheet is used here, to determine to run this function or not
        LinearBaselineCorrectorSemiAutomatic(ExperimentData, G.baslineType, G.massesToBackgroundCorrect, G.earlyBaselineTimes, G.lateBaselineTimes)
        print('Linear Baseline Correction, Semiautomatic, Complete')
        ExperimentData.ExportCollector("LinearBaselineCorrectorSemiAutomatic")
        
    if G.lowerBoundThresholdChooser == 'yes':#calls above function if option is desired in data edit file
        LowerBoundThresholdFilter (ExperimentData, G.massesToLowerBoundThresholdFilter, G.lowerBoundThresholdPercentage, G.lowerBoundThresholdAbsolute)
        print('LowerBound Threshold Filter Complete')
        ExperimentData.ExportCollector("LowerBoundThresholdFilter")
        
    #displays mid-preprocessing graph    
    if G.grapher == 'yes':
        print("Pre-marginalChangeRestrictor Graph")
        Draw(ExperimentData.times, ExperimentData.workingData, ExperimentData.mass_fragment_numbers, 'no', 'Amp', graphFileName ='midProcessingGraph', fileSuffix = G.iterationSuffix)

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
#exports the user input file so that it can be used in the next iteration
def ExportUserInputFile(fileName):
    
    #Creating an updated UI file
    globalsFE_saveFile = fileName 
    globalsFE_loadFile = fileName 
    #create export object
    globalsFE_object = ei.module_export_import(globalsFE_saveFile,globalsFE_loadFile,G)
    
    #save variables to the text file 
    globalsFE_object.save_params()
    
#this function is used to append any list to a file in an executable fashion
def AppendList(listVariableName, List, FileName, entriesPerLine): 
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
           #if the previous entry was the last one in that line, add a newline character
           if (listIndex + 1) % entriesPerLine == 0: 
               f.write('\n')
        #write in the closing bracket for the list
        f.write(']')
    return None
    
def StringSearch(string, keyword = '', antikeyword = ''):
    if keyword in string and not antikeyword in string:
        return True
    else:
        return False

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
    if not suffixlist == []:
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
    #if the user has entered an iteration name
    if iterativeAnalysis == False or iterativeAnalysis == True:
         #create the default directory
        iterationDirectoryName = '.\\_iter_%s' %str(iterationNumber)
    else:
        #set that name to be the directory along with the correct number 
        iterationDirectoryName = '.\\' + iterativeAnalysis + '_iter_%s' %str(iterationNumber)
    return iterationDirectoryName

def IterationDirectoryPreparation(iterativeAnalysis, iterationNumber, iterate = False):
    #implied arguments for this function are G.referenceFileName and G.dataToAnalyzeFileName
    if iterate:
        iterationNumber += 1
    G.iterationNumber = iterationNumber
    
    iterationDirectoryName = SpecificIterationName(iterativeAnalysis, iterationNumber)
       
    #confirm that the directory exists
    EnsureDirectory(iterationDirectoryName)
    
    #Change the working directory to the new directory name. 
    'THIS IS A HIGHLY SIGNIFICANT LINE, because it redirects all of the output for the rest of the program run'
    os.chdir(iterationDirectoryName) #this will be changed back at the end of the program
    
    if not iterate:
        #record the old file names 
        G.oldReferenceFileName = G.referenceFileName
        G.oldDataToAnalyzeFileName = G.dataToAnalyzeFileName
        
        #construct the file names for the current run of the program
        referenceFileNameTemp = G.referenceFileName[:-18] +  str(G.iterationSuffix) + G.referenceFileName[-4:]
        dataToAnalyzeFileNameTemp = G.dataToAnalyzeFileName[:-21] +  str(G.iterationSuffix) + G.dataToAnalyzeFileName[-4:]
        
        #copy the experimental and reference files into new names for this iterative run
        shutil.copy(G.referenceFileName, referenceFileNameTemp)
        shutil.copy(G.dataToAnalyzeFileName, dataToAnalyzeFileNameTemp)
        
        #change the globals to reflect the renaming of the ref and exp files
        G.referenceFileName =  referenceFileNameTemp
        G.dataToAnalyzeFileName =  dataToAnalyzeFileNameTemp
        
        #construct file names for the next run of the program 
        G.nextRefFileName = G.referenceFileName[:-11] +  str('_unused') + G.referenceFileName[-11:]
        G.nextExpFileName = G.dataToAnalyzeFileName[:-11] +  str('_remaining') + G.dataToAnalyzeFileName[-11:]
    
    return None
    #implied returns: G.oldReferenceFileName, G.oldDataToAnalyzeFileName, G.referenceFileName,G.dataToAnalyzeFileName, G.nextRefFileName, G. nextExpFileName, G.iterationNumber 

def IterationFirstDirectoryPreparation(iterativeAnalysis,iterationNumber):
    #implied arguments for this function are G.referenceFileName and G.dataToAnalyzeFileName
    #this global value is set so that each export statement can label the output files correctly
    G.iterationNumber = iterationNumber
    
    iterationDirectoryName = SpecificIterationName(iterativeAnalysis, iterationNumber)
       
    #confirm that the directory exists
    EnsureDirectory(iterationDirectoryName)
    
    #Change the working directory to the new directory name. 
    'THIS IS A HIGHLY SIGNIFICANT LINE, because it redirects all of the output for the rest of the program run'
    os.chdir(iterationDirectoryName) #this will be changed back at the end of the program
    
    #copy the first UserInputFile into the first iteration directory
    ExportUserInputFile("UserInput_iter_1.py")
    #append the variable list to the user input file
    AppendList("__var_list__", G.__var_list__, "UserInput_iter_1.py", 5)
    
    #record the old file names 
    G.oldReferenceFileName = G.referenceFileName
    G.oldDataToAnalyzeFileName = G.dataToAnalyzeFileName
    #construct the file names for the first run of the program
    G.referenceFileName = G.referenceFileName[:-4] +  str(G.iterationSuffix) + G.referenceFileName[-4:]
    G.dataToAnalyzeFileName = G.dataToAnalyzeFileName[:-4] +  str(G.iterationSuffix) + G.dataToAnalyzeFileName[-4:]
    
    #construct file names for the second run of the program 
    G.nextRefFileName = G.referenceFileName[:-11] + '_unused_iter_1' + G.referenceFileName[-4:]
    
    G.nextExpFileName = G.dataToAnalyzeFileName[:-11] + '_remaining_iter_1' + G.dataToAnalyzeFileName[-4:]
    
    return None 
    #implied returns: G.oldReferenceFileName, G.oldDataToAnalyzeFileName, G.referenceFileName,G.dataToAnalyzeFileName, G.nextRefFileName, G. nextExpFileName, G.iterationNumber 

#The IterativeAnalysisDirectory and Variable Population function is used to shrink the size of the program analysis and redirect the output. 
def IADirandVarPopulation(iterativeAnalysis, chosenMassFragments, chosenMolecules, ExperimentData, ReferenceData, ReferenceDataFullCopy):
    #implied arguments: G.dataSimulation, G.referenceFileName, G.dataToAnalyzeFileName, G.nextRefFileName, G.oldReferenceFileName, G.chosenMolecules, G.iterationNumber
    #override data simulation to yes if it was not selected
    if G.dataSimulation != 'yes':
        print("Iterative analysis cannot find the remaining signals in the experiment without signal simulation being run.")
        print("User selection to skip signal simulation has been overridden. ")
        G.dataSimulation = 'yes'
    #Warn the user if they are trying to run an iteration that has no molecules to solve. (This would cause a complex error further on in the program if allowed to run.)    
    if len(ReferenceData.molecules) == 0:
        print("Warning Message: There are inadequate molecules to perform another iteration. Please confirm that there are still remaining molecules to solve.")
        sys.exit()
         
    #Selecting unused Reference Data
    unusedMolecules = []
    for molecule in ReferenceDataFullCopy.molecules:
        if not molecule in G.chosenMolecules:
            unusedMolecules.append(molecule)
            
    #Make into a global variable for future use.     
    G.unusedMolecules = unusedMolecules
    
    #Export current Reference Data  
    #Reference data is trimmed prior to this function
    ExportXYYYData(G.referenceFileName, ReferenceData.provided_reference_intensities, ReferenceData.molecules, abscissaHeader = 'M/Z')
    
    #Export current Experimental Data
    #Experimental data is trimmed prior to this function, but it still needs to be exported  
    ExportXYYYData(G.dataToAnalyzeFileName, ExperimentData.workingData, ExperimentData.mass_fragment_numbers,
              abscissaHeader = ExperimentData.abscissaHeader, dataType = 'preProcessed', rowIndex = ExperimentData.times)
   
    #export reference data for next iteration
    if G.iterationNumber == 1: #first iteration files aren't in standard locations
        DataFunctions.TrimReferenceFileByMolecules(unusedMolecules, "..\\%s" %G.oldReferenceFileName, unusedReferenceFileName = G.nextRefFileName)
    else: #not first iteration
    #generate unused reference data
        DataFunctions.TrimReferenceFileByMolecules(unusedMolecules, G.oldReferenceFileName, unusedReferenceFileName = G.nextRefFileName)
    
    return None
    #implied returns: G.unusedMolecules

def IterativeAnalysisPostProcessing(ExperimentData, simulateddata, mass_fragment_numbers,ExperimentDataFullCopy, times, concdata, molecules):
    #implied arguments: G.iterationSuffix, G.nextRefFileName, G.nextExpFileName, G.iterativeAnalysis, G.unusedMolecules, G.iterationNumber
    #remove the signals that have already been solved for from the data set
    ExperimentData.workingData = DataFunctions.RemoveSignals(ExperimentDataFullCopy.workingData, ExperimentDataFullCopy.mass_fragment_numbers, simulateddata, mass_fragment_numbers)
    
    #Export the remaining experimental signals
    ExportXYYYData(G.nextExpFileName, ExperimentDataFullCopy.workingData, ExperimentDataFullCopy.mass_fragment_numbers, abscissaHeader = ExperimentData.abscissaHeader, dataType = 'Experiment', rowIndex = ExperimentData.times)
    
    #update the suffix number and create the next user input file
    G.iterationSuffix = '_iter_%s' %str(G.iterationNumber + 1)
    nextUserInputFileName = 'UserInput%s.py' %G.iterationSuffix 
    
    #revert to the parent directory
    os.chdir('..')
    #create the next iteration directory and change the cwd into the next iteration directory
    IterationDirectoryPreparation(G.iterativeAnalysis, G.iterationNumber, iterate = True)
     
    #save the new file name for the next user input file 
    G.dataToAnalyzeFileName = G.nextExpFileName 
    G.referenceFileName = G.nextRefFileName
    #updating the selected molecules for the next user input file
    G.chosenMolecules = G.unusedMolecules
    
    #export the user input specifications 
    ExportUserInputFile(nextUserInputFileName)
    #append the variable list to the user input file
    AppendList("__var_list__", G.__var_list__, nextUserInputFileName, 5)
    
    if G.iterativeAnalysis == True:
        iterationDirectoryName = '_iter_%s' %(str(G.iterationNumber - 1))
    if not G.iterativeAnalysis == True:
        iterationDirectoryName = '%s_iter_%s' %(G.iterativeAnalysis, str(G.iterationNumber - 1))
    #copy the experimental signals to the next iteration
    shutil.copy("..\%s\%s" %(iterationDirectoryName, G.dataToAnalyzeFileName), os.getcwd())
    #copy the next reference file from the previous iteration folder to the next iteration folder
    shutil.copy("..\%s\%s" %(iterationDirectoryName, G.referenceFileName), os.getcwd())
    
    #returning to the parent directory
    os.chdir('..')
    
    #Adding the Additional concentrations to the overall concentration results
    moleculeConcLabels = ['%s Concentration Relative to CO' % molecule for molecule in molecules] 
    DataFunctions.AppendColumnsToCSV(G.TotalConcentrationsOutputName, concdata, moleculeConcLabels, times, ["Time"])
    
    return None
     #implied returns: G.referenceFileName, G.dataToAnalyzeFileName, G.nextRefFileName, G.chosenMolecules, G.iterationSuffix
###############################################################################
#########################  Functions to read data files #######################
###############################################################################
#These functions read in the experimental data file and the reference file. The
#returned variables can then be used to initialize the respective classes.

def readDataFile(dataToAnalyzeFileName):

 #read the csv file into a dataframe.  dataFrame means "dataframe" and is a pandas object.
    dataFrame = pandas.read_csv('%s' %dataToAnalyzeFileName, header=None)
    ''' generate mass fragment list'''
    #select only the 2nd row down, all columns except for the first. 
		#"iloc" is a pandas dataframe function. All it does is select a portion of the data.
    dfmass = dataFrame.iloc[1][1:]
    #convert to matrix
    masses = dfmass.as_matrix()
    #sort through the matrix and remove labels
    for i in range(0,len(masses)):
        masses[i] = masses[i].replace('mass','')
        masses[i] = masses[i].replace('m','')
    #convert the matrix to integers 
    mass_fragment_numbers = masses.astype(numpy.float32)
            
    '''generate time list'''
    # set abscissa header (time or Temp, etc.)
    abscissaHeader = dataFrame.iloc[1][0]
    #select column of times
    dftimes = dataFrame.iloc[2:][0]
    #convert to matrix
    times = dftimes.as_matrix()
    #save as class object with type float
    times = times.astype(numpy.float)
    #if the user wants to analyze one point, the data is doubled in length
    #to prevent future index problems
    if len(times) == 1:
        times = numpy.append(times,times)
        times[1] = times[0]*1.1
   
    '''generate collected data'''
    #select matrix of raw signals
    dfcollected = dataFrame.iloc[2:,1:]
    #convert to matrix
    collected = dfcollected.as_matrix()
    #save as class object with type float
    rawCollectedData = collected.astype(numpy.float)
    #if the user wants to analyze one point, the data is doubled in length
    #to prevent future index problems
    if len(rawCollectedData) == 1:
        rawCollectedData = numpy.vstack((rawCollectedData,rawCollectedData))

        
    return mass_fragment_numbers, abscissaHeader, times, rawCollectedData, dataToAnalyzeFileName

#readReferenceFile is a helper function that reads the reference file in a certain form and returns the
#variables and data that are used to initialize the class. It can read files both in XYYY and XYXY form.
def readReferenceFile(referenceFileName, form):        
     #This function converts the XYXY data to an XYYY format
    def FromXYXYtoXYYY(provided_reference_intensities):
        masslists = [] #future lists must be must empty here to append in the for loops
        relativeintensitieslists = [] #future list
        #this loops gathers all the mass fragment numbers for each molecule in one list of arrays, while a second
        #list is made, gathering the relative intensities so that they were indexed the same as their mass fragment
        #numbers in the other list
        #this for loop grabs stuff from the reference array, whose orientation and identity is shown in the flow chart arrays document
        for referenceBy2Index in range(0,len(provided_reference_intensities[0,:]),2):#array-indexed for loop, only gets every other value, as half the indexes are mass lists, and the other half are relative intensity
            masslists.append(provided_reference_intensities[:,referenceBy2Index])#these are lists of arrays
            relativeintensitieslists.append(provided_reference_intensities[:,referenceBy2Index+1])#the relative intensities are after every counter, so there is a +1 (it is array indexed so since the first column is a mass list all the +1's are relative intensities)
        masslist = [] #future list
        #This for loop gets all of the mass fragments from the first index of the list, basically by not adding the 
        #'nan's or empty spaces after the numbers
        for referenceIndex in range(len(provided_reference_intensities[:,0])): #array-indexed for loop
            if str(masslists[0][referenceIndex]) != 'nan': #we do not want nan's in our array, the genfromtxt function calls empty boxes in excel (might be in .csv as well)'nan'.
                masslist.append(masslists[0][referenceIndex])
        #this second nested for loop gathers all the other mass fragment numbers that have not already been added to
        #the masslist, basically obtaining all the masses in the reference data and then after the loop they are sorted
        #using .sort, then an empty array of zeros is made to accommodate the output array
        for masslistIndex in range(1,len(masslists)):#array-indexed for loop, starts at one because it's checking against all arrays besides itself
            for referenceIndex in range(len(provided_reference_intensities[:,0])):#array-indexed for loop
                if str(masslists[masslistIndex][referenceIndex]) != 'nan':
                    if sum(masslists[masslistIndex][referenceIndex] == numpy.array(masslist)) == 0:#if the value being looked at is not equal to anything in our masslist already
                        masslist.append(masslists[masslistIndex][referenceIndex])
        masslist.sort()#puts the list in order
        reference_holder = numpy.zeros([len(masslist),len(provided_reference_intensities[0,:])/2+1])#makes an array that is full of zeros to hold our future reference array
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
        provided_reference_intensities = reference_holder
        return provided_reference_intensities
    
     #read the csv file into a dataframe
    dataFrame = pandas.read_csv('%s' %referenceFileName, header = None)
    
    if form == 'xyyy':
        ''' generate reference matrix'''
        #remove top 4 rows
        dfreference = dataFrame.iloc[4:][:]
        #convert to matrix
        reference = dfreference.as_matrix()
        #convert the matrix to floats
        provided_reference_intensities = reference.astype(numpy.float)
        #clear rows of zeros
        provided_reference_intensities=DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(provided_reference_intensities,startingRowIndex=1)
    
        '''generate electron number list'''
        #select row of electron numbers
        dfelectronnumbers = dataFrame.iloc[2][1:]
        #convert to matrix
        electronnumbers = dfelectronnumbers.as_matrix()
        #save as class object with type int
        electronnumbers = electronnumbers.astype(numpy.int32)
   
        '''generate list of molecule names'''
        #select row of names
        dfmolecules = dataFrame.iloc[1][1:]
        #convert to matrix
        molecules = dfmolecules.as_matrix()
        #save as class object with type string
        molecules = molecules.astype(numpy.str)
        
        '''generate list of molecular weights'''
        #select row of names
        dfmolecularWeights = dataFrame.iloc[3][1:]
        #convert to matrix
        molecularWeights = dfmolecularWeights.as_matrix()
        #save as class object with type float
        molecularWeights = molecularWeights.astype(numpy.float)
        
        '''generate list of source information'''
        #select row of names
        dfsourceInfo = dataFrame.iloc[0][1:]
        #convert to matrix
        sourceInfo = dfsourceInfo.as_matrix()
        #save as class object with type string
        sourceInfo = sourceInfo.astype(numpy.str)
        
        '''generate list of massfragments monitored '''
        mass_fragment_numbers_monitored = provided_reference_intensities[:,0]
        
    elif form == 'xyxy':
        '''generate reference matrix'''
        #remove top 4 rows
        dfreference = dataFrame.iloc[4:][:]
        #convert to matrix
        reference = dfreference.as_matrix()
        #convert the matrix to floats 
        provided_reference_intensities = reference.astype(numpy.float)
        #convert reference from XYXY to XYYY
        provided_reference_intensities=FromXYXYtoXYYY(provided_reference_intensities)
        #clear rows of zeros
        provided_reference_intensities=DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(provided_reference_intensities,startingRowIndex=1)
        
        '''generate electron numbers list'''
        #create data frame of electron numbers
        dfelectronnumbers = dataFrame.iloc[2,1::2]
        #convert to matrix
        electronnumbers = dfelectronnumbers.as_matrix()
        #save as class object with type int
        electronnumbers = electronnumbers.astype(numpy.int32)
        
        '''generate list of molecule names'''
        #select matrix of names
        dfmolecules = dataFrame.iloc[1,1::2]
        #convert to matrix
        molecules = dfmolecules.as_matrix()
        #save as class object with type string
        molecules = molecules.astype(numpy.str)
        
        '''generate list of molecular weights'''
        #select row of names
        dfmolecularWeights = dataFrame.iloc[3][1::2]
        #convert to matrix
        molecularWeights = dfmolecularWeights.as_matrix()
        #save as class object with type float
        molecularWeights = molecularWeights.astype(numpy.float)
        
        '''generate list of source information'''
        #select row of names
        dfsourceInfo = dataFrame.iloc[0][1::2]
        #convert to matrix
        sourceInfo = dfsourceInfo.as_matrix()
        #save as class object with type string
        sourceInfo = sourceInfo.astype(numpy.str)

        '''generate list of massfragments monitored '''
        mass_fragment_numbers_monitored = provided_reference_intensities[:,0]
        
    return provided_reference_intensities, electronnumbers, molecules, molecularWeights, sourceInfo, mass_fragment_numbers_monitored, referenceFileName, form      

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
		
    def __init__(self, mass_fragment_numbers, abscissaHeader, times, rawCollectedData, dataToAnalyzeFileName=None):
        
        [self.mass_fragment_numbers, self.abscissaHeader, self.times, self.rawCollectedData, self.dataToAnalyzeFileName]=[mass_fragment_numbers, abscissaHeader, times, rawCollectedData, dataToAnalyzeFileName]
        
        
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
        
        if G.ExportAtEachStep == 'yes':
            #record data of experiment
            self.dataToExport.append(self.workingData.copy())
            #record times from the data of the experiment
            self.experimentTimes.append(self.times.copy())
            
    def ExportMSData(self):
        print("\n Collected Export List:")
        for savePoint in range(len(self.runTimeAtExport)):
            print(self.labelToExport[savePoint])
            print(self.runTimeAtExport[savePoint])
            if G.ExportAtEachStep == 'yes':
                #inserting the data for a particular savePoint
                filename = 'Exported%s%s.csv'%(savePoint, self.labelToExport[savePoint]) #FIXME: for DataSmoother, and some others, the debug output has a "Time" header but the time is not exported.
                data = self.dataToExport[savePoint]
                abscissa = self.experimentTimes[savePoint]
                colIndex = ['m%s'% int(y) for y in self.mass_fragment_numbers]
                DataFunctions.MSDataWriterXYYY(filename, data, abscissa, colIndex, self.abscissaHeader)
                                        
class MSReference (object):
    def __init__(self, provided_reference_intensities, electronnumbers, molecules, molecularWeights, sourceInfo, mass_fragment_numbers_monitored, referenceFileName=None, form=None):
        [self.provided_reference_intensities, self.electronnumbers, self.molecules, self.molecularWeights, self.sourceInfo, self.mass_fragment_numbers_monitored, self.referenceFileName, self.form]=[provided_reference_intensities, electronnumbers, molecules, molecularWeights, sourceInfo, mass_fragment_numbers_monitored, referenceFileName, form]
            
        '''Initializing Export Collector Variables'''
        #start the timer function
        self.previousTime = timeit.default_timer()
        #initalize debugging lists
        #These lists are appended in parallel so a variable of the same index from each list will be related
        self.runTimeAtExport = []
        self.labelToExport = []
        self.dataToExport = []
        #self.experimentTimes = []
            
    def ExportCollector(self, callingFunction, use_provided_reference_intensities = False):
        #record current time
        currentTime = timeit.default_timer()
        #add net time to list of run times
        self.runTimeAtExport.append(currentTime - self.previousTime)
        #record current time for next function's use
        self.previousTime = currentTime
        #add the name of the calling function to mark its use
        self.labelToExport.append(callingFunction) 
        
        if G.ExportAtEachStep == 'yes':
            #record data of experiment
            if use_provided_reference_intensities:
                self.dataToExport.append(self.provided_reference_intensities.copy())
            elif not use_provided_reference_intensities:
                self.dataToExport.append(self.standardized_reference_intensities.copy())
            
    def ExportFragmentationPatterns(self):
        print("\n Reference Debugging List:")
        for savePoint in range(len(self.runTimeAtExport)):
            print(self.labelToExport[savePoint])
            print(self.runTimeAtExport[savePoint])
            if G.ExportAtEachStep == 'yes':
                #inserting the data for a particular savePoint
                filename = 'Exported%s%s.csv'%(savePoint, self.labelToExport[savePoint])
                data = self.dataToExport[savePoint]
                colIndex = ['%s'% y for y in self.molecules]
                ExportXYYYData(filename,data,colIndex, fileSuffix = G.iterationSuffix)

    # This class function removes all rows of zeros from
    # the XYYY sorted reference data.
    def ClearZeroRows(self):
        zeros = True
        placeHolder = 0
        for i in self.provided_reference_intensities:
            for j in i[1:]:
                if j != 0:
                    zeros = False
            if(zeros == True):
                self.provided_reference_intensities = numpy.delete(self.provided_reference_intensities, placeHolder,0)
                placeHolder = placeHolder -1
            zeros= True    
            placeHolder = placeHolder + 1
            
#This class function converts the XYXY data to an XYYY format
    def FromXYXYtoXYYY(self):
        masslists = [] #future lists must be must empty here to append in the for loops
        relativeintensitieslists = [] #future list
        #this loops gathers all the mass fragment numbers for each molecule in one list of arrays, while a second
        #list is made, gathering the relative intensities so that they were indexed the same as their mass fragment
        #numbers in the other list
        #this for loop grabs stuff from the reference array, whose orientation and identity is shown in the flow chart arrays document
        for referenceBy2Index in range(0,len(self.provided_reference_intensities[0,:]),2):#array-indexed for loop, only gets every other value, as half the indexes are mass lists, and the other half are relative intensity
            masslists.append(self.provided_reference_intensities[:,referenceBy2Index])#these are lists of arrays
            relativeintensitieslists.append(self.provided_reference_intensities[:,referenceBy2Index+1])#the relative intensities are after every counter, so there is a +1 (it is array indexed so since the first column is a mass list all the +1's are relative intensities)
        masslist = [] #future list
        #This for loop gets all of the mass fragments from the first index of the list, basically by not adding the 
        #'nan's or empty spaces after the numbers
        for referenceIndex in range(len(self.provided_reference_intensities[:,0])): #array-indexed for loop
            if str(masslists[0][referenceIndex]) != 'nan': #we do not want nan's in our array, the genfromtxt function calls empty boxes in excel (might be in .csv as well)'nan'.
                masslist.append(masslists[0][referenceIndex])
        #this second nested for loop gathers all the other mass fragment numbers that have not already been added to
        #the masslist, basically obtaining all the masses in the reference data and then after the loop they are sorted
        #using .sort, then an empty array of zeros is made to accommodate the output array
        for masslistIndex in range(1,len(masslists)):#array-indexed for loop, starts at one because it's checking against all arrays besides itself
            for referenceIndex in range(len(self.provided_reference_intensities[:,0])):#array-indexed for loop
                if str(masslists[masslistIndex][referenceIndex]) != 'nan':
                    if sum(masslists[masslistIndex][referenceIndex] == numpy.array(masslist)) == 0:#if the value being looked at is not equal to anything in our masslist already
                        masslist.append(masslists[masslistIndex][referenceIndex])
        masslist.sort()#puts the list in order
        reference_holder = numpy.zeros([len(masslist),len(self.provided_reference_intensities[0,:])/2+1])#makes an array that is full of zeros to hold our future reference array
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
        self.provided_reference_intensities = reference_holder
############################################################################################################################################
###############################################Algorithm Part 2: Analysing the Processed Data###############################################
############################################################################################################################################
    
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
    moleculenum = len(reciprocal_matching_correction_values[0,:])
    row_number = len(reciprocal_matching_correction_values[:,0])
    import itertools 
    combinations = list(itertools.combinations(list(range(row_number)),moleculenum)) 
    if combinations == []:#This function will not work without enough mass fragments, so the user must know the problem
        print('****************************************')
        print('Not enough matching mass fragments input')
        print("This means that at some point in the analysis, there were not enough masses in the reference file to apply the inverse method. It could mean you have too many overlapping masses for the molecules you are trying to resolve.  You can get around this by using the '#//Reference Mass Fragmentation Threshold//' feature to exclude tiny fragementation peaks. This would be done by setting the value to 'yes' for  applyReferenceMassFragmentsThresholds feature with referenceMassFragmentFilterThreshold, such as referenceMassFragmentFilterThreshold = 5.0 .  Alternatively, to be more targeted, if you know *which* fragmentation patterns could be overlapping, you could set those minor fragments to 0 in your reference pattern csv file. TODO: Print out the relevant masses here. This requires keeping track of when they are selected prior to combination maker, and possibly passing them as an additional argument.")
        print('****************************************')
    combinations_len = len(combinations) 
    correctionarray = numpy.zeros([1,moleculenum])
    intensityarray = numpy.zeros([1,moleculenum])
    rawsignalarray = numpy.zeros([1,1])
    reciprocal_correction_Values_list = []
    intensitylist = []
    rawsignallist = [] 
    massfragrow = 'yuh,'
    massfraglist = []
    #this loop gets two lists that contain al the arrays of possible answers and inverse correction values,first rows are made,
    #then they are stacked into an array, this array then, has the first row of zeros deleted, then it is added to an indice
    #of a list where it will later be accessed by the combination solver. three of these lists are output for signals,corrections
    #and relative intensities
    for combinationnum in range(combinations_len): #array-indexed for loop
        for moleculecounter in range(moleculenum):    #array-indexed for loop
            correctionrow = reciprocal_matching_correction_values[combinations[combinationnum][moleculecounter],:] 
            intensityrow = monitored_reference_intensities[combinations[combinationnum][moleculecounter],:]
            rawsignalrow = rawsignalsarrayline[combinations[combinationnum][moleculecounter],:]
            massfragrow = massfragrow + str(mass_fragment_numbers[combinations[combinationnum][moleculecounter]]) + ','
            correctionarray = numpy.vstack([correctionarray,correctionrow]) 
            intensityarray = numpy.vstack([intensityarray,intensityrow])
            rawsignalarray = numpy.vstack([rawsignalarray,rawsignalrow])
            if moleculecounter == moleculenum-1:#the end of the nested loop: the rows just made are entered into arrays
                correctionarray = numpy.delete(correctionarray,(0),axis=0)
                intensityarray = numpy.delete(intensityarray,(0),axis=0)
                rawsignalarray = numpy.delete(rawsignalarray,(0),axis=0)
                massfragrow1 = massfragrow.split(',',1)[1]
                massfraglist.append(massfragrow1)
                massfragrow = 'yuh,'
                reciprocal_correction_Values_list.append(correctionarray)
                intensitylist.append(intensityarray)
                rawsignallist.append(rawsignalarray)
                correctionarray = numpy.zeros([1,moleculenum])
                intensityarray = numpy.zeros([1,moleculenum])
                rawsignalarray = numpy.zeros([1,1])
    combinations_len = len(combinations)
    return [combinations_len,rawsignallist,reciprocal_correction_Values_list,intensitylist,massfraglist]


#This function simply solves each of the combinations, drawing the respective values out of the lists and uses numpy.linalg
def CombinationSolver(combinations_len,rawsignallist,reciprocal_correction_Values_list,molecules,massfraglist):
    compositions = []
    for combinationcounter in range (combinations_len):  #array-indexed for loop
        if numpy.linalg.det(reciprocal_correction_Values_list[combinationcounter]) != 0:#if the determinant is zero, then doing the linalg.solve function will stop the entire script- so you must use this method
            solutions = numpy.linalg.solve(reciprocal_correction_Values_list[combinationcounter], rawsignallist[combinationcounter])
            composition = solutions 
            compositions.append(composition)
    return[compositions]
    

#compresses the data into an array of averages and standard deviations- then prints the results
def DataCompressor(signals,molecules,type):
    moleculenum = len(molecules)
    averagegroup = []
    stddevgroup = []
    average = []
    stddev = []
    for moleculecounter in range(moleculenum): #this part of the code is new for this version (3) and it gets avg and stddev
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
    [combinations_len,rawsignallist,reciprocal_correction_Values_list,intensitylist,massfraglist] = CombinationMaker (reciprocal_matching_correction_values,rawsignalsarrayline,monitored_reference_intensities,mass_fragment_numbers)
    [compositions] = CombinationSolver (combinations_len,rawsignallist,reciprocal_correction_Values_list,molecules,massfraglist)
    [averagecomp,stddevcomp] = DataCompressor(compositions,molecules,'composition')
    return averagecomp
    
    

#this function finds the significance of a specified value in an array to the array as a whole
def IndElemSignificanceCalculator(dataArray, column, moleculesLikelihood):
    length = len(dataArray)
    #variable to hold the terms in the summation
    numbers = []
    #for each value in the array
    for moleculecounter in range(length):
        if dataArray[moleculecounter] != 0: #if the value is zero, then the final value will be zero as well: it doesn't have to be changed
            #calculates the unweighted ratio of each value, scaled by the likelihood of that molecule 
            summationTerm = abs((moleculesLikelihood[moleculecounter]*dataArray[moleculecounter])**float(-1)*(dataArray[column]*moleculesLikelihood[column]-1))
            numbers.append(summationTerm)
    #the following line can be replace with code such as "significance = (sum(numbers)**SumCoeffient)*(array[column]**ValueCoefficent)"
    # if you would like to add coefficents to increase or decrease the weighting of each term
    significance = sum(numbers)*dataArray[column]*moleculesLikelihood[column]
    return significance

#This function compiles a list of the significances of each row to a particular column 
def ElemSignificanceCalculator(anArray,columncounter, moleculesLikelihood):
    #find the number of rows
    row_num = len(anArray)
    #empty list to store values
    sigValuesList = []
    #for each row...
    for rowcounter in range(row_num):
        # the "Significance" of that row to the column is calculated
        sigValue = IndElemSignificanceCalculator(anArray[rowcounter], columncounter, moleculesLikelihood)
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
        if overallOrder == []:
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
#in that row and that molecules own value
def DistinguishedArrayChooser(refMassFrags,correctionValues,rawSignals,moleculeLikelihoods,sensitivityValues):
    #the shape of the referenceData is found 
    row_num = len(refMassFrags[:,0])
    column_num = len(refMassFrags[0,:])
    
    #The acceptable threshold is determined by the SensitivityValue function
    sensitivityValues = ListLengthChecker(sensitivityValues, column_num, 1)
   
    #the moleculesLikelihood is corrected if it wasn't entered by the use.
    moleculeLikelihoods = ListLengthChecker(moleculeLikelihoods, column_num, 1)
    
    #all values below the specified relative intensity must be set the minThreshold value
    #This is because a subfunction attempts to divide by each value
    for columncounter in range(column_num):
        for rowcounter in range(row_num):
            if refMassFrags[rowcounter,columncounter] < sensitivityValues[columncounter]: 
                refMassFrags[rowcounter,columncounter] = 0 #sensitivityThresholdValue[0]
                
    #The correct order of the needed rows is determined by this function
    order = ImportantAbscissaIdentifier(refMassFrags,moleculeLikelihoods)
    
    #empty lists to store results i.e. shortened arrays
    shortRefMassFrags = []
    shortCorrectionValues = []
    shortRawSignals = []
    
    #add the correct row to each list
    for row_num in order:
        shortRefMassFrags.append(refMassFrags[row_num])
        shortCorrectionValues.append(correctionValues[row_num])
        shortRawSignals.append(rawSignals[row_num])
              
    #This section stacks the chosen rows from lists into arrays
    shortRefMassFrags = numpy.asarray(shortRefMassFrags)
    shortCorrectionValues = numpy.asarray(shortCorrectionValues)
    shortRawSignals = numpy.asarray(shortRawSignals)
    
    #finding the size of the new array
    row_num = len(shortRefMassFrags[:,0])
    column_num = len(shortRefMassFrags[0,:])
    #This section replaces the minThreshold's that were chosen as threshold values with 0s
    for columncounter in range(column_num):
        for rowcounter in range(row_num):
            if shortRefMassFrags[rowcounter,columncounter] <= sensitivityValues[columncounter]: 
                shortRefMassFrags[rowcounter,columncounter] = 0
                
    #The shortened arrays are finally returned to the Inverse Method solver                
    return shortRefMassFrags,shortCorrectionValues,shortRawSignals
    
    
#this function takes the data from important abscissa identifier and 
def InverseMethodDistinguished(monitored_reference_intensities,reciprocal_matching_correction_values,rawsignalsarrayline):
    monitored_reference_intensities,reciprocal_matching_correction_values,rawsignalsarrayline = DistinguishedArrayChooser (monitored_reference_intensities,reciprocal_matching_correction_values,rawsignalsarrayline, G.moleculeLikelihoods,G.sensitivityValues)
    if numpy.linalg.det(reciprocal_matching_correction_values) != 0:#only solves if determinant is not equal to zero
        solutions = numpy.linalg.solve(reciprocal_matching_correction_values,rawsignalsarrayline)
    else:
        print('The Array Chosen is Singular')
        solutions = numpy.zeros(len(rawsignalsarrayline)) # the solutions are made into all zeros if the chosen array is singular
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
def SignalSimulation(sampleparameterpoints,*list):
    rawsignalsarrayline = list[0]#the item in the list is the raw signal arrayline
    reciprocal_matching_correction_values = list[1] #the second item is the matching correction values
    objectiveFunctionType = list[2]#the input argument contains both the objectiveFunctionType and the two arrays: raw signals array line and matching correction values
    xyyData = numpy.zeros([3,len(rawsignalsarrayline)]) #a three line array is made that will be entered into the function below
    xyyData[1:2,:] = numpy.hstack(rawsignalsarrayline) #this second line is the raw signals
    xyyData[2:3,:] = numpy.hstack(numpy.array(numpy.matrix(reciprocal_matching_correction_values)*numpy.matrix(numpy.vstack(sampleparameterpoints)))) #the third row is the calculated signals
    objectiveFunctionDictionary = ObjectiveFunctionGenerator(xyyData,0)
    if objectiveFunctionType == 'weightedSAR': #based on the choice given the output will be chosen from this called functions dictionary
        objective_function = objectiveFunctionDictionary['weightedSAR']
    if objectiveFunctionType == 'weightedSSR':
        objective_function = objectiveFunctionDictionary['weightedSSR']
    if objectiveFunctionType == 'ssr':
        objective_function = objectiveFunctionDictionary['ssr']
    if objectiveFunctionType == 'sar':
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
def BruteForce(molecules,specifications,reciprocal_matching_correction_values,rawsignalsarrayline,objectiveFunctionType,maxPermutations=100001):
    list = []#a list is made in order to use inside the forward problem since we are not defining anything globally
    list.append(rawsignalsarrayline)
    list.append(reciprocal_matching_correction_values)
    list.append(objectiveFunctionType)
    from scipy import optimize
    summation = numpy.zeros(len(specifications))
    
    quotient = numpy.zeros(len(specifications))
    
    product = 1
    for specificationnumber in range(len(specifications)):
        summation[specificationnumber] = specifications[specificationnumber][1]-specifications[specificationnumber][0]
        quotient[specificationnumber] = summation[specificationnumber]/float(specifications[specificationnumber][2])
        product = quotient[specificationnumber]*product

    if product < maxPermutations:
        resbrute = optimize.brute(SignalSimulation, specifications, args = list, full_output=True, finish=None)#calls our forward function SignalSimulation

    else:
        print('Error: Too many Permutations')
        sys.exit()
    answers = resbrute[0]

    return answers
    
    
#this function is a path is sequential linear subtraction, which can be used alongside the inverse
#method or as opposed to it. Either way, it only produces one set of values, so it has no need for the 
#data compressor function and starts right after the correction values are obtained
def SLSUniqueFragments(molecules,monitored_reference_intensities,reciprocal_matching_correction_values,rawsignalsarrayline, timeIndex, time):

    # This is creating a local copy of 'monitored_reference_intensities' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reference_intensities_SLS = copy.deepcopy(monitored_reference_intensities)

    # This is creating a local copy of 'reciprocal_matching_correction_values' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_correction_values_SLS = copy.deepcopy(reciprocal_matching_correction_values)

    # This is creating a local copy of 'rawsignalsarrayline' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_rawsignals_SLS = copy.deepcopy(rawsignalsarrayline)

    # This is creating a local copy of 'molecules' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_molecules_SLS = copy.deepcopy(molecules)
    
    selectedvariable = []
    rawsignals1 = []
    solutions1 = numpy.zeros([1,len(remaining_correction_values_SLS[0,:])])
    solutions = solutions1[0]
    moleculenum = len(remaining_correction_values_SLS[0,:]) 
    row_num = len(remaining_correction_values_SLS[:,0])
    stable_row_num = row_num 
    stable_molecule_num = moleculenum 
    stablecorrection_values = remaining_correction_values_SLS
    molecules_unedited = copy.deepcopy(molecules)
    usedmolecules = numpy.zeros(len(remaining_molecules_SLS))
    reminder2 = 0 
    #These nested for loops make the whole function run - they go through all the rows
    #and all the columns, but they do this as many times as there are rows, so that all
    #the values that can be found using this method will be. The values for row_num and
    #moleculenum are reset every cycle in order to run the other loops right
    for molNumIndex in range(moleculenum):#array-indexed for loop
        row_num = len(remaining_correction_values_SLS[:,0])
        moleculenum = len(remaining_correction_values_SLS[0,:])
        reminder = 1  
        for rowcounter in range(row_num):#array-indexed for loop
            for moleculecounter in range(moleculenum):#array-indexed for loop
                if reminder == 1:#if the reminder is zero, the loop can't go through, so the solution is found, the the loop restarts with another index.
                    if remaining_correction_values_SLS[rowcounter,moleculecounter] != 0:#If the value in the correction_values is not zero, it is kept
                        selectedvariable.append(remaining_correction_values_SLS[rowcounter,moleculecounter]) 
                        rawsignals1.append(remaining_rawsignals_SLS[rowcounter]) 
                        whichmolecule = moleculecounter 
                    #once all the molecules have been looked at for a certain row of correction values, if there is
                    #a row with only one mass fragment, then that row is kept, so that it can be used to find the 
                    #compositions of the molecules
                    if moleculecounter == moleculenum-1: #if you get to the end of the column
                        if len(selectedvariable) == 1: #if there was only one value that was non zero
                            for moleculechooser in range(stable_molecule_num): #array-indexed for loop
                                for rowchooser in range(stable_row_num):#array-indexed for loop
                                    if stablecorrection_values[rowchooser,moleculechooser] == selectedvariable:#finds index of stable correction values where the selected variable is
                                        whichmolecule2 = moleculechooser 
                            #this for loop multiplies the relative amount that was just found by the correction values in order to determine the amount of signal that the 
                            #certain molecule in question accounted for so that right after the for loop, this signal could be erased from the full signal, and thus the other
                            #molecules relative amounts could be solved for
                            rawsignalsused = numpy.zeros([row_num,1])
                            for rowcounter2 in range(row_num):#array-indexed for loop
                                if stablecorrection_values[rowcounter2,whichmolecule2]!= 0: #used to find the raw signals that are made b the molecules that are being deleted
                                    rawsignalsused[rowcounter2] = stablecorrection_values[rowcounter2,whichmolecule2] * (float(rawsignals1[0][0]))/float(selectedvariable[0])
                            solutions[whichmolecule2] = ((float(rawsignals1[0][0]))/float(selectedvariable[0]))
                            remaining_correction_values_SLS = numpy.delete(remaining_correction_values_SLS,(whichmolecule),axis = 1)
                            remaining_reference_intensities_SLS = numpy.delete(remaining_reference_intensities_SLS,(whichmolecule),axis = 1)
                            remaining_molecules_SLS = numpy.delete(remaining_molecules_SLS,(whichmolecule))
                            usedmolecules[whichmolecule2] = 1
                            remaining_rawsignals_SLS = remaining_rawsignals_SLS - rawsignalsused
                            reminder = reminder-1 #if an answer is found the reminder is changed so that the loop has to be reset before it can search again
                            reminder2 = reminder2 +1 #know the number of solutions
                        if reminder2 == len(solutions): #if all of the solutions have been found
                            reminder = 2
                        selectedvariable = []
                        rawsignals1 = []

        # This block of code is a printing statement to show the user what order the molecules are being solved in
        # This is a csv file so should be delimited with commas
        if G.SLSUniquePrint == 'yes':
            with open(G.SLSUniqueExport, 'a') as f:
                f.write('%s,' %timeIndex)
                f.write('%s,' %time)
                for x in range(len(usedmolecules)):
                    f.write('%s,' %usedmolecules[x])
                f.write("UsedMolecules \n")
        
    if remaining_correction_values_SLS.size > 0:#if there are correction values left (i.e. not all the solutions have been found)
        #this for loop is used to delete any rows entirely composed of zeros, since the molecules percentages are found
        #based on the column index, the rest of the loop still works fine, this is just so that there are not singular
        #matrices made by combination maker. (originally there would be up to hundreds of simply useless matrices solved)
        place_holder = 0
        for correctionIndex in range(len(remaining_correction_values_SLS[:,0])):#array-indexed for loop
            if any(remaining_correction_values_SLS[correctionIndex-place_holder,:]) == 0:#all rows of only zeros are deleted
                remaining_correction_values_SLS = numpy.delete(remaining_correction_values_SLS,correctionIndex-place_holder,axis = 0)
                remaining_rawsignals_SLS = numpy.delete(remaining_rawsignals_SLS,correctionIndex-place_holder,axis = 0)
                remaining_reference_intensities_SLS = numpy.delete(remaining_reference_intensities_SLS,correctionIndex-place_holder,axis = 0)
                place_holder = place_holder + 1#since the arrays are being deleted, this keeps the indexing correct
    if sum(usedmolecules) == 0:#if none of the solutions have been found
        solutions = []
        usedmolecules = []

    return [remaining_molecules_SLS,remaining_reference_intensities_SLS,remaining_correction_values_SLS,remaining_rawsignals_SLS,solutions,molecules_unedited,usedmolecules]
    
#this sls method cuts smaller, solvable arrays out of the large array and uses numpy.linalg.solve to find the signals
#relative to Co for those molecules, then via the SLSMethod function below the function works in conjunction with the 
#sls unique fragments in order to solve all the molecules in the array that can be solved via arrays of ize one, two or
#three (ultimately, any size array can be solved via this method (if the method went to infinity rather than stopping at
#3) provided that the array were square- for, in fact, none of these method (save for the brute method) can solve for the 
#signals if the first array is not square
def SLSCommonFragments(reciprocal_matching_correction_values,rawsignalsarrayline,monitored_reference_intensities,molecules,scaledConcentrationsarray,molecules_unedited,conversionfactor,datafromcsv,DataRangeSpecifierlist,objectiveFunctionType,counterforspecifications,maxPermutations = 100001):

    #TODO: #FIXME: It seems 
    # like this function forces Brute towards the end, rather than checking 
    # what the user defined finisher is. This should be checked, eventually.

    # This is creating a local copy of the monitored reference intensities which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reference_intensities_SLS = copy.deepcopy(monitored_reference_intensities)

    # This is creating a local copy of the reciprocal_matching_correction_values which will become
    # truncated as the molecules are solved and masses are removed
    remaining_correction_values_SLS = copy.deepcopy(reciprocal_matching_correction_values)

    # This is creating a local copy of the rawsignalsarrayline which will become
    # truncated as the molecules are solved and masses are removed
    remaining_rawsignals_SLS = copy.deepcopy(rawsignalsarrayline)

    # This is creating a local copy of 'molecules' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_molecules_SLS = copy.deepcopy(molecules)
    #FIXME investigate why these two values are hardcoded and whether they can be made into user variables (Charles 180730)
    sizedecide = 5
    reminder3 = 0
    #molecules_unedited = remaining_molecules_SLS
    selectedvariablesholder = []
    selectedcorrectionsholder = []
    selectedmassfragsholder = []
    selectedrawsigsholder = []
    reminder = 1
    #this for loop goes through the whole nested loop as many times as the user would like, with the default setting being 3,
    #It resets the selected variable holder, which finds each row with the same number of variables 
    for sizecounter in range(2,sizedecide+1):#array-indexed for loop
    #if the array is empty here then the values are set to zero using the else statement
    #the original indexeers will cause an error because they require a two dimensional array
        if remaining_correction_values_SLS.shape[0] != 0:
            row_num = len(remaining_correction_values_SLS[:,0])
            moleculenum = len(remaining_correction_values_SLS[0,:])
        else:
            row_num = len(remaining_correction_values_SLS)
            moleculenum = len(remaining_correction_values_SLS)
        selectedvariables = []
        selectedcorrections = []
        selectedmassfrags = []
        selectedrawsigs = []
        useablecorrections = []
        useablemassfrags = []
        useablerawsigs = []
        #this loop looks at each row for the correction value array, resets the selected variable value
        for rowcounter in range(row_num):#array-indexed for loop
            selectedvariable = []
            reminder2 = 0
            reminder4 = 0
            #this loop looks at each molecule's correction value for each row, only if the reminder is equal to zero(will be
            #used later)
            for moleculecounter in range(moleculenum):#array-indexed for loop
                if reminder == 1:#will only run if no row has been chosen already
                    #the code below this only happens if the reminder is equal to one
                    #here a zero or a one is added to the selected variable, based on whether the correction value array has a
                    #value present or not, this way the different rows can be compared
                    if remaining_correction_values_SLS[rowcounter,moleculecounter] != 0:#chooses all values that are non zero- inputs ones
                        selectedvariable.append(1)  
                    else:
                        selectedvariable.append(0) 
                    #once the selected variable has been completed, this if statement is called
                    if moleculecounter == moleculenum-1: #at the end of the row
                        #this if statement gets all the rows with the needed 
                        if sum(selectedvariable) == sizecounter: #if the row has the correct number of variables
                            selectedvariables.append(selectedvariable)
                            selectedcorrections.append(remaining_correction_values_SLS[rowcounter])
                            selectedmassfrags.append(remaining_reference_intensities_SLS[rowcounter])
                            selectedrawsigs.append(remaining_rawsignals_SLS[rowcounter])
            #If the loop has finished making the selected variables list, then this second loop starts, and checks each of these
            #lists to see if this variable has enough rows with the same columns being used to be used as our array to solve
            if rowcounter == row_num - 1:#after all the selected variables have been put in arrays
                for rowcounter2 in range(row_num):#array-indexed for loop
                    selectedvariable = []
                    reminder2 = 0
                    for moleculecounter in range(moleculenum):#array-indexed for loop
                        if reminder == 1:#will only run if no row has been chosen already
                            if remaining_correction_values_SLS[rowcounter2,moleculecounter] != 0:#makes the selected variable again
                                selectedvariable.append(1)  
                            else:
                                selectedvariable.append(0)
                            if moleculecounter == moleculenum-1: #once it has been made
                            #this for loop iterates across all of the values that are the same length as the variable we're 
                            #checking for, it then makes an array using all of the rows that match (it will at least match with
                            #itself, because itself is in the selected variables)
                                for selectedIndex in range(len(selectedvariables)):#array-indexed for loop
                                    if selectedvariables[selectedIndex] == selectedvariable:#adds itself to the list of matching rows,checks other rows of same size
                                        reminder2 = reminder2 + 1
                                        if reminder2 == 1:#keeps the values that will be used in solving; if this row ends up being used
                                            useablecorrections = numpy.array(selectedcorrections[selectedIndex])
                                            useablemassfrags = numpy.array(selectedmassfrags[selectedIndex])
                                            useablerawsigs = numpy.array(selectedrawsigs[selectedIndex])
                                        else:
                                            useablecorrections = numpy.vstack([useablecorrections,selectedcorrections[selectedIndex]])
                                            useablemassfrags = numpy.vstack([useablemassfrags,selectedmassfrags[selectedIndex]])
                                            useablerawsigs = numpy.vstack([useablerawsigs,selectedrawsigs[selectedIndex]])
                                #if you are checking for arrays larger than size two then the size three array will check all the size 
                                #two arrays for lines that only contain values in the columns that have values for this size 3 row; 
                                #adding the matching rows to an array, the same one as before
                                if reminder4 == 0 and sizecounter > 2:#checks against past sizes for matching within one (has only values where the chosen row has values),reminder4
                                                                      #says that this has already been done so that the function will not add a bunch of extra rows and therefore
                                                                      #make the rest of the array unsolvable
                                    for rowIndex in range(len(selectedvariablesholder)):#this goes through each layer of chosen rows
                                        for colIndex in range(len(selectedvariablesholder[rowIndex])):#this goes through the rows of each chosen layer previous to the one being looked at
                                            #the discrepancy between the rows 
                                            if (len(selectedvariable) - sum(numpy.array(selectedvariable) == numpy.array(selectedvariablesholder[rowIndex][colIndex]))) <= sizecounter - (rowIndex + 2):#if they are good, then it adds all the required fields from that row
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
                                    for deleter in range(len(remaining_correction_values_SLS[:,0])):#array-indexed for loop
                                        place_holder2 = 0
                                        for checker in range(len(useablecorrections[:,0])):#array-indexed for loop
                                            if place_holder2 == 0:
                                                if all(remaining_correction_values_SLS[deleter-place_holder,:] == useablecorrections[checker,:]):#gets index of correction values to delete- the ones chosen
                                                    remaining_correction_values_SLS = numpy.delete(remaining_correction_values_SLS,deleter-place_holder,axis = 0)
                                                    remaining_rawsignals_SLS = numpy.delete(remaining_rawsignals_SLS,deleter-place_holder,axis = 0)
                                                    remaining_reference_intensities_SLS = numpy.delete(remaining_reference_intensities_SLS,deleter-place_holder,axis = 0)
                                                    place_holder = place_holder + 1
                                                    place_holder2 = place_holder2 + 1
                                    place_holder = 0
                                    place_holder2 = 0
                                    usedmolecules = selectedvariable
                                    for selectedcounter in range(len(selectedvariable)):#array-indexed for loop
                                        if selectedvariable[selectedcounter - place_holder] == 0:#all of the zeros in the chosen rows are eliminated, because they will only make solving more difficult
                                            selectedvariable = numpy.delete(selectedvariable,selectedcounter - place_holder)
                                            useablecorrections = numpy.delete(useablecorrections,selectedcounter - place_holder,axis = 1)
                                            useablemassfrags = numpy.delete(useablemassfrags,selectedcounter - place_holder,axis = 1)
                                            place_holder = place_holder + 1 #saves place when deleting
                                        if usedmolecules[selectedcounter] == 1:# in same loop, the values that are going to be used will be deleted from the arrays getting passed on- that are not solved yet
                                            remaining_correction_values_SLS = numpy.delete(remaining_correction_values_SLS,selectedcounter - place_holder2,axis = 1)
                                            remaining_reference_intensities_SLS = numpy.delete(remaining_reference_intensities_SLS,selectedcounter - place_holder2,axis = 1)
                                            remaining_molecules_SLS = numpy.delete(remaining_molecules_SLS,selectedcounter - place_holder2)
                                            place_holder2 = place_holder2 + 1 #saves place when deleting
                                    useablemassfrags,useablecorrections,useablerawsigs = DistinguishedArrayChooser (useablemassfrags,useablecorrections,useablerawsigs, G.moleculeLikelihoods,G.sensitivityValues)
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
                                        if useablemolecules == []:
                                            useablemolecules = molecules_unedited
                                        #this works in the same way as the brute force method for the remaining matrix at the end of the sls method, where when
                                        #the counter is equal to zero the moleculespecifier chooses all of the array line. This is because when the counter is zero
                                        #the input is always the predicted values according to the inverse method, while afterwards it is the actual signals array
                                        #in this way, for when the counter is zero, the program changes the value of the array to the predicted solutions
                                        if counterforspecifications == 0:
                                            scaledConcentrationsarray = numpy.linalg.solve(useablecorrections,useablerawsigs)
                                        specifications = DataRangeSpecifier(useablemolecules,counterforspecifications,molecules_unedited,conversionfactor,datafromcsv,DataRangeSpecifierlist,scaledConcentrationsarray)
                                        solutions = BruteForce(useablemolecules,specifications,useablecorrections,useablerawsigs,objectiveFunctionType,maxPermutations)
                                    else:
                                        print('The Array Chosen is Singular')
                                        solutions = numpy.zeros(len(remaining_rawsignals_SLS)) # the solutions are made into all zeros if the chosen array is singular
                                    reminder = 0
                                    reminder3 = 1
        selectedvariablesholder.append(selectedvariables)
        selectedcorrectionsholder.append(selectedcorrections)
        selectedmassfragsholder.append(selectedmassfrags)
        selectedrawsigsholder.append(selectedrawsigs)
    if reminder3 == 0:#if no solutions have been found
        solutions = []
        usedmolecules = []
    else: #sets solutions based on the used molecules array that was made
        solutionsholder = numpy.zeros(len(usedmolecules))
        place_holder = 0 
        
        for usedmoleculesIndex in range(len(usedmolecules)):#array-indexed for loop

            if usedmolecules[usedmoleculesIndex] == 1:#when there is a molecule in that position
                solutionsholder[usedmoleculesIndex] = solutions[place_holder]
                place_holder = place_holder + 1 #helps add the two arrays together
        solutions = solutionsholder
    return [molecules_unedited,remaining_reference_intensities_SLS,remaining_correction_values_SLS,remaining_rawsignals_SLS,solutions,remaining_molecules_SLS,usedmolecules]


#this function simply calls the other functions to be used, based on the user input pathway, that means that this
#function can send the sls to unique or common fragments, to inverse or brute method after, and sends the data back 
#and forth between the unique and common fragments for the common fragments method
def SLSMethod(molecules,monitored_reference_intensities,reciprocal_matching_correction_values,rawsignalsarrayline,timeIndex,conversionfactor,datafromcsv,molecules_copy,DataRangeSpecifierlist,SLSChoices,mass_fragment_numbers,permutationNum,scaledConcentrationsarray,objectiveFunctionType,time,maxPermutations=100001):

    # This is creating a local copy of the monitored_reference_intensities which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reference_intensities_SLS = copy.deepcopy(monitored_reference_intensities)

    # This is creating a local copy of the reciprocal_matching_correction_values which will become
    # truncated as the molecules are solved and masses are removed
    remaining_correction_values_SLS = copy.deepcopy(reciprocal_matching_correction_values)

    # This is creating a local copy of the rawsignalsarrayline which will become
    # truncated as the molecules are solved and masses are removed
    remaining_rawsignals_SLS = copy.deepcopy(rawsignalsarrayline)

    # This is creating a local copy of 'molecules' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_molecules_SLS = copy.deepcopy(molecules)
    
    # Read sls options from input files
    [uniqueOrCommon,slsFinish,distinguished] = SLSChoices
                   
    if uniqueOrCommon == 'unique': #user input
        [remaining_molecules_SLS,remaining_reference_intensities_SLS,remaining_correction_values_SLS,remaining_rawsignals_SLS,solutions,molecules_unedited,usedmolecules] = SLSUniqueFragments(remaining_molecules_SLS,remaining_reference_intensities_SLS,remaining_correction_values_SLS,remaining_rawsignals_SLS, timeIndex, time)
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
                if remaining_correction_values_SLS.size != 0: #if the array exists
                    for rowcounter in range(len(remaining_correction_values_SLS[:,0])):#array-indexed for loop
                        zeros = 0
                        for columncounter in range(len(remaining_correction_values_SLS[0,:])):#array-indexed for loop
                            if remaining_correction_values_SLS[rowcounter,columncounter] == 0:#if there is a zero in the array
                                zeros = zeros + 1
                            if columncounter == len(remaining_correction_values_SLS[0,:])-1:#if the whole row has been checked
                                if zeros == len(remaining_correction_values_SLS[0,:])-(desirednum):#if the number of ones if the row is equal to the size being checked for
                                    zerosrows = zerosrows + 1
                    if zerosrows >= desirednum: #if the number of matching rows is equal to or greater than the size being looked for
                        if desirednum == 1: #unique method used if there are rows with only one correction value
                            [remaining_molecules_SLS,remaining_reference_intensities_SLS,remaining_correction_values_SLS,remaining_rawsignals_SLS,solutions,molecules_unedited,usedmolecules] = SLSUniqueFragments (remaining_molecules_SLS,remaining_reference_intensities_SLS,remaining_correction_values_SLS,remaining_rawsignals_SLS, timeIndex, time)
                            if len(usedmolecules) != 0:
                                unique.append([solutions,molecules_unedited,usedmolecules])
                                order.append('unique')
                            place_holder = place_holder +1 #if this happens, than place_holder is made equal to 1, so that the loop must run all the way through again
                            if len(solutions) == 0: #in case no solutions were found
                                place_holder = place_holder - 1 #while loop can stop- so it does not get stuck
                        else: #common fragments method 
                            
                            [molecules_unedited,remaining_reference_intensities_SLS,remaining_correction_values_SLS,remaining_rawsignals_SLS,solutions,remaining_molecules_SLS,usedmolecules] = SLSCommonFragments (remaining_correction_values_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_SLS,remaining_molecules_SLS,scaledConcentrationsarray,molecules_copy,conversionfactor,datafromcsv,DataRangeSpecifierlist,objectiveFunctionType,timeIndex, maxPermutations)
                            if len(usedmolecules) != 0:
                                common.append([solutions,molecules_unedited,usedmolecules])
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
                    usedmolecules = unique[0][2] #usedmolecules is set based on the first one used as well
                if order[ordercounter] == 'common':#if the last item of the order list was common
                    solutions = common[0][0]
                    molecules_unedited = common[0][1]
                    usedmolecules = common[0][2]
    else: # throw error because uniqueOrCommon is not set properly
        raise ValueError("The value of 'uniqueOrCommon' is {}, it should be either 'common' or 'unique'".format(uniqueOrCommon)) 
        
    #if the sls method does not solve for all the molecules, then the rest are sent to the inverse method 
    #where the remaining matrix is solved
    if remaining_correction_values_SLS.size != 0:#if everything hasn't already been solved
        if slsFinish == 'inverse':#if inverse finish is chosen
             if distinguished == 'yes':#user input, choosing between distinguished inverse method or combinations method
                 concentrationsFromFinisher = InverseMethodDistinguished (remaining_reference_intensities_SLS,remaining_correction_values_SLS,remaining_rawsignals_SLS)
             else:
                 concentrationsFromFinisher = InverseMethod (remaining_correction_values_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_SLS,mass_fragment_numbers,remaining_molecules_SLS,'composition')
           
        if slsFinish == 'brute':#if brute method is chosen
            if timeIndex == 0:#the first time is always run through the inverse method, where the ranges can use this information the loops afterwards
                concentrationsFromFinisher = InverseMethod (remaining_correction_values_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_SLS,mass_fragment_numbers,remaining_molecules_SLS,'composition')
                
                if math.log(permutationNum,5) >= len(remaining_molecules_SLS):
                    #Ashi believes the above comparision is to ensure each molecule has at least 5 concentrations checked
                    specifications = DataRangeSpecifier(remaining_molecules_SLS,timeIndex,molecules_copy,conversionfactor,datafromcsv,DataRangeSpecifierlist,concentrationsFromFinisher)
                    concentrationsFromFinisher = BruteForce(remaining_molecules_SLS,specifications,remaining_correction_values_SLS,remaining_rawsignals_SLS,objectiveFunctionType,maxPermutations)
                else:
                    print("Warning: The number of permutations requested is too small to allow for 5 possibilities per molecule."
                          + "Switching to use Inverse instead of Brute for slsFinish for this datapont." 
                          + "If you wish to use Brute, increase the size of permutations in the user input file.")
                    concentrationsFromFinisher = InverseMethod (remaining_correction_values_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_SLS,mass_fragment_numbers,remaining_molecules_SLS,'composition')
            else: #after the first time these functions are called
                if math.log(permutationNum,5) >= len(remaining_molecules_SLS):
                    #Ashi believes the above comparision is to ensure each molecule has at least 5 concentrations checked
                    specifications = DataRangeSpecifier(remaining_molecules_SLS,timeIndex,molecules_copy,conversionfactor,datafromcsv,DataRangeSpecifierlist,
                                                        scaledConcentrationsarray)
                    concentrationsFromFinisher = BruteForce(remaining_molecules_SLS,specifications,remaining_correction_values_SLS,remaining_rawsignals_SLS,objectiveFunctionType, maxPermutations)
                else:
                    concentrationsFromFinisher = InverseMethod (remaining_correction_values_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_SLS,mass_fragment_numbers,remaining_molecules_SLS,'composition')
        # the concentrations that were solved for by the Finisher are stored as a list
        # to make them easier to use and then discard in the for loop as they are added to solutions
        remainingMolecules = list(concentrationsFromFinisher.copy())
        #print remainingMolecules
        #adds the finisher solutions with the inital analysis solutions 
        for moleculecounter in range(len(molecules_unedited)):
            if len(usedmolecules) == 0:
                print("Warning: if you have chosen to use unique fragment SLS and your data has no unique fragments (not unique to any molecule), "\
                "then the program may be about to crash. If it does, change SLS method to common.")
            #if the molecule wasn't solved for in the inital analysis
            if usedmolecules[moleculecounter] == 0:
                # then add the appropriate Finisher concentration for that molecule  
                solutions[moleculecounter] = remainingMolecules.pop(0) 
            
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
                              maxPermutations = 100001):

    # This is creating a local copy of the monitored_reference_intensities which will become
    # truncated as the molecules are solved and masses are removed
    remaining_reference_intensities_filter = copy.deepcopy(monitored_reference_intensities)

    # This is creating a local copy of the reciprocal_matching_correction_values which will become
    # truncated as the molecules are solved and masses are removed
    remaining_correction_values_SLS = copy.deepcopy(reciprocal_matching_correction_values)

    # This is creating a local copy of the rawsignalsarrayline which will become
    # truncated as the molecules are solved and masses are removed
    remaining_rawsignals_SLS = copy.deepcopy(rawsignalsarrayline)

    # This is creating a local copy of 'molecules' which will become
    # truncated as the molecules are solved and masses are removed
    remaining_molecules_SLS = copy.deepcopy(molecules)
    
    absentmolecules1 = numpy.ones([1,len(remaining_correction_values_SLS[0,:])])
    absentmolecules = absentmolecules1[0]
    molecules_unedited = copy.deepcopy(molecules)
    place_holder = 0
    place_holder2 = 0
    monitored_reference_intensities_copy = monitored_reference_intensities
    rawsignalsarray_copy = remaining_rawsignals_SLS
    summation = sum(rawsignalsarray_copy)
    [applyRawSignalThresholds,rawSignalThresholdValue,sensitivityThresholdValue,rawSignalThresholdDivider,rawSignalThresholdLimit,rawSignalThresholdLimitPercent] = ThresholdList
    #this section of the code enables the function to eliminate from the raw signal array (for the purpose of this function only)
    #the highest value in the array, if it makes up over 90 percent of the raw signals present. This is useful because if one of 
    #the signals does become this great then it will eliminate all the other signals present when it becomes very high
    if rawSignalThresholdLimitPercent == []:#if no user input
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
        sensitivityThresholdValue = 5
    elif len(sensitivityThresholdValue) > 1:#user input
        sensitivityThresholdValue = sensitivityThresholdValue[0]
    #This for loop goes through all of the rows of the remaining_rawsignals_SLS and finds the values that are lower than the given 
    #threshold as defined above, and these rows are deleted from the remaining_rawsignals_SLS, the remaining_correction_values_SLS and the matching
    #mass fragments array, all using a place hold so that even when the rows are deleted, the counting will still work for 
    #the loop. The a nested for loop looks in the matching mass fragments array (a copy, that is not changed in the previous
    #loop), and finds which molecules from each of those deleted rows have a relative intensity above the given threshold 
    #value in their respective mass fragmentation pattern. These molecules are saved as zeros in the absentmolecules array
    for rawcounter in range(len(remaining_rawsignals_SLS)):#array-indexed for loop
        if remaining_rawsignals_SLS[rawcounter-place_holder] < rawSignalThresholdValue:#any values below threshold are deleted
            remaining_rawsignals_SLS = numpy.delete(remaining_rawsignals_SLS,rawcounter-place_holder,axis = 0)
            remaining_correction_values_SLS = numpy.delete(remaining_correction_values_SLS,rawcounter-place_holder,axis = 0)
            remaining_reference_intensities_filter = numpy.delete(remaining_reference_intensities_filter,rawcounter-place_holder,axis = 0)
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
            remaining_correction_values_SLS = numpy.delete(remaining_correction_values_SLS,absentmoleculescounter-place_holder2,axis = 1)
            place_holder2 = place_holder2 + 1    
    #this last section of code just adds the absent molecules together with the solutions to get all of the molecules needed, and
    #prints off the results using a copy of the molecules array made earlier; this part of the function is very similar to what the
    #sls function has to add its array to that of what it gets back from the inverse function that it uses.
    if any(absentmolecules) == 0:#everything is below threshold-not probable
        print('No Significant Peaks')
    else:
        if answer == 'inverse':#if the inverse method is wanted
            if distinguished == 'yes':#distinguished method
                solutions = InverseMethodDistinguished(remaining_reference_intensities_filter,remaining_correction_values_SLS,remaining_rawsignals_SLS)
            else:#combinations method
                solutions = InverseMethod(remaining_correction_values_SLS,remaining_rawsignals_SLS,remaining_reference_intensities_filter,mass_fragment_numbers,remaining_molecules_SLS,'composition')
        if answer == 'sls':#sls method
            solutions = SLSMethod(remaining_molecules_SLS,remaining_reference_intensities_filter,remaining_correction_values_SLS,remaining_rawsignals_SLS,timeIndex,conversionfactor,datafromcsv,molecules_unedited,DataRangeSpecifierlist,SLSChoices,mass_fragment_numbers,permutationNum,scaledConcentrationsarray,objectiveFunctionType, time, maxPermutations)
        timeIndex = 0
        for moleculecounter in range(len(molecules_unedited)):#array-indexed for loop
            if absentmolecules[moleculecounter] == 1:#gets index for present molecule
                absentmolecules[moleculecounter] = solutions[moleculecounter-timeIndex]
            else:#keeps index in check
                timeIndex = timeIndex+1
    solutions = absentmolecules
    return solutions
    
#This little bit of code enables the user to input a couple values so that the conversion between signal relative to CO and 
#concentration can be found, and as such the rest of the concentrations can be found as well. If no data is input here, then
#there are no concentrations printed out, only signals. (if they are both printed, it is done on separate excel sheets)
def RatioFinder (ReferenceData, ExperimentData, concentrationFinder,molecule,moleculeConcentration,massNumber,moleculeSignal,units):
    #By putting these declarations first, the variables will always be initilized
    #preventing errors further on in the code
    ExperimentData.conversionfactor = 0
    ExperimentData.units = 'torr'
    
    if concentrationFinder == 'yes':#user input
        ExperimentData.units = units
        for moleculecounter in range(len(ReferenceData.molecules)):#array-indexed for loop
            for masscounter in range(len(ExperimentData.mass_fragment_numbers)):#array-indexed for loop
                if molecule == ReferenceData.molecules[moleculecounter]:#gets molecule index
                    if massNumber == ExperimentData.mass_fragment_numbers[masscounter]:#gets index
                        ExperimentData.conversionfactor = (moleculeConcentration*ReferenceData.reciprocal_matching_correction_values[masscounter,moleculecounter])/float(moleculeSignal)
    return ExperimentData
    
    
#this function is going to be rather simple, but it will be the forward function, that simulates raw signals from the calculated
#signals that we acquired, which will later be printed out in the main() function. The inputs for this function are the signals 
#array as well as the correction values, in order to simulate the raw signals
def RawSignalsSimulation (scaledConcentrationsarray,reciprocal_matching_correction_values):
    simulateddata = numpy.zeros([len(scaledConcentrationsarray[:,0]),len(reciprocal_matching_correction_values[:,0])+1])#a simulated data array made the height of the signals array and the width equal to the correction arrays height
    times = scaledConcentrationsarray[:,0]#the first row of the signals array is the times
    scaledConcentrationsarray = scaledConcentrationsarray[:,1:] #because one of the lines here is the times, and it need to be removed
    for scaledtimeIndex in range(len(scaledConcentrationsarray[:,0])):#array-indexed for loop
        simulateddata[scaledtimeIndex:scaledtimeIndex+1,1:] = numpy.transpose(numpy.matrix(reciprocal_matching_correction_values) * numpy.matrix(numpy.vstack(scaledConcentrationsarray[scaledtimeIndex,:])))#the data is simulated by multiplying the matrix of correction values by the raw signals for each row
    simulateddata[:,0] = times #the times are added back in so they can be printed more easily
    return simulateddata
    
    
#this function is a necessity if there are negatives in your answer, it finds those negatives and sends those negatives along
#with the molecule that affects them the most and sends them both to the brute method so that they can both be solved again
# the molecule with a larger amount is checked for signals near its original signal, while the other molecule is checked for
#data from zero, up to the bigger molecule's signal
def NegativeAnalyzer (solutionsline,reciprocal_matching_correction_values,rawsignalsarrayline,molecules,objectiveFunctionType,maxPermutations=100001):
    solutionslinedata = solutionsline[1:]# gets rid of the times for our data array
    negatives = []
    indexes = []
    for solutionsIndex in range(len(solutionslinedata)): #looks through the line
        if solutionslinedata[solutionsIndex] < 0: #if there is a value below zero it keeps the value and the index
            negatives.append(solutionslinedata[solutionsIndex])
            indexes.append(solutionsIndex)
    NGstart = timeit.default_timer()
    if len(negatives) > 0:#if there are negatives then the function runs
        for negativesIndex in range(len(negatives)):#does this for each negative
            for matchCorrIndexCol in range(len(reciprocal_matching_correction_values[:,0])):#looks through the correction values
                if reciprocal_matching_correction_values[matchCorrIndexCol,indexes[negativesIndex]] == max(reciprocal_matching_correction_values[:,indexes[negativesIndex]]):#finds the index of the negative molecule's largest correction value 
                    correction1index = matchCorrIndexCol
            presentmoleculeslist = []
            for matchCorrIndexRow in range(len(reciprocal_matching_correction_values[0,:])):#goes through the correction values
                if reciprocal_matching_correction_values[correction1index,matchCorrIndexRow] != 0:#if the molecule has a relative intensity (other than zero) at the mass fragment chosen (by the last loop)
                    presentmoleculeslist.append(1)
                else: #if there is no molecule a zero is appended to the list
                    presentmoleculeslist.append(0)
            presentmoleculesarray = numpy.array(presentmoleculeslist)
            solutionslinepresentarray = solutionslinedata*presentmoleculesarray #the ones and zeros list is multiplied by the solutions, allowing only molecules with the mass fragment desired to be selected later
            for solutionsIndex2 in range(len(solutionslinedata)):#goes through the solution line
                if max(solutionslinepresentarray) != 0:#if there are any above zero
                    if solutionslinepresentarray[solutionsIndex2] == max(solutionslinepresentarray):#the highest value is used
                        correction2index = solutionsIndex2
                else:
                    if solutionslinedata[solutionsIndex2] == max(solutionslinedata):#if there are none with that mass fragment, the highest solution is chosen
                        correction2index = solutionsIndex2
            arrayamalgam = reciprocal_matching_correction_values[:,indexes[negativesIndex]],reciprocal_matching_correction_values[:,correction2index]#an array amalgam is made with  two columns for the two chosen molecules
            arrayamalgam = numpy.array(arrayamalgam)
            arrayamalgam = numpy.transpose(arrayamalgam) #the array is transposed so it can be used in matrix multiplication
            solutionslinedata[indexes[negativesIndex]] = 0#the index of the molecule chosen is made it a zero
            maximum = solutionslinedata[correction2index] #the second molecule chosen
            solutionslinedata[correction2index] = 0#the second value is made zero too
            reciprocal_matching_correction_values_copy = numpy.array(reciprocal_matching_correction_values)
            reciprocal_matching_correction_values_copy[:,indexes[negativesIndex]] = 0#the two columns in the correction values array are made into zeros
            reciprocal_matching_correction_values_copy[:,correction2index] = 0
            rawsignalsubtractionvalue = numpy.matrix(reciprocal_matching_correction_values_copy)*numpy.matrix(numpy.vstack(solutionslinedata))#The raw signals are simulated from the correction values and raw signals containing all molecules except those not chosen
            rawsignalsarraylinecopy = rawsignalsarrayline - numpy.array(rawsignalsubtractionvalue)#The simulated raw signals are subtracted from the actuals and the left over raw signals are due to only the molecules left
            ranges = numpy.linspace(0,maximum/float(10),50)#the negative molecule is checked for between the higher molecule's signal/10, in 100ths of the range
            userange = ranges[1]-ranges[0] #the increments are calculated here
            specifications = [(0,maximum/float(10),userange),(maximum/float(2),maximum*2,maximum*0.15)]#the specifications array is made here, with the higher molecule being checked in ten places within a factor of 2 of itself
            if sum(specifications[0]) == 0 and sum(specifications[1] == 0):
                return solutionsline # This means that all of the values were negative or zero
            answers = BruteForce(molecules,specifications,arrayamalgam,rawsignalsarraylinecopy,objectiveFunctionType,maxPermutations)#brute method used- 200 permutations- 20*10 from the increments above
            solutionslinedata[indexes[negativesIndex]] = answers[0]#sets the first solution
            solutionslinedata[correction2index] = answers[1]#sets the second
    solutionsline[1:] = solutionslinedata #sets the new data, with the times
    print((timeit.default_timer() - NGstart))
    return solutionsline

'''This function exports all XYYY Data to the User specified document (usually a CSV)'''
#TODO Future Development: This function could benefit from creating folders to store the different
#runs of the program if it is run in the same directory as before
## NOTE: This function leaves a tailing comma on each line of the
## created csv file. When pandas is used to read the csv this
## will result in the creation of a column of 'nan'
## ImportWorkingData() has been modified to remove this column
def ExportXYYYData(outputFileName, data, colIndex, abscissaHeader = 'Mass', fileSuffix = '', dataType = None, rowIndex = [], units = None): 
    formatedColIndex = colIndex
    if dataType == 'preProcessed' or dataType == 'simulated' or dataType == 'Experiment':
        formatedColIndex = ['m%s' % MFNumber for MFNumber in colIndex]
    if dataType == 'scaled':
        formatedColIndex = ['%s Concentration Relative to CO' % molecule for molecule in colIndex]
    if dataType == 'concentration':
        label = ' Concentration(%s)' % units
        formatedColIndex = [molecule + label for molecule in colIndex]
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
        fullArrayToExport = numpy.vstack((formatedColIndex,data))
    #occasionally, abscissaHeader needs to be inserted without rowIndex being used
    except ValueError: 
        formatedColIndex = numpy.hstack((abscissaHeader,formatedColIndex))
        fullArrayToExport = numpy.vstack((formatedColIndex,data))
        
    #if the row index isn't included in the data, then add it 
    if rowIndex != []:    
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
    numpy.savetxt(filename, fullArrayToExport, delimiter = ',', fmt ="%s")
  

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
'''
def Draw(times, data, molecules, concentrationFinder, units, graphFileName = '', fileSuffix = ''):
    
    import matplotlib.pyplot as plt
    
    colormap = plt.cm.gist_ncar
    plt.gca().set_color_cycle([colormap(i) for i in numpy.linspace(0, 0.9, len(data[0,:]))])
    ax = plt.subplot(111)
    for x in range(len(data[0,:])):
        ax.plot(times,data[:,x])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(molecules,loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('Time (minutes)')
    if concentrationFinder == 'yes':
        plt.ylabel('Concentration (%s)'%(units))
    else:
        plt.ylabel('Concentration Relative to CO')

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

    plt.show()

        

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
def createSLSUniqueOrderFile(abscissaHeader, molecules):
    outputFileName = 'SLSUniqueOrder.csv'
    if G.iterativeAnalysis:
        #then the filename will have a suffix attached
        outputFileName = outputFileName[:-4] + '_iter_%s' %G.iterationNumber + outputFileName[-4:]
    
    with open(outputFileName,'w') as fp:
        # Headers
        fp.write('{},{}'.format(
            'Data Point',abscissaHeader))
        # in writing molecule names in the header here
        # I assume the actual numbers are written in the same
        # order later 
        for molecule in molecules:
            fp.write(",{}".format(molecule))
        fp.write('\n')

'''
This function takes in the end result of everything and exports it to the 
log file. 
'''
def PopulateLogFile():
    filename6 = 'LogFile.txt' #the log file is printed here
    f6 = open(filename6,'a')
    f6.write('\n')
    f6.write('referenceFileName = %s \n'%(G.referenceFileName))
    f6.write('form  = %s \n'%(G.form))
    f6.write('dataToAnalyzeFileName = %s \n'%(G.dataToAnalyzeFileName ))
    if G.timeRangeLimit == 'yes':#some of the lines in the backgroundinput file don't need to be printed unless a selection is made, so the if statements here make that happen
        f6.write('timeRangeLimit = %s \n'%(G.timeRangeLimit))
        f6.write('timeRangeStart = %s \n'%(G.timeRangeStart))
        f6.write('timeRangeFinish = %s \n'%(G.timeRangeFinish))
    if G.backgroundMassFragment != []:#some of the lines in the backgroundinput file don't need to be printed unless a selection is made, so the if statements here make that happen
        f6.write('backgroundMassFragment = %s \n'%(G.backgroundMassFragment))
        f6.write('backgroundSlopes = %s \n'%(G.backgroundSlopes))
        f6.write('backgroundIntercepts = %s \n'%(G.backgroundIntercepts))
    if G.tuningCorrection == 'yes':
        f6.write('tuningCorrection = %s \n'%G.tuningCorrection)
        f6.write('referenceCorrectionCoefficientA = %s \n'%(G.referenceCorrectionCoefficients['A']))
        f6.write('referenceCorrectionCoefficientB = %s \n'%(G.referenceCorrectionCoefficients['B']))
        f6.write('referenceCorrectionCoefficientC = %s \n'%(G.referenceCorrectionCoefficients['C']))
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
        f6.write('moleculesRange = %s \n'%(G.moleculesRange))
        if G.csvFile == 'yes':
            f6.write('csvFileName = %s \n'%(G.csvFileName))
        else:
            f6.write('lowerBound = %s \n'%(G.dataLowerBound))
            f6.write('higherBound = %s \n'%(G.dataUpperBound))
            f6.write('increments = %s \n'%(G.increments))
        f6.write('permutationNum = %s \n'%(G.permutationNum))
    if G.extractReferencePatternFromDataOption == 'yes':
        f6.write('extractReferencePatternFromDataOption = %s \n'%(G.extractReferencePatternFromDataOption))
        f6.write('rpcMoleculesToChange = %s \n'%(G.rpcMoleculesToChange))
        f6.write('rpcMoleculesToChangeMF = %s \n'%(G.rpcMoleculesToChangeMF))
        f6.write('rpcTimeRanges = %s \n'%(G.rpcTimeRanges))
    if G.applyReferenceMassFragmentsThresholds == 'yes':
        f6.write('applyReferenceMassFragmentsThresholds = %s \n'%(G.applyReferenceMassFragmentsThresholds))
        f6.write('referenceMassFragmentFilterThreshold = %s \n'%(G.referenceMassFragmentFilterThreshold))
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
    if G.applyRawSignalThresholds == 'yes':
        f6.write('applyRawSignalThresholds = %s \n'%(G.applyRawSignalThresholds))
        f6.write('rawSignalThresholdValue = %s \n'%(G.rawSignalThresholdValue))
        f6.write('sensitivityThresholdValue = %s \n'%(G.sensitivityThresholdValue))
        f6.write('rawSignalThresholdDivider = %s \n'%(G.rawSignalThresholdDivider))
    if G.rawSignalThresholdLimit  == 'yes':
        f6.write('rawSignalThresholdLimit = %s \n'%(G.rawSignalThresholdLimit))
        f6.write('rawSignalThresholdLimitPercent  = %s \n'%(G.rawSignalThresholdLimitPercent))
    if G.negativeAnalyzerYorN == 'yes':
        f6.write('negativeAnalyzerYorN = %s \n'%(G.negativeAnalyzerYorN))
    if G.dataAnalysis == 'yes':
        f6.write('answer = %s \n'%(G.solverChoice))
        if G.solverChoice == 'sls':
            f6.write('uniqueOrCommon = %s \n'%(G.uniqueOrCommon))
            f6.write('slsFinish = %s \n'%(G.slsFinish))
            if G.slsFinish == 'brute':
                f6.write('objectiveFunctionType = %s \n'%(G.objectiveFunctionType))
            if G.slsFinish == 'inverse':
                f6.write('distinguished = %s \n'%(G.distinguished))
        if G.solverChoice == 'inverse':
            f6.write('distinguished = %s \n'%(G.distinguished))
    if G.concentrationFinder == 'yes':
        f6.write('concentrationFinder = %s \n'%(G.concentrationFinder))
        f6.write('molecule = %s \n'%(G.molecule))
        f6.write('moleculeSignal = %s \n'%(G.moleculeSignal))
        f6.write('massNumber = %s \n'%(G.massNumber))
        f6.write('moleculeConcentration = %s \n'%(G.moleculeConcentration))
        f6.write('units = %s \n'%(G.units))
    f6.write('resolvedScaledConcentrationsOutputName  = %s \n'%(G.resolvedScaledConcentrationsOutputName ))
    f6.write('concentrationsOutputName = %s \n'%(G.concentrationsOutputName))
    f6.write('simulatedSignalsOutputName = %s \n'%(G.simulatedSignalsOutputName))
    f6.write('Run Time %.5f seconds \n'%(G.checkpoint-G.start))
    f6.write('###################################################################### \n')
    f6.close()#once the log file is printed the program is finished
    return None

##################################################################################################################
###############################################Algorithm Part 3: Main Control Function ###################################
##################################################################################################################
def main():
     #This section is to overwrite the UI if iterative analysis is in the process of being run. 
    highestIteration = int(FindHighestDirNumber("_iter_"))
    iterationDirectorySuffix = '_iter_%s' %str(highestIteration)
    for directoryName in os.listdir():
        if iterationDirectorySuffix in directoryName:
            userInputName = 'UserInput%s' %iterationDirectorySuffix
            userInputPath = '%s.%s' %(directoryName, userInputName)
            global G
            UserInput2 = importlib.import_module('..%s' %userInputName, '%s' %userInputPath)
            G = UserInput2
            break
    if G.iterativeAnalysis:    
        G.iterationNumber = highestIteration
        G.iterationSuffix = iterationDirectorySuffix
    elif not G.iterativeAnalysis:
        G.iterationSuffix = ''
    
    #Record the time
    G.start = timeit.default_timer()
    G.checkpoint = timeit.default_timer()
    
    #if this is not the first iterative run, then the required files are all stored in the highest iteration directory
    if G.iterativeAnalysis and G.iterationNumber != 1:
        #implied arguments for this function are G.referenceFileName and G.dataToAnalyzeFileName
        IterationDirectoryPreparation(G.iterativeAnalysis, G.iterationNumber)
    
    #initalize the data classes with the data from given Excel files
    #These are being made into globals primarily for unit testing and that functions are expected to receive the data as arguments rather than accessing them as globals
    global ReferenceDataList
    global ExperimentData
    global currentReferenceData
    
    [exp_mass_fragment_numbers, exp_abscissaHeader, exp_times, exp_rawCollectedData, exp_dataToAnalyzeFileName]=readDataFile(G.dataToAnalyzeFileName)
    ExperimentData = MSData(exp_mass_fragment_numbers, exp_abscissaHeader, exp_times, exp_rawCollectedData, dataToAnalyzeFileName=exp_dataToAnalyzeFileName)
    ReferenceDataList = GenerateReferenceDataList(G.referenceFileName,G.form)
    currentReferenceData = ReferenceDataList[0]
    #Prints a warning if the user has more reference files than specified time ranges
    if len(G.referenceFileName) > len(G.referencePatternTimeRanges):
        print("WARNING: There are more reference files given than time ranges")

    #if this is the first iterative run, then the reference and experimental files need to have been imported before the iteration can begin
    if G.iterativeAnalysis and G.iterationNumber == 1 :
        #implied arguments for the following function are G.referenceFileName and G.dataToAnalyzeFileName
        IterationFirstDirectoryPreparation(G.iterativeAnalysis, G.iterationNumber)

    # Skip preProcessing all together if we are loading analyzed data
    if(G.dataAnalysis == 'load'):
        print("DataAnalysis set to 'load': skipping preprocessing")

        # Skip preprocessing
        G.preProcessing = 'skip'
        
        
    if(G.preProcessing == 'yes'):
        
        if G.iterativeAnalysis:
            #create a copy of the Reference Data
             ReferenceDataFullCopy = copy.deepcopy(currentReferenceData)
             
        # Trim the reference data according to the selected molecules list
        if G.specificMolecules == 'yes' or G.iterativeAnalysis:
            (currentReferenceData.provided_reference_intensities, currentReferenceData.electronnumbers, currentReferenceData.molecules, currentReferenceData.mass_fragment_numbers_monitored) = TrimDataMolecules(currentReferenceData, G.chosenMolecules) 
             
        # Perform the actual data preprocessing on ExperimentData
        ExperimentData = DataInputPreProcessing(ExperimentData)
        print("Data PreProcessing completed")
        
        if G.iterativeAnalysis:
            #make a copy of the experimental data
            ExperimentDataFullCopy = copy.deepcopy(ExperimentData)
            
        #Trim the experimental data according to the mass fragments in G.chosenMassFragments and the reference data
        (ExperimentData.workingData, ExperimentData.mass_fragment_numbers) = trimDataMasses(ExperimentData, currentReferenceData) 
        
        #This graph call is graphing fully preprocessed data.
        if G.grapher == 'yes':
            print('PreProcessed Data Graph')
            Draw(ExperimentData.times, ExperimentData.workingData, ExperimentData.mass_fragment_numbers, 'no', 'Amp', graphFileName = 'PreprocessingAfterSmoothing', fileSuffix = G.iterationSuffix)

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
    
        #The iterative analysis preprocessing creates the proper export folder and exports the unused reference data
        if G.iterativeAnalysis:
            IADirandVarPopulation(G.iterativeAnalysis, G.chosenMassFragments, G.chosenMolecules, ExperimentData, currentReferenceData, ReferenceDataFullCopy)
    
    elif(G.preProcessing == 'skip'):

        # Even if we skip the real preProcessing some manipulations must be performed on the
        # MS data

        # Trim the data according to the mass fragments in G.chosenMassFragments and the reference data
        (ExperimentData.workingData, ExperimentData.mass_fragment_numbers) = trimDataMasses(ExperimentData, currentReferenceData)
        
        # Output to make sure user knows we are skipping Preprocessing
        G.timeSinceLastCheckPoint = timeit.default_timer() - G.checkpoint
        G.checkpoint = timeit.default_timer()
        print("Preprocessing set to 'skip'. Performed only mandatory PreProcessing with Time: {:.4f}".format(G.timeSinceLastCheckPoint))
            
    elif(G.preProcessing == 'load' and G.dataAnalysis == 'yes'):
            
        #This function call loads the preprocessed data
        print("Loading the preprocessed data from file '{}'".format(G.preProcessedDataOutputName))
        ExperimentData.workingData, ExperimentData.mass_fragment_numbers, ExperimentData.times = ImportWorkingData(G.preProcessedDataOutputName)

        # Trim the data according to the mass fragments in G.chosenMassFragments and the reference data
        (ExperimentData.workingData, ExperimentData.mass_fragment_numbers) = trimDataMasses(ExperimentData, currentReferenceData)
             
        # Output to make sure user knows we are loading Preprocessing
        G.timeSinceLastCheckPoint = timeit.default_timer() - G.checkpoint
        G.checkpoint = timeit.default_timer()
        print("Preprocessed data loaded. Time: {:.4f}".format(G.timeSinceLastCheckPoint))

    else:
        # if we are here then 'G.preProcessing' != ('yes' or 'skip' or 'load')
        raise ValueError("The value of preProcessing is not set appropriately, it should be 'yes', 'skip' or 'load'." +
                         "Or you are attempting to load pre-processed data without running data analysis")

    #TODO make a variable allMoleculesAnalyzed that is a list containing all the molecules analyzed so far
    ## Here perform the ReferenceData preprocessing that is required regardless of the selection for 'G.preProcessing'
    # and needed if G.dataAnalysis == 'load' or 'yes'
    if (G.dataAnalysis == 'yes' or G.dataAnalysis =='load'):
        #for loop to preprocess all MSReference objects
        for i in range(len(ReferenceDataList)):
            # Reference Pattern Changer
            if G.extractReferencePatternFromDataOption == 'yes':
                ReferenceDataList[i].provided_reference_intensities = ExtractReferencePatternFromData(ExperimentData, ReferenceDataList[i], G.rpcMoleculesToChange, G.rpcMoleculesToChangeMF, G.rpcTimeRanges)
                ReferenceDataList[i].ExportCollector('ExtractReferencePatternFromData',use_provided_reference_intensities = True)
                print('ReferencePatternChanger complete')
                    
            # Some initial preprocessing on the reference data
            ReferenceDataList[i] = ReferenceInputPreProcessing(ReferenceDataList[i])
            # Set the ReferenceData.monitored_reference_intensities and
            # ReferenceData.reciprocal_matching_correction_values fields
            # based on the masses in ExperimentData.mass_fragment_numbers
            ReferenceDataList[i] = Populate_reciprocal_matching_correction_values(ExperimentData.mass_fragment_numbers,ReferenceDataList[i])
            # Remove reference species that have no mass fragment data
            # from the ReferenceData fields monitored_reference_intensities, reciprocal_matching_correction_values and molecules
            ## TODO: Consider changing this function to take the array directly i.e.
            ## (monitored_reference_intensities) so that it can potentially be applied to other arrays
            ## like ReferenceData.standardized_reference_intensities
            ReferenceDataList[i] = UnnecessaryMoleculesDeleter(ReferenceDataList[i])
            ReferenceDataList[i].ExportCollector('UnnecessaryMoleculesDeleter')
    
            # Export the reference data files that have been stored by ReferenceData.ExportCollector
            ReferenceDataList[i].ExportFragmentationPatterns()

            
    if (G.dataAnalysis == 'yes'):
                
        # Reset the checkpoint timer for the data analysis section
        G.checkpoint = timeit.default_timer()
        #check to make sure that there are enough mass fragments to solve for each variable. 
        if len(ExperimentData.mass_fragment_numbers) < len(currentReferenceData.molecules):
            raise SystemError("There are too few mass fragments to solve for the number of molecules provided. \nData Analysis has been ended.")
            
        ##Start: Preparing data for data analysis based on user input choices
        # If we are only interested in a subset of the MS data
        # remove the irrelevant mass data series from ExperimentData.mass_fragment_numbers
        # and the corresponding colums from ExperimentData.workingData
        if G.specificMassFragments == 'yes':
            print("MassFragChooser")
            (ExperimentData.workingData, ExperimentData.mass_fragment_numbers) = DataFunctions.KeepOnlySelectedYYYYColumns(ExperimentData.workingData,
                                                                                                                           ExperimentData.mass_fragment_numbers,
                                                                                                                           G.chosenMassFragments)
            ExperimentData.ExportCollector("MassFragChooser")

        
        # Since SlSUniqueFragments is potentially used in a number of the analysis options
        # set up the 'SLSUniqueOrder.csv' file headers here
        # note that this open overwrites previous file contents
        # while subsequent opens to this file will append
        if G.SLSUniquePrint == 'yes':
            createSLSUniqueOrderFile(ExperimentData.abscissaHeader,
                                     currentReferenceData.molecules)
            
        #this numpy.zeros line is going to be the array that holds all of the answers before they are printed out, which
        #is done in order to save time and decrease expense
        concentrationsScaledToCOarray = numpy.zeros(len(currentReferenceData.molecules)+1)
        concentrationsarray = numpy.zeros(len(currentReferenceData.molecules)+1)

        # Loading user choices for data analysis
        DataRangeSpecifierlist = [G.dataRangeSpecifierYorN, G.signalOrConcentrationRange,
                                  G.csvFile, G.moleculesRange, G.csvFileName,G.dataUpperBound,
                                  G.dataLowerBound, G.increments, G.permutationNum]
        SLSChoices = [G.uniqueOrCommon, G.slsFinish, G.distinguished]
        ThresholdList = [G.applyRawSignalThresholds, G.rawSignalThresholdValue, G.sensitivityThresholdValue,
                         G.rawSignalThresholdDivider, G.rawSignalThresholdLimit, G.rawSignalThresholdLimitPercent]
        
        # Calculate a coefficient for doing a unit conversion on concentrations
        ExperimentData = RatioFinder(currentReferenceData, ExperimentData, G.concentrationFinder,
                                      G.molecule, G.moleculeConcentration, G.massNumber, G.moleculeSignal, G.units)
	##End: Preparing data for data analysis based on user input choices
        
        for timeIndex in range(len(ExperimentData.workingData[:,0])):#the loop that runs the program to get a set of signals/concentrations for each time   
            #This print statement was used to track the progress of the program during long analysis runs
            if ((timeIndex % 100) == 0 and timeIndex != 0):
                print(timeIndex)
           

            ## TODO: Find out why RawSignalsArrayMaker() takes longer to run when preprocessed data is
            # computed directly rather than loaded. It doesn't seem to effect rawsignalsarrayline in
            # either case so not a priority. 
            rawsignalsarrayline = RawSignalsArrayMaker(currentReferenceData.mass_fragment_numbers_monitored,
                                                       ExperimentData.mass_fragment_numbers,ExperimentData.workingData,
                                                       timeIndex,currentReferenceData.referenceabscissa)#gets the collected values that will be solved
            
            
            if G.applyRawSignalThresholds == 'yes':#user input, this function calls either sls or inverse, deletes thresholds
                    solutions =RawSignalThresholdFilter(G.distinguished, currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline,
                                                         currentReferenceData.monitored_reference_intensities,currentReferenceData.molecules,timeIndex,ExperimentData.mass_fragment_numbers,
                                                         ThresholdList,G.solverChoice,ExperimentData.times[timeIndex],ExperimentData.conversionfactor,ExperimentData.datafromcsv,
                                                         DataRangeSpecifierlist,SLSChoices,G.permutationNum,concentrationsScaledToCOarray,G.objectiveFunctionType, G.maxPermutations)
            else:#otherwise the main analysis functions are called
                
                if G.solverChoice == 'inverse':#user input, the inverse method
                    if G.distinguished == 'yes':#user input, choosing between distinguished inverse method or combinations method

                        solutions = InverseMethodDistinguished(currentReferenceData.monitored_reference_intensities,currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline)
                    else:
                        solutions = InverseMethod(currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline,currentReferenceData.monitored_reference_intensities,ExperimentData.mass_fragment_numbers,currentReferenceData.molecules,'composition')
    
                elif G.solverChoice == 'sls':#user input, the SLS method is chosen)
                    solutions = SLSMethod(currentReferenceData.molecules,currentReferenceData.monitored_reference_intensities,currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline, timeIndex, ExperimentData.conversionfactor, ExperimentData.datafromcsv,currentReferenceData.molecules,DataRangeSpecifierlist,SLSChoices,ExperimentData.mass_fragment_numbers,G.permutationNum,concentrationsScaledToCOarray,G.objectiveFunctionType,ExperimentData.times[timeIndex],G.maxPermutations)
            
            arrayline = []

            for moleculecounter in range(len(currentReferenceData.molecules)):#array-indexed for loop, this is the same data structure as the inverse method above once the line of solutions is found, see above for comments
                if moleculecounter == 0:#only for the first loop will times be added to new collected data
                    arrayline.append(ExperimentData.times[timeIndex])
                arrayline.append(solutions[moleculecounter])
            arrayline = numpy.array(arrayline)
            
#            if G.fullBrute == 'yes':
#                specifications = 
#                arrayline = BruteForce(ReferenceData.molecules,  G.objectiveFunctionType)
            
            if G.negativeAnalyzerYorN == 'yes':
                arrayline = NegativeAnalyzer(arrayline,currentReferenceData.reciprocal_matching_correction_values,rawsignalsarrayline,currentReferenceData.molecules,G.objectiveFunctionType)
            
            concentrationline = numpy.zeros(len(arrayline))
            concentrationline[0] = arrayline[0]
            concentrationline[1:] = arrayline[1:]*ExperimentData.conversionfactor

            
            if timeIndex > 1:#all additons after second index
                concentrationsScaledToCOarrayholder = numpy.zeros([len(concentrationsScaledToCOarray[:,0])+1,len(currentReferenceData.molecules)+1])
                concentrationsScaledToCOarrayholder[0:len(concentrationsScaledToCOarray[:,0]),:] = concentrationsScaledToCOarray
                concentrationsScaledToCOarray = concentrationsScaledToCOarrayholder
                concentrationsScaledToCOarray[len(concentrationsScaledToCOarray[:,0])-1,:] = arrayline
                concentrationsarrayholder = numpy.zeros([len(concentrationsarray[:,0])+1,len(currentReferenceData.molecules)+1])
                concentrationsarrayholder[0:len(concentrationsarray[:,0]),:] = concentrationsarray
                concentrationsarray = concentrationsarrayholder
                concentrationsarray[len(concentrationsarray[:,0])-1,:] = concentrationline
            if timeIndex == 1:#additions at second index index
                concentrationsScaledToCOarrayholder = numpy.zeros([2,len(currentReferenceData.molecules)+1])
                concentrationsScaledToCOarrayholder[0,:] = concentrationsScaledToCOarray
                concentrationsScaledToCOarray = concentrationsScaledToCOarrayholder
                concentrationsScaledToCOarray[len(concentrationsScaledToCOarray[:,0])-1,:] = arrayline
                concentrationsarrayholder = numpy.zeros([2,len(currentReferenceData.molecules)+1])
                concentrationsarrayholder[0,:] = concentrationsarray
                concentrationsarray = concentrationsarrayholder
                concentrationsarray[len(concentrationsarray[:,0])-1,:] = concentrationline
            if timeIndex == 0:#first index
                concentrationsScaledToCOarray = arrayline
                concentrationsarray = concentrationline

            
        print('Data Analysis Finished.')
        #show net time for Data Analysis
        G.timeSinceLastCheckpoint = timeit.default_timer() - G.checkpoint
        G.checkPoint = timeit.default_timer()
        print('Data Analysis Time: ', (G.timeSinceLastCheckpoint))
        
        #this section exports and graphs the analyzed signals 
        if G.generatePercentages == 'yes':
            percentagesOutputArray = GeneratePercentages(concentrationsScaledToCOarray)
            ExportXYYYData(G.scaledConcentrationsPercentages, percentagesOutputArray, currentReferenceData.molecules, fileSuffix = G.iterationSuffix)
        ExportXYYYData(G.resolvedScaledConcentrationsOutputName, concentrationsScaledToCOarray, currentReferenceData.molecules, abscissaHeader = "Time", fileSuffix = G.iterationSuffix, dataType = str('scaled'))
        times = concentrationsScaledToCOarray[:,0]#the times are just the first column of the array
        data = concentrationsScaledToCOarray[:,1:]#the data is the whole array except the first column, which is the times
        
        if G.concentrationFinder == 'yes':
            ExportXYYYData(G.concentrationsOutputName, concentrationsarray, currentReferenceData.molecules, abscissaHeader = "Time", fileSuffix = G.iterationSuffix, dataType = 'concentration', units = G.units)
            times = concentrationsarray[:,0]
            data = concentrationsarray[:,1:]
        
        #Graph the concentration/relative signal data
        if G.grapher == 'yes':
            Draw(times, data, currentReferenceData.molecules, G.concentrationFinder, G.units, graphFileName='graphAfterAnalysis', fileSuffix = G.iterationSuffix)

            
    if G.dataSimulation =='yes':
        
        if G.dataAnalysis == 'skip':
            print("Warning: Simulation cannot be run without data")
            sys.exit()
        if G.dataAnalysis == 'load':
            #This function call loads the preprocessed data
            concentrationsScaledToCOarray = ImportAnalyzedData(G.resolvedScaledConcentrationsOutputName)
        #reset timer so that data simiulation can be timed
        G.checkpoint = timeit.default_timer()
        #now the simulated data function uses the answer array and finds what the raw signals would be produced with their signals
        simulateddata = RawSignalsSimulation (concentrationsScaledToCOarray, currentReferenceData.reciprocal_matching_correction_values)
        #Exporting the simulated signals data
        ExportXYYYData(G.simulatedSignalsOutputName, simulateddata, ExperimentData.mass_fragment_numbers, abscissaHeader = "Time", fileSuffix = G.iterationSuffix, dataType = 'simulated')
        #show net time for simulation
        G.timeSinceLastCheckPoint = timeit.default_timer() - G.checkpoint
        G.checkPoint = timeit.default_timer()
        print("Simulation Finished")
        print('Simulation Time: ', (G.timeSinceLastCheckPoint))
        
    CreateLogFile()    
    PopulateLogFile()
    print("LogFile complete")
    
    if G.iterativeAnalysis:
        IterativeAnalysisPostProcessing(ExperimentData, simulateddata, ExperimentData.mass_fragment_numbers, ExperimentDataFullCopy, times, data, currentReferenceData.molecules)

if __name__ == '__main__':
    main()

