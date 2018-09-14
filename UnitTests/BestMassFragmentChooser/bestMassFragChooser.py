"""
Created on Tue Jun 26 14:56:24 2018

@author: Andrea
"""

import MSRESOLVE
import numpy
import XYYYDataFunctionsSG as DataFunctions
import UserInput as G
import itertools
import copy
import time
import bisect

#TODO finish ranking systems to use update or insert. This will keep a list 
#containing 2 tuples: 1) the mass fragment combinaiton and 2) the values of the
#Objective functions. The values of the objective functions will contain 3
#elementss: 1) the number of remaining molecules left from SLS, 2) the negative
#significance sum and 3) the rough uniqueness value. After the loop though SLS,
#the mass fragments will be sorted by objective values in an ascending order.
#The value to insert tuple will contian the mass fragment and the objective
#functions tuple.

def updateOrInsert(storeInObjectiveFunctionValuesList, valueToInsertTuple):
    #If the mass fragment combination is not already present in the list insert
    #it
    if valueToInsertTuple[0] not in storeInObjectiveFunctionValuesList[:][0]:
        storeInObjectiveFunctionValuesList=bisect.insort(storeInObjectiveFunctionValuesList,valueToInsertTuple)
    else: #If the mass fragment combinaiton is already in the list
        #Loop throught the storedInObjectiveFunctionValuesList to determine where
        #The mass fragment combination was
        for storeInObjectiveFunctionValuesIndex, storedInObjectiveFunctionValues in enumerate(storeInObjectiveFunctionValuesList):
            #If the mass fragment combinaiton was in the list, set the
            if valueToInsertTuple[1]==storedInObjectiveFunctionValues[1]:
                storeInObjectiveFunctionValuesList[storeInObjectiveFunctionValuesIndex]=valueToInsertTuple
                
                
#The best mass frag chooser, uses reference patterns for selected molecules to
#deteremine which mass fragements would be the best to monitor in order to 
#correctly identify the composition and concentration of the chosen molecules
#in an MS sample. The best mass fragment chooser will run the reference data
#through a series of checks and the SLS method to determine which combinations
#of mamss fragments will result in a solution via SLS. It will return the N 
#most significant possible combinaitons (if any). Currently, the feature 
#does not work with SLS common due to errors in SLS common.
def bestMassFragChooser(moleculesToMonitor, moleculesLikelihood, numberOfMassFragsToMonitor, referenceFileName, referenceForm,significantPeakThreshold=5, onTheFlySLS=False, keep_N_ValuesInRoughUniquenessCheck=1000, keep_N_ValuesInSignificanceFactorCheck=1000, finalNumberOfCombinationsToKeep=10):
    
    #Initialize a timer
    start=time.time()
    
    #A solution cannot be found if the number of fragmemnts to monitor is less 
    #than the number of molecules to monitor
    if len(moleculesToMonitor)>numberOfMassFragsToMonitor: 
        raise ValueError('The number of mass fragments to monitor must be larger than or equal to the number of molecules to monitor')
    
    #The 4th "theoretical" case of having no limits post-loop SLS does not make sense, because there is a risk of too many combinations (which is why we make limits).
    #So if someone tries to choose that, we force it into on the fly SLS (without limits).
    if (keep_N_ValuesInRoughUniquenessCheck == False and keep_N_ValuesInSignificanceFactorCheck == False) and onTheFlySLS == False:
        onTheFlySLS == True
    
    #Initiialize the data range specifier list. This list contains variables
    #defined in the user input file. It is not of importnace for the best mass
    #frag choser, but is required as an arguement that needs to be passed to
    #the SLS method.
    DataRangeSpecifierlist = [G.dataRangeSpecifierYorN, G.signalOrConcentrationRange,
						  G.csvFile, G.moleculesToRestrict, G.csvFileName,G.dataUpperBound,
						  G.dataLowerBound, G.bruteIncrements, G.permutationNum]
    
    #Initializes the SLS choices variable. As of now, the the best mass frag
    #chooser can only run with unique SLS do to errors in the SLS common
    #method.
    SLSChoices = ['unique', G.slsFinish, G.distinguished]
    
    G.uniqueOrCommon='unique'
    
    G.excludeMoleculesIfSignificantFragmentNotObserved=G.rawSignalThresholdMethod
    #Initialize a ReferenceData class object
    [provided_reference_patterns, electronnumbers, molecules, molecularWeights, sourceInfo, mass_fragment_numbers_monitored, referenceFileName, form]=MSRESOLVE.readReferenceFile(referenceFileName, referenceForm)
    ReferenceData = MSRESOLVE.MSReference(provided_reference_patterns, electronnumbers, molecules, molecularWeights, sourceInfo, mass_fragment_numbers_monitored, referenceFileName=referenceFileName, form=form)

    truncatedReferenceData=copy.deepcopy(ReferenceData)
    
    #Truncate the data so that it matches the molecules chosen for consideration. #TODO: should probably change moleculesToMonitor to be G.chosenMolecules
    
    truncatedReferenceData = MSRESOLVE.trimDataMoleculesToMatchChosenMolecules(truncatedReferenceData, moleculesToMonitor)

    
    # standardize the reference data columns such that the maximum value is 
    #100 and everything else is linearly scaled according that the maximum 
    #value scaling
    truncatedReferenceData.standardized_reference_patterns=MSRESOLVE.StandardizeReferencePattern(truncatedReferenceData.provided_reference_patterns,len(truncatedReferenceData.molecules))
    truncatedReferenceData.standardized_reference_patterns = MSRESOLVE.CorrectionValueCorrector(truncatedReferenceData.standardized_reference_patterns, G.referenceCorrectionCoefficients,
                                               G.referenceLiteratureFileName, G.referenceMeasuredFileName,
                                               G.measuredReferenceYorN)  
    
    #Removing entire rows of data for mass fragments with all reference
    #intensities below the threshold.
    truncatedReferenceData.standardized_reference_patterns=DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(truncatedReferenceData.standardized_reference_patterns,startingRowIndex=1,threshold=significantPeakThreshold)

    #we need to repopulate the mass fragments since we're removing columns/rows.
    truncatedReferenceData.provided_mass_fragments=truncatedReferenceData.standardized_reference_patterns[:,0]
    
    #Set all values below the significance threshold to zero (This is in case one molecule had a value that was bigger, while another molecule didn't. Since we want an array with only the significant values, we're removing the rest.
    #First create a new array, then do the filtering. The threshold filtering line has =0 at the end.  the syntax is a bit strange, but it works.
    truncatedReferenceData.standardized_reference_patterns_thresholdFiltered = truncatedReferenceData.standardized_reference_patterns*1.0 #lazy copy. 
    truncatedReferenceData.standardized_reference_patterns_thresholdFiltered[truncatedReferenceData.standardized_reference_patterns_thresholdFiltered<significantPeakThreshold]=0
       
    #FIXME: Currently the rest of the code uses truncatedReferenceData.standardized_reference_patterns, so now we're going to set it equal to do that before the correction values.
    #it seems like this is not the correct order to do things, since now the correction values are distorted. But correctionValues obtain uses the data objects inside the Reference object.
    # so until someone thinks further about the algorithm, it might have to stay this way.
    truncatedReferenceData.standardized_reference_patterns = truncatedReferenceData.standardized_reference_patterns_thresholdFiltered
  
    truncatedReferenceData.correction_values = MSRESOLVE.CorrectionValuesObtain(truncatedReferenceData) 
    #Create the correction values to be used in the SLS method
    truncatedReferenceData.correction_values = MSRESOLVE.CorrectionValuesObtain(truncatedReferenceData)

    
    #Need to reorder the list of molecular likelihoods so they match the class
    #molecules object. In order to reorder the molecular likelihoods, a 
    #dictionary was created to keep track of the molecular likelihoods for each
    #of the molecules
    moleculesLikelihoodDict=dict([(molecule,moleculesLikelihood[moleculeIndex]) for moleculeIndex, molecule in enumerate(moleculesToMonitor)])
    reorderedMoleculesLikelihood=[]
    #A list of the reordered molecular likelihoods had to be created in order to
    #allow for multiplication in the following step.
    for molecule in truncatedReferenceData.molecules:
        reorderedMoleculesLikelihood.append(moleculesLikelihoodDict[molecule])

        
    #Fabricate a data for a single abscissa value based on the probability of 
    #each molecule being present and the reference intensities.
    fabricatedDataTemp=numpy.dot(truncatedReferenceData.standardized_reference_patterns[:,1:], reorderedMoleculesLikelihood)
    
    #For SLS to work properly, the data cannot be a 1-D array. A duplicte row 
    #is used to prevent any problems
    fabricatedData=numpy.ones((2,len(fabricatedDataTemp)))
    fabricatedData[0],fabricatedData[1]=fabricatedDataTemp,fabricatedDataTemp
    
    #fabricates abscissa values to be used in the SLS method. The numbers 1 and
    #2 are arbitrary times and are used later on for SLS consistency
    fabricatedAbscissa=numpy.array([1,2])
    
    #The parallel reference data array is made to represent the locations of 
    #any data and zeros. A zero is in the location of zeros in the reference
    #data. 1s represent the location of any non-zero value
    truncatedReferenceDataas1sand0s=(truncatedReferenceData.standardized_reference_patterns!=0)/1.0

    #A parallel array to the mass fragment numbers monitored is initialized as 
    #zeros. This array is altered during the loop through all possible mass 
    #fragment combinations so that a 1 reprsents a chosen molecule.
    currentMassFragmentsas1sand0s = numpy.zeros(len(truncatedReferenceData.provided_mass_fragments))    

    #If any mass fragment combinaiton results in a completely filled reference
    #pattern, this mass frgment pattern could potentially be used for an 
    #inverse method. These combinations will be stored in the following list.
    allOverlappingPatterns=[]
    
    #Initialize lists to be used to store mass frgament combinaitons and 
    #corresponding list of objetcive fuction ranking values
    topRoughUniquenessSumsList=[]
    topMassFragCombinationsRUList=[]
    largestMagnitudeSigFactorSumsList=[]
    topMassFragCombinationsSFList=[]
    topBestMassFragments=[]
    topSignificanceSumsForMassFrags=[]
    bestMassFragReference=None
    
    MSRESOLVE.currentReferenceData=copy.deepcopy(truncatedReferenceData)
    
    #Loop through all of the possible mass fragment combinations
    for massFragCombination in itertools.combinations(truncatedReferenceData.provided_mass_fragments,numberOfMassFragsToMonitor):    
        
        #Alter the currentMassFragmentsas1sand0s array so that 1s represent the
        #locaiton of any current mass fragments. Loops through the monitored
        #mass fragments for those found in the mass fragment combination.
        for index, massFragment in enumerate(truncatedReferenceData.provided_mass_fragments):
            #If the mass fragment is found in the mass fragment combinaiton, 
            #that locaiton in currentMassFragmentsas1sand0s needs to be set to 1
            if massFragment in massFragCombination:
                currentMassFragmentsas1sand0s[index]=1
            else: #For the use of the same variable through multiple iterations, the other values in the array need to be set back to zero. 
                currentMassFragmentsas1sand0s[index]=0
        MSRESOLVE.currentReferenceData.mass_fragment_numbers_monitored=massFragCombination
        #Reference data arrays need to be multiplied by the 
        #currentMassFragmentsas1sand0s in order to only keep the data that is
        #representative of the current mass fragments.
        currentRefDataas1sand0s=numpy.multiply(truncatedReferenceDataas1sand0s.T,currentMassFragmentsas1sand0s).T
        currentFragReferencePattern=numpy.multiply(truncatedReferenceData.standardized_reference_patterns.T,currentMassFragmentsas1sand0s).T

        #Sums are taken across each molecule to make checks and additional 
        #caluculations more effecient. The axis=0 argument in the sum funciton
        #creates a 1-D array contianing the sums for each molecule
        rowSumsList=numpy.sum(currentRefDataas1sand0s, axis=0)
       
        
        #The if statement passes the mass combinations that have at least one
        #signal for each molecule, but are not entirely filled with signals.
        #If the function fails, nothing is done and there is no else statment.
        [passesRowsSumChecks,allOverlappingPatterns]= MSRESOLVE.passesRowsSumChecks(rowSumsList, massFragCombination, allOverlappingPatterns)
            
        if passesRowsSumChecks:    
            #We start as false, and change to True after a value has been stored in a top list.
            valueStoredInTopList=False 
            
            #If the rough uniqueness check is desired, run the rough uniqueness
            #check.
            if keep_N_ValuesInRoughUniquenessCheck!=False:
                [topRoughUniquenessSumsList,topMassFragCombinationsRUList,valueStoredInRUTopList] = MSRESOLVE.roughUniquenessCheck(rowSumsList, topRoughUniquenessSumsList,topMassFragCombinationsRUList, keep_N_ValuesInRoughUniquenessCheck, massFragCombination)
                #For on the fly SLS, we need to know if the if a combinaiton 
                #was stored in the top mass fragment combinaitons of the rough
                #uniquness check
                if valueStoredInRUTopList: valueStoredInTopList=True

            #If the significance factor check is desired, run the signnificance
            #factor check.
            if keep_N_ValuesInSignificanceFactorCheck!=False:
                [largestMagnitudeSigFactorSumsList,topMassFragCombinationsSFList, valueStoredInSFTopList]=MSRESOLVE.significanceFactorCheck(currentFragReferencePattern[:,1:],largestMagnitudeSigFactorSumsList,topMassFragCombinationsSFList, massFragCombination, keep_N_ValuesInSignificanceFactorCheck, reorderedMoleculesLikelihood)
                #For on the fly SLS, we need to know if the if a combinaiton 
                #was stored in the top mass fragment combinaitons of the 
                #significance factor check
                if valueStoredInSFTopList: valueStoredInTopList=True
            
            #Need to only use the fabricated data of the current mass fragments.
            #This just sets reference data intensites for non-current mass
            #fragments to zero. Truncation of these zero rows is not 
            #necessay due to the raw signals array maker used later.                
            fabricatedDatacurrent=numpy.multiply(fabricatedData, currentMassFragmentsas1sand0s)

            #The correction value matrix is multiplied by the 
            #currentMassFragmentsas1sand0s. The correction values of the non-
            #current mass fragments are set to zeros.
            correctionValuescurrent=numpy.multiply(truncatedReferenceData.correction_values, currentMassFragmentsas1sand0s)
            
            #For SLS to funciton properly, the rows of non-current mass 
            #fragments must be removed from the correction values and the 
            #current reference.
            correctionValuescurrentTruncated=(DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(correctionValuescurrent.T))
            
            currentFragReferencePatternTruncated=DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(currentFragReferencePattern)
            #there are two cases of on the fly sls. One with limits, and one without. #This is cases 1 and 3.
            #on the fly SLS with no limiting (checks *every* combination, since our limits are not rigorous)
            ##This is the case of onTheFlySLS with limits. We run SLSpossible if the "value" was stored in either list.
            #This case is advantageous compared to post-loop SLS because fewer combinations get stored if SLSpossible rejects some/many combinations.
            if onTheFlySLS:
                #This if statement is not included in the pseudocode. It was split into two parts in the pseudocode, but can be combinded into one to avoid funcitonalizing run SLS
                #Case 3 is represented by the first half of the if statement and case 1 is represented by the second portion
                if (keep_N_ValuesInRoughUniquenessCheck==False and keep_N_ValuesInSignificanceFactorCheck==False) or valueStoredInTopList:
                    
               #Loop through the the first value in the fabricatedAbscissa. Only the first row of the first abscissa value is of interest sice the second row is just a duplicate.
               #There is no need to run through the duplicate. The loop form is kept only for consistency with the actual running of the SLS method in MSRESOLVE and for the 
               #needed arguments of the SLS method.
                    for timeIndex in range(len(fabricatedAbscissa)-1):
                       #Create a raw signal array that is used for the SLS method. This is mimicking the same requirements as in MSRESOLVE
                        rawsignalsarrayline = MSRESOLVE.RawSignalsArrayMaker(massFragCombination,
                                                                   currentFragReferencePattern[:,0],fabricatedDatacurrent,
                                                                   timeIndex,currentFragReferencePattern[:,0])
                        
                        #The SLS method with the best mass frag chooser 
                        #variable set to True will return the unsolved 
                        #molecules
                        SLSReturnRemainingMolecules=MSRESOLVE.SLSMethod(truncatedReferenceData.molecules,currentFragReferencePatternTruncated[:,1:],correctionValuescurrentTruncated,rawsignalsarrayline,timeIndex,[],[],truncatedReferenceData.molecules,DataRangeSpecifierlist,SLSChoices,massFragCombination,G.permutationNum,[],G.bruteOption,fabricatedAbscissa[timeIndex],maxPermutations=100001,bestMassFragChooser=True)

                        #If there is no unsolved molecules, append this combinaiton to the list of top best mass fragments.
                        #If order to prevent the storage of a large amount of combinations, only an N number of fragments are stored.
                        #N is set from the finalNumberOfCombinationsToKeep user specified variable. The N number are stored based on
                        #the largest sum of significance values.
                        #FIXME: The 1 is present instead of a 0 due to a bug in the SLS unique method
                        if len(SLSReturnRemainingMolecules)==0:
                            [topSignificanceSumsForMassFrags, topBestMassFragments,valueStored]=MSRESOLVE.significanceFactorCheck(currentFragReferencePattern[:,1:] ,topSignificanceSumsForMassFrags,topBestMassFragments, massFragCombination, finalNumberOfCombinationsToKeep, reorderedMoleculesLikelihood)
                            if massFragCombination==topBestMassFragments[0]: bestMassFragReference=currentFragReferencePatternTruncated
                            

    if not onTheFlySLS:
        #Create a set of mass fragments to run through from those stored by the limiting checks. There are no duplicate mass fragment combinations in the list.
        limitingChecksMassFrags=list(set(topMassFragCombinationsRUList+topMassFragCombinationsSFList))
        
        #Iterate through the set of mamss fragment combinaitons stored by the limiting checks
        for massFragCombination in limitingChecksMassFrags:
            #Alter the currentMassFragmentsas1sand0s array so that 1s represent the
            #locaiton of any current mass fragments. Loops through the monitored
            #mass fragments for those found in the mass fragment combination.
            for index, massFragment in enumerate(truncatedReferenceData.provided_mass_fragments):
                #If the mass fragment is found in the mass fragment combinaiton, 
                #that locaiton in currentMassFragmentsas1sand0s needs to be set to 1
                if massFragment in massFragCombination:
                    currentMassFragmentsas1sand0s[index]=1
                else: #For the use of the same variable through multiple iterations, the other values in the array need to be set back to zero. 
                    currentMassFragmentsas1sand0s[index]=0

            #Generate the necessary data arrays
            #Reference data arrays need to be multiplied by the 
            #currentMassFragmentsas1sand0s in order to only keep the data that is
            #representative of the current mass fragments.
            currentRefDataas1sand0s=numpy.multiply(truncatedReferenceDataas1sand0s.T,currentMassFragmentsas1sand0s).T
            
            #For SLS to funciton properly, the rows of non-current mass 
            #fragments must be removed from the correction values and the 
            #current reference.
            currentFragReferencePattern=numpy.multiply(truncatedReferenceData.standardized_reference_patterns.T,currentMassFragmentsas1sand0s).T
            currentFragReferencePatternTruncated=DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(currentFragReferencePattern)
            
            #This just sets reference data intensites for non-current mass
            #fragments to zero. Truncation of these zero rows is not 
            #necessay due to the raw signals array maker used later.      
            fabricatedDatacurrent=numpy.multiply(fabricatedData, currentMassFragmentsas1sand0s)     
            #The correction value matrix is multiplied by the 
            #currentMassFragmentsas1sand0s. The correction values of the non-
            #current mass fragments are set to zeros.
            correctionValuescurrent=numpy.multiply(truncatedReferenceData.correction_values, currentMassFragmentsas1sand0s)
            correctionValuescurrentTruncated=(DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(correctionValuescurrent.T))
            
            #Loop through the the first value in the fabricatedAbscissa. Only the first row of the first abscissa value is of interest sice the second row is just a duplicate.
           #There is no need to run through the duplicate. The loop form is kept only for consistency with the actual running of the SLS method in MSRESOLVE and for the 
           #needed arguments of the SLS method.
            for timeIndex in range(len(fabricatedAbscissa)-1):
               #Create a raw signal array that is used for the SLS method. This is mimicking the same requirements as in MSRESOLVE
                rawsignalsarrayline = MSRESOLVE.RawSignalsArrayMaker(massFragCombination,
                                                           currentFragReferencePattern[:,0],fabricatedDatacurrent,
                                                           timeIndex,currentFragReferencePattern[:,0])    
                #The SLS method with the best mass frag chooser 
                #variable set to True will return the unsolved 
                #molecules
                SLSReturnRemainingMolecules=MSRESOLVE.SLSMethod(truncatedReferenceData.molecules,currentFragReferencePatternTruncated[:,1:],correctionValuescurrentTruncated,rawsignalsarrayline,timeIndex,[],[],truncatedReferenceData.molecules,DataRangeSpecifierlist,SLSChoices,massFragCombination,G.permutationNum,[],G.bruteOption,fabricatedAbscissa[timeIndex],maxPermutations=100001,bestMassFragChooser=True)

                #If there is no unsolved molecules, append this combinaiton to the list of top best mass fragments.
                #If order to prevent the storage of a large amount of combinations, only an N number of fragments are stored.
                #N is set from the finalNumberOfCombinationsToKeep user specified variable. The N number are stored based on
                #the largest sum of significance values.
                #FIXME: The 1 is present instead of a 0 due to a bug in the SLS unique method
                if len(SLSReturnRemainingMolecules)==0:
                    [topSignificanceSumsForMassFrags, topBestMassFragments,valueStored]=MSRESOLVE.significanceFactorCheck(currentFragReferencePattern[:,1:] ,topSignificanceSumsForMassFrags,topBestMassFragments, massFragCombination, finalNumberOfCombinationsToKeep, reorderedMoleculesLikelihood)
                    if massFragCombination==topBestMassFragments[0]: bestMassFragReference=currentFragReferencePatternTruncated
                  
    end=time.time()
    #The time is kept for printing purposes
    totalTime=end-start
#If there is a solvable set of mass fragments, a reference data set containing only the best mass fragment combination and the selected molecules will be exported. 
    if len(topBestMassFragments)!=0:
        sourceHeader=numpy.append('Source:',truncatedReferenceData.sourceInfo)
        moleculesHeader=numpy.append('Molecules',truncatedReferenceData.molecules)
        electronHeader=numpy.append('Electron Numbers',truncatedReferenceData.electronnumbers)
        massHeader=numpy.append('Molecular Mass', truncatedReferenceData.molecularWeights)
        #The header will be stacked before saving the data array. In order to get the multiple headers present, an array will be used to contain all of the information to allow for stacking.
        fullHeaderArray=numpy.array([sourceHeader,moleculesHeader,electronHeader,massHeader])
    	#In order to export the proper header using the stacking, the abscissa header is set as blank so it doesn't cause anything to be written where the full header should go.
        MSRESOLVE.ExportXYYYData('bestMassFragReference.csv', bestMassFragReference,fullHeaderArray, abscissaHeader='')    
    else:
        print('There are no SLS Solvable mass fragment combinations')
    return topBestMassFragments, allOverlappingPatterns
