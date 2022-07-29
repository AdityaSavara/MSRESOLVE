"""
Created on Tue Jun 26 14:56:24 2018

@author: Andrea
"""


import numpy
import XYYYDataFunctionsSG as DataFunctions
import itertools
import copy
import time
import bisect
#Implict Import: Near the top of the bestMassFragChooser function there is a section which tries to import MSRESOLVE, and if not then separapte imports external dependences such as the Reference Data class.
#The idea of permanently separating that class out of the MSRESOLVE module was considered, e.g. to make a module called MSDataStructures or something like that, but separating has been decided-against at this time.

            
                
#The best mass frag chooser, uses reference patterns for selected molecules to
#deteremine which mass fragements would be the best to monitor in order to 
#correctly identify the composition and concentration of the chosen molecules
#in an MS sample. The best mass fragment chooser will run the reference data
#through a series of checks and the SLS method to determine which combinations
#of mamss fragments will result in a solution via SLS. It will return the N 
#most significant possible combinaitons (if any). Currently, the feature 
#does not work with SLS common due to errors in SLS common.
def bestMassFragChooser(chosenMolecules,
    moleculesLikelihood,
    numberOfMassFragsToMonitor,
    referenceFileName,
    referenceForm,
    significantPeakThreshold=5,
    onTheFlySLS=False,
    keep_N_ValuesInRoughUniquenessCheck=1000,
    keep_N_ValuesInSignificanceFactorCheck=1000,
    finalNumberOfCombinationsToKeep=10,
    exportToFile=True,
    useExtentOfSLSUniqueSolvable=True, printProgress=True,
    minimizeDependencies = False):
    
    #The below code is to allow replacing MSRESOLVE with the "standalone" dependencies file.
    global MSRESOLVE
    if minimizeDependencies == False:
        import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)
        MSRESOLVE.G.calculateUncertaintiesInConcentrations = False #need to turn the uncertainties off for bestMassFragChooser.
    elif minimizeDependencies == True:
        import MSResolveDependenciesForExtentOfSLSsolvable
        MSRESOLVE = MSResolveDependenciesForExtentOfSLSsolvable
        #The minimizeDependencies case does not have any calculateUncertaintiesInConcentrations feature, so that does not need to be turned off for this elif block. 
        
    
    
    progressCounter = 0 #just initializing here.
    
    #Initialize a timer
    start=time.time()
    
    #A solution cannot be found if the number of fragments to monitor is less 
    #than the number of molecules to monitor
    if len(chosenMolecules)>numberOfMassFragsToMonitor: 
        raise ValueError('The number of mass fragments to monitor must be larger than or equal to the number of molecules to monitor')
    

    #The below options are only relevant if SLS is going to be used.  So they do not matter if the person is using useExtentOfSLSUniqueSolvable.
    if useExtentOfSLSUniqueSolvable == False:
        import UserInput as G
        G.applyReferenceMassFragmentsThresholds = 'no' #The applyReferenceMassFragmentsThresholds must be turned off for BestMassFragChooser to work in some cases. significantPeakThreshold should be used instead.
        G.uncertainties_dict = {} #FIXME: This is only a temporary solution.
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
        
        G.excludeMoleculesIfSignificantFragmentNotObserved=G.applyRawSignalThresholds
    #Initialize a ReferenceData class object
    [provided_reference_patterns, electronnumbers, molecules, molecularWeights, 
        SourceOfFragmentationPatterns, sourceOfIonizationData,
        relativeIonizationEfficiencies, 
        moleculeIonizationType, mass_fragment_numbers_monitored, 
        referenceFileName, form] = MSRESOLVE.readReferenceFile(
        referenceFileName, referenceForm)
        
    ReferenceData = MSRESOLVE.MSReference(provided_reference_patterns, 
        electronnumbers, molecules, molecularWeights,
        SourceOfFragmentationPatterns, sourceOfIonizationData,
        relativeIonizationEfficiencies, moleculeIonizationType, 
        mass_fragment_numbers_monitored, referenceFileName=referenceFileName, 
        form=form)

    truncatedReferenceData=copy.deepcopy(ReferenceData)
    
    #Truncate the data so that it matches the molecules chosen for consideration. 
    truncatedReferenceData = MSRESOLVE.trimDataMoleculesToMatchChosenMolecules(
        truncatedReferenceData, chosenMolecules)

    
    # standardize the reference data columns such that the maximum value is 
    #100 and everything else is linearly scaled according that the maximum 
    #value scaling
    truncatedReferenceData.standardized_reference_patterns=MSRESOLVE.StandardizeReferencePattern(
        truncatedReferenceData.provided_reference_patterns,len(truncatedReferenceData.molecules))
    
    #if using the TuningCorrector, then there is another external dependency
    if minimizeDependencies == False:
        import UserInput as G
        if type(G.referenceCorrectionCoefficients) == type({}):# if it's a dictionary we need to make it into a list.
            G.referenceCorrectionCoefficients = [G.referenceCorrectionCoefficients['A'],G.referenceCorrectionCoefficients['B'],G.referenceCorrectionCoefficients['C']]
        truncatedReferenceData.standardized_reference_patterns, ReferenceData.standardized_reference_patterns_tuning_uncertainties = MSRESOLVE.TuningCorrector(
            truncatedReferenceData.standardized_reference_patterns, G.referenceCorrectionCoefficients, G.referenceCorrectionCoefficients_cov,
            G.referenceFileDesiredTuningAndForm, G.referenceFileExistingTuningAndForm, G.tuningCorrection)  
    
    #Removing entire rows of data for mass fragments with all reference
    #intensities below the threshold.
    truncatedReferenceData.standardized_reference_patterns=DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(
        truncatedReferenceData.standardized_reference_patterns,
        startingRowIndex=1,threshold=significantPeakThreshold)

    #we need to repopulate the mass fragments since we're removing columns/rows.
    truncatedReferenceData.provided_mass_fragments=truncatedReferenceData.standardized_reference_patterns[:,0]
    
    #Set all values below the significance threshold to zero (This is in case 
    #one molecule had a value that was bigger, while another molecule didn't. 
    #Since we want an array with only the significant values, we're removing the rest.
    #First create a new array, then do the filtering. The threshold filtering 
    #line has =0 at the end.  the syntax is a bit strange, but it works.
    truncatedReferenceData.standardized_reference_patterns_thresholdFiltered = truncatedReferenceData.standardized_reference_patterns*1.0 #lazy copy. 
    truncatedReferenceData.standardized_reference_patterns_thresholdFiltered[
        truncatedReferenceData.standardized_reference_patterns_thresholdFiltered<significantPeakThreshold]=0
       
    #FIXME: Currently the rest of the code uses truncatedReferenceData.standardized_reference_patterns, 
    #so now we're going to set it equal to do that before the correction values.
    #it seems like this is not the correct order to do things, since now the 
    #correction values are distorted. But correctionValues obtain uses the data 
    #objects inside the Reference object.
    #So until someone thinks further about the algorithm, it might have to stay this way.
    truncatedReferenceData.standardized_reference_patterns = truncatedReferenceData.standardized_reference_patterns_thresholdFiltered
    
    #The below options are only relevant if SLS is going to be used.  So they do not matter if the person is using useExtentOfSLSUniqueSolvable.
    if minimizeDependencies == False:
        #Create the correction values to be used in the SLS method
        truncatedReferenceData.correction_values, truncatedReferenceData.correction_values_relative_uncertainty = MSRESOLVE.CorrectionValuesObtain(truncatedReferenceData)

    
    #Need to reorder the list of molecular likelihoods so they match the class
    #molecules object. In order to reorder the molecular likelihoods, a 
    #dictionary was created to keep track of the molecular likelihoods for each
    #of the molecules
    moleculesLikelihoodDict=dict([(moleculeName,moleculesLikelihood[moleculeIndex]) for moleculeIndex, moleculeName in enumerate(chosenMolecules)])
    reorderedMoleculesLikelihood=[]
    #A list of the reordered molecular likelihoods had to be created in order to
    #allow for multiplication in the following step.
    for moleculeName in truncatedReferenceData.molecules:
        reorderedMoleculesLikelihood.append(moleculesLikelihoodDict[moleculeName])

        
    #Fabricate a data for a single abscissa value based on the probability of 
    #each molecule being present and the reference intensities.
    fabricatedDataTemp=numpy.dot(
        truncatedReferenceData.standardized_reference_patterns[:,1:],
        reorderedMoleculesLikelihood)
    
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
    truncatedReferenceDataas1sand0s=(
        truncatedReferenceData.standardized_reference_patterns!=0)/1.0

    #A parallel array to the mass fragment numbers monitored is initialized as 
    #zeros. This array is altered during the loop through all possible mass 
    #fragment combinations so that a 1 reprsents a chosen molecule.
    currentMassFragmentsas1sand0s = numpy.zeros(len(truncatedReferenceData.provided_mass_fragments))    

    #If any mass fragment combination results in a completely filled reference
    #pattern, this mass fragment pattern could potentially be used for an 
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
 
    #Trim the mass number column from 
    #truncatedReferenceData.standardized_reference_patterns
    intensityMatrix = truncatedReferenceData.standardized_reference_patterns[:,1:]

    #Using intensityMatrix calculate the signficance of every 
    #intensity value. These are stored in significance matrix which has the 
    #same shape as intensityMatrix
    significanceMatrix = generateSignificanceMatrix(
        intensityMatrix, moleculesLikelihood, minThreshold=significantPeakThreshold)
    
    #we need a version that is pristine before the later parts of the algorithm.
    MSRESOLVE.currentReferenceData=copy.deepcopy(truncatedReferenceData)
    
    #ccalculate the number of combinations to consider, because if it is large we will force the progress bar to turn on.
    import math
    k_c = numberOfMassFragsToMonitor
    n_c = len(truncatedReferenceData.provided_mass_fragments)
    combinationsToConsider = int(math.factorial(n_c)/  (math.factorial(n_c-k_c)*math.factorial(k_c)))
    if combinationsToConsider > 1E6:
        print("There will be more than 1 million combinations explored, print progress is being set to true to show a progress bar.")
        printProgress=True
    if printProgress==True:
            try:
                from tqdm import tqdm #this is a module for a progress bar. It's a bit weird that one imports the function with same name as module.
                t = tqdm(total=combinationsToConsider, position = 0) #position = 0 shouldn't be required, but helps prevent buggy bar behaviour.
                tqdm_Failed = False
            except:
                tqdm_Failed = True
                print("Warning: tqdm import or initialization failed. Progress will be printed in lines instead of a progress bar. Use pip tqdm to install the dependency required for the progress bar.")
                

    
    if useExtentOfSLSUniqueSolvable:
        
        #List store the best solvabilities, same shape and with
        #entries to correspond to topBestMassFragments
        #will store 2-tuples of (extentSolvable, summedSignificance)
        #for every combination in topBestMassFragments
        topSolvabilitiesAndSignificance = []
        
        #Loop through all of the possible mass fragment combinations
        #Assume that the masses within each combination are ordered
        for massFragCombination in itertools.combinations(
            truncatedReferenceData.provided_mass_fragments,
            numberOfMassFragsToMonitor):
            if printProgress==True:
                #combinations are “k items from a set of n... Find all the ways to pick (k) numberOfMassFragsToMonitor people from (n) provided_mass_fragments, and divide by the k! (numberOfMassFragsToMonitor!) variants”
                #general formula  is n!/(n-k)!k!
                progressCounter += 1
                if tqdm_Failed == False:
                    t.update(1) #this is updating the progress bar.
                else: 
                    print(progressCounter, "out of", combinationsToConsider)
                
            
            #Alter the currentMassFragmentsas1sand0s array so that 1s represent the
            #location of any current mass fragments. Loops through the monitored
            #mass fragments for those found in the mass fragment combination.
            for index, massFragment in enumerate(
                truncatedReferenceData.provided_mass_fragments):
                #If the mass fragment is found in the mass fragment combination, 
                #that location in currentMassFragmentsas1sand0s 
                #needs to be set to 1
                if massFragment in massFragCombination:
                    currentMassFragmentsas1sand0s[index]=1
                #For the use of the same variable through
                #multiple iterations, the other values in the array 
                #need to be set back to zero. 
                else:         
                    currentMassFragmentsas1sand0s[index]=0

            MSRESOLVE.currentReferenceData.mass_fragment_numbers_monitored=(
                massFragCombination)
            
            #Reference data arrays need to be multiplied by the 
            #currentMassFragmentsas1sand0s in order to only keep the data that is
            #representative of the current mass fragments.
            currentRefDataas1sand0s=numpy.multiply(
                truncatedReferenceDataas1sand0s.T,
                currentMassFragmentsas1sand0s).T
                
            currentFragReferencePattern=numpy.multiply(
                truncatedReferenceData.standardized_reference_patterns.T,
                currentMassFragmentsas1sand0s).T   
            
            #Sums are taken across each molecule to make checks and additional 
            #caluculations more effecient. The axis=0 argument in the sum funciton
            #creates a 1-D array contianing the sums for each molecule
            rowSumsList=numpy.sum(currentRefDataas1sand0s, axis=0)
            
            #The if statement passes the mass combinations that have at least one
            #signal for each molecule, but are not entirely filled with signals.
            #If the function fails, nothing is done and there is no else statment.
            [passesRowsSumChecks,
                allOverlappingPatterns]= MSRESOLVE.passesRowsSumChecks(
                    rowSumsList, massFragCombination, allOverlappingPatterns, numberOfMassFragments=numberOfMassFragsToMonitor)
            
            if passesRowsSumChecks:
                
                #1-D array of mass fragment indices in 
                #currentMassFragmentsas1sand0s that correspond to 
                #those in massFragCombination
                relevantMassFragsIdx = numpy.where(
                    currentMassFragmentsas1sand0s != 0)[0]

                #We need the matrix of intensities as ones and zeros, that 
                #corresponds to the mass fragments and species we are currently 
                #investigating
                relevantIntensities1sand0s = currentRefDataas1sand0s[
                    relevantMassFragsIdx,1:]
                
                #We also need to remove the irrelevant mass fragment rows 
                #from the significanceMatrix
                relevantSignificanceMatrix = significanceMatrix[
                    relevantMassFragsIdx,:]
                
                #Perform a simplified SLS to check how well this combination 
                #performs.
                (extentSolvable, summedSignificance,
                solvedSpecies, massFragsUsed) = ExtentOfSLSUniqueSolvable(
                    massFragCombination, chosenMolecules,
                    relevantIntensities1sand0s, relevantSignificanceMatrix)
                    
                #Pass to storeAndPop to see how this combination compares 
                #to the others.
                [topSolvabilitiesAndSignificance,
                    topBestMassFragments, listUpdated] = MSRESOLVE.storeAndPop(
                    topSolvabilitiesAndSignificance,
                    (extentSolvable, summedSignificance),
                    topBestMassFragments,
                    massFragCombination,
                    finalNumberOfCombinationsToKeep,
                    optimumType="Maximum"
                    )
                #print(massFragCombination, extentSolvable, summedSignificance)

                #Keep track of the best combination
                if massFragCombination == topBestMassFragments[0]:
                    bestMassFragReference = currentFragReferencePattern[
                        relevantMassFragsIdx,:]
        
        #End for loop over all combinations of mass frags
        print("finished the for loop")
    elif useExtentOfSLSUniqueSolvable == False: #when useExtentOfSLSUniqueSolvable is False, actual SLS is performed.
        #Loop through all of the possible mass fragment combinations            
        for massFragCombination in itertools.combinations(truncatedReferenceData.provided_mass_fragments,numberOfMassFragsToMonitor):    
            if printProgress==True:
                progressCounter += 1
                if tqdm_Failed == False:
                    t.update(1) #this is updating the progress bar.
                else: 
                    print(progressCounter, "out of", combinationsToConsider)
            #Alter the currentMassFragmentsas1sand0s array so that 1s represent the
            #location of any current mass fragments. Loops through the monitored
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
            [passesRowsSumChecks,allOverlappingPatterns]= MSRESOLVE.passesRowsSumChecks(rowSumsList, massFragCombination, allOverlappingPatterns, numberOfMassFragments=numberOfMassFragsToMonitor)
                
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
                            SLSReturnRemainingMolecules=MSRESOLVE.SLSMethod(truncatedReferenceData.molecules,currentFragReferencePatternTruncated[:,1:],correctionValuescurrentTruncated,rawsignalsarrayline,timeIndex,[],[],truncatedReferenceData.molecules,DataRangeSpecifierlist,SLSChoices,massFragCombination,G.permutationNum,[],G.objectiveFunctionType,fabricatedAbscissa[timeIndex],maxPermutations=100001,bestMassFragChooser=True)

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
                    SLSReturnRemainingMolecules=MSRESOLVE.SLSMethod(
                        truncatedReferenceData.molecules,
                        currentFragReferencePatternTruncated[:,1:],
                        correctionValuescurrentTruncated,rawsignalsarrayline,
                        timeIndex,[],[],truncatedReferenceData.molecules,
                        DataRangeSpecifierlist,SLSChoices,massFragCombination,
                        G.permutationNum,[],G.objectiveFunctionType,
                        fabricatedAbscissa[timeIndex],maxPermutations=100001,
                        bestMassFragChooser=True)

                    #If there is no unsolved molecules, append this combinaiton to the list of top best mass fragments.
                    #If order to prevent the storage of a large amount of combinations, only an N number of fragments are stored.
                    #N is set from the finalNumberOfCombinationsToKeep user specified variable. The N number are stored based on
                    #the largest sum of significance values.
                    #FIXME: The 1 is present instead of a 0 due to a bug in the SLS unique method
                    if len(SLSReturnRemainingMolecules)==0:
                        [topSignificanceSumsForMassFrags, topBestMassFragments,valueStored]=MSRESOLVE.significanceFactorCheck(currentFragReferencePattern[:,1:] ,topSignificanceSumsForMassFrags,topBestMassFragments, massFragCombination, finalNumberOfCombinationsToKeep, reorderedMoleculesLikelihood)
                        if massFragCombination==topBestMassFragments[0]: 
                            bestMassFragReference=currentFragReferencePatternTruncated

    #Everything below is exporting etc. and happens regardless of useExtentOfSLSUniqueSolvable.
    if printProgress==True:
        if tqdm_Failed == False:
            t.close() #this will close the progress bar.
    end=time.time()
    #The time is kept for printing purposes
    totalTime=end-start
    ''' Here is where the time is printed '''
    print("Time Taken:", totalTime)
    
    #If there is a solvable set of mass fragments, 
    #a reference data set containing only the best mass fragment 
    #combination and the selected molecules will be exported. 
    if exportToFile==True:
        if len(topBestMassFragments)!=0:
            commentsLine = copy.copy(truncatedReferenceData.SourceOfFragmentationPatterns)
            for elemIndex in range(len(commentsLine)):
                commentsLine[elemIndex] = ''
            commentsHeader = numpy.append('#ComentsLine',commentsLine)        
            moleculesHeader=numpy.append('Molecules',truncatedReferenceData.molecules)
            if minimizeDependencies == False: #below are only populated if we are using the full MRESOLVE
                electronHeader=numpy.append('Electron Numbers',truncatedReferenceData.electronnumbers)
                knownMoleculeIonizationTypeHeader = numpy.append('moleculeIonizationType',truncatedReferenceData.moleculeIonizationType)
                relativeIonizationEfficienciesHeader = numpy.append('relativeIonizationEfficiencies',truncatedReferenceData.relativeIonizationEfficiencies)
                fragmentationSourceHeader = numpy.append('SourceOfFragmentationPatterns',truncatedReferenceData.SourceOfFragmentationPatterns)
                ionizationSourceHeader = numpy.append("sourceOfIonizationData",truncatedReferenceData.sourceOfIonizationData)        
                massHeader=numpy.append('Molecular Mass', truncatedReferenceData.molecularWeights)
                #The header will be stacked before saving the data array. 
                #In order to get the multiple headers present, an array will be used 
                #to contain all of the information to allow for stacking.
                fullHeaderArray=numpy.array([
                    commentsHeader, moleculesHeader, electronHeader,
                    knownMoleculeIonizationTypeHeader,
                    relativeIonizationEfficienciesHeader, fragmentationSourceHeader,
                    ionizationSourceHeader, massHeader])
            if minimizeDependencies == True:
                fullHeaderArray=numpy.array([
                    commentsHeader, moleculesHeader])                    
            #In order to export the proper header using the stacking, the abscissa header is set as blank so it doesn't cause anything to be written where the full header should go.
            MSRESOLVE.ExportXYYYData('bestMassFragReference.csv', 
                bestMassFragReference,
                fullHeaderArray, 
                abscissaHeader='')   
        with open('bestMassFragCombinations.txt', 'w') as bmfcTop:
            if useExtentOfSLSUniqueSolvable == True:
                ordering = "The setting of useExtentOfSLSUniqueSolvable was used. The best combination is the last in this file. \n"
            if useExtentOfSLSUniqueSolvable == False:
                ordering = "The best combination is the first in this file. \n"
            bmfcTop.write(ordering)
            bmfcTop.write("Time taken:" +str(totalTime) + '\n')
            bmfcTop.write(str(topBestMassFragments))
            

    else:
        print('There are no SLS Solvable mass fragment combinations')
    return topBestMassFragments, allOverlappingPatterns


def ExtentOfSLSUniqueSolvable(massFrags, species, 
    onesAndZeros, significanceMatrix):
    """
    Arguments:
        massFrags - An n by 1 numpy array that holds the mass fragment 
                    masses. This is an index on the vertical axis of 
                    the onesAndZeros matrix.
    
        species - An m by 1 numpy array that holds the names of the
                  species being considered as strings. This is an index
                  for the horizontal axis of the onesAndZeros matrix.
                    
        onesAndZeros - This is the n by m 2-d numpy array of one and zeros.
                       Each row corresponds to a particular mass fragment
                       and each column to a species. Should be 
                       cleaned by this point to remove irrelevant
                       mass fragment rows.
                       
        significanceMatrix - A n by m 2-d numpy array that mirrors onesAndZeros.
                             The elements of this matrix indicate the 
                             significance of the corresponding intensity. This 
                             is used to select the best unique value/row in 
                             situations where more than one are available.
                       
    Returns a 3-tuple (extentSolvable, solvedSpecies, massFragsUsed):
    
        extentSolvable - An integer that indicates how many species were 
                         solved through unique SLS. It is the primary 
                         metric to judge the suitability of a particular
                         mass fragment combination
                         
        summedSignificance - An scalar that is the sum of the significances 
                             of the intensities used for SLS. It serves
                             as a secondary metric for selecting mass 
                             fragement combinations. 
                             
        
        solvedSpecies - 1-d numpy Array, shadows the species argument to this 
                        function. Initially set to all zeros. When a species is 
                        solvable through unique SLS we change the appropriate 
                        element from 0 to 1 to indicate that fact.
                         
        massFragsUsed - A list that specifies, in order, which 
                        mass fragments were used for the SLS unique 
                        solutions.
    """
        
    extentSolvable = 0 
    summedSignificance = 0
    solvedSpecies = numpy.zeros(len(species))
    massFragsUsed = []
    
    for molecule in range(len(species) + 1):
            
        rowSum = numpy.sum(onesAndZeros, axis=1)
    
        #Is there a 1 in the rowSum numpy array?
        #Note: See https://stackoverflow.com/questions/432112/is-there-a-numpy-function-to-return-the-first-index-of-something-in-an-array
        #For info on numpy.where function
        onesIndices = numpy.where(rowSum == 1)[0]
        if len(onesIndices) > 0: #If there is at least 1, then there is at least one unique row.
                        
            #get a list of the unique value coordinates within onesAndZeros
            #where those coordinates are stored as tuples (row,col)
            #Below, for (rowOfUnique, colOfUnique), the rowOfUnique is the index of the unique mass fragment *within* the massFrags combination provided, and colOfUnique is the index of them molecule.
            uniqueCoordinates = []
            for rowOfUnique in onesIndices:
                colOfUnique = numpy.where(onesAndZeros[rowOfUnique,:] == 1)
                uniqueCoordinates.append((rowOfUnique, colOfUnique))
            
            #If there is more than one unique coordinate
            #pick the best one based on the significanceMatrix
            #Below, for (r,c), the r is the index of the unique mass fragment *within* the massFrags combination provided, and c is the index of them molecule.
            if len(uniqueCoordinates) > 1:
                coordSignificances = [significanceMatrix[r,c] for 
                    (r,c) in uniqueCoordinates]
                
                #reorder uniqueCoordinates by significance
                orderedUniqueCoordinates = [
                    x for (molecule,x) in sorted(zip(
                        coordSignificances, uniqueCoordinates),
                        key=lambda pair:pair[0])]
                
                #Take the coordinates with the highest significance
                bestUniqueRowIdx = orderedUniqueCoordinates[-1][0]
                bestUniqueColIdx = orderedUniqueCoordinates[-1][1]
            
            # if there is only a single unique value then it's the best
            else:
                bestUniqueRowIdx = uniqueCoordinates[0][0]
                bestUniqueColIdx = uniqueCoordinates[0][1]
            
            
            
            #indicate that the mass fragment corresponding to this row
            #was used. 
            massFragsUsed.append(massFrags[bestUniqueRowIdx])
            
            #indicate that this species is now solved for
            solvedSpecies[bestUniqueColIdx] = 1
            
            #Set the column that corresponds to the species we have just 
            #solved to zeros. Now we can sum rows again
            onesAndZeros[:,bestUniqueColIdx] = 0
        
            extentSolvable += 1
            summedSignificance += float(
                significanceMatrix[bestUniqueRowIdx, bestUniqueColIdx])
            #print(massFrags)
            #print(extentSolvable, solvedSpecies, summedSignificance, massFragsUsed)
                
        #Else there was no unique one in the rowSum. 
        #We can proceed no further with unique SLS.
        else:
            #print(extentSolvable, solvedSpecies, summedSignificance, massFragsUsed)
            return (extentSolvable, summedSignificance,
                    solvedSpecies, massFragsUsed)   
                      
                      
def generateSignificanceMatrix(intensityMatrix, moleculesLikelihood=None, minThreshold=None, maxIntensityPossible=100): 
    """
    Given a matrix of reference intensities, with a horizontal
    axis of species and a vertical axis of mass fragments (excluding 
    label columns/rows), this function
    determines the significance of each intensity in the matrix
    
    Arguments:
        intensityMatrix - An n by m 2-d numpy array that contains the 
                          reference intensities corresponding to 
                          n mass fragments and m species.
                          
      moleculesLikelihood - this is a one-dimensional array of values between 0 and 1
                         which are basically the probabilities of each of the m species.
                         
        minThreshold -  this is basically the lower cutoff for intensity of a fragment, and is necessary
                         to have in order to make the best possible significance matrix.
                         
     minThreshold -  this is basically the maximum intensity possible of a fragment, and is necessary
                         to have in order to make the best possible significance matrix.
                          
    Returns:
        significanceMatrix - An n by m 2-d numpy array that mirrors 
                             intensityMatrix. Each element in significanceMatrix
                             is the significance of the corresponding intensity 
                             from intensityMatrix.
    """

    numberOfMolecules = len(intensityMatrix[0])
    if moleculesLikelihood == None: #if it is not provided, assume a likelihood of 1 for each molecule.
        numpy.ones(numberOfMolecules)
    
    #We need to form  2-d array that shadows the intensity data
    #from intensityMatrix.
    #Initially we will set it to zeros and then populate in loop.
    significanceMatrix = numpy.zeros(
        intensityMatrix.shape)
        
    #loop through columns (i.e. species)
    ##TODO this loop was designed to work with the ElemSignificanceCalculator
    ## and IndElemSignificanceCalculator functions in MSRESOLVE that are
    ## apparently designed for a different application.
    ## If it proves to be a performance bottleneck (which I doubt since 
    ## we aren't looping here), we can rewrite those to cut down on the 
    ## required number of row sums. As it is now we sum a row for every element
    ## in each column, and we do this for every column.
    for (columnIdx, _) in enumerate(significanceMatrix.T):
        #determine, for every element in this column (i.e. for this
        #species) the element's significance
        colSignificances = MSRESOLVE.ElemSignificanceCalculator(
            intensityMatrix, columnIdx, moleculesLikelihood, 
            minThreshold=minThreshold, maxIntensityPossible=maxIntensityPossible)
            
        #Now record these significances
        significanceMatrix[:,columnIdx] = colSignificances
        
    return significanceMatrix

    
#Optional todo item (probably to never finish)
# finish ranking systems to use update or insert. This will keep a list 
#containing 2 tuples: 1) the mass fragment combinaiton and 2) the values of the
#Objective functions. The values of the objective functions will contain 3
#elementss: 1) the number of remaining molecules left from SLS, 2) the negative
#significance sum [or significance sum in case of ExtentOfSLSUniqueSolvable] 
#and 3) the rough uniqueness value. After the loop though SLS,
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
