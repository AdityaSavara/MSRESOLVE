import ParsingFunctions as parse
import copy 
'''
parseUserInput parses the variables in the user input file
It passes in G as an argument
This function is designed to serve as a standard for parsing particular variables
#### WARNING: Anything changed to parallelVectorize for chosenMolecules length, or rather chosenMoleculesForParsing, needs to be added to delimitedStringOfVariablesToUnparallelize insde MSRESOLVE in IterativePrepareNextIterationInputFiles. ####
'''
def parseUserInput(currentUserInput):
    #Input Files
    currentUserInput.referencePatternsFileNamesList = parse.listCast(currentUserInput.referencePatternsFileNamesList) #referenceFileName needs to be a list
    currentUserInput.referencePatternsFileNamesList = parse.stripListOfStrings(currentUserInput.referencePatternsFileNamesList)
    currentUserInput.referencePatternsFormsList = parse.listCast(currentUserInput.referencePatternsFormsList) #form needs to be a list
    currentUserInput.referencePatternsFormsList = parse.stripListOfStrings(currentUserInput.referencePatternsFormsList)
    currentUserInput.referencePatternsFormsList = parse.parallelVectorize(currentUserInput.referencePatternsFormsList,len(currentUserInput.referencePatternsFileNamesList)) #form needs to be a list of the same length as referenceFileName
    currentUserInput.referencePatternTimeRanges = parse.listCast(currentUserInput.referencePatternTimeRanges) #RefPatternTimeRanges needs to be a list
    parse.strCheck(currentUserInput.dataToAnalyzeFileName,'dataToAnalyzeFileName') #dataToAnalyzeFileName must be a string
    currentUserInput.dataToAnalyzeFileName = currentUserInput.dataToAnalyzeFileName.strip()
 
    
    #preProcessing, dataAnalysis, dataSimulation, grapher
    parse.strCheck(currentUserInput.preProcessing,'preProcessing')
    parse.strCheck(currentUserInput.dataAnalysis,'dataAnalysis')
    parse.strCheck(currentUserInput.dataSimulation,'dataSimulation')
    parse.strCheck(currentUserInput.grapher,'grapher')
            
    #Time Range
    parse.strCheck(currentUserInput.timeRangeLimit,'timeRangeLimit')
    #Time Ranges are both floats
    if currentUserInput.timeRangeLimit == 'yes':
        currentUserInput.timeRangeStart = float(currentUserInput.timeRangeStart) 
        currentUserInput.timeRangeFinish = float(currentUserInput.timeRangeFinish)  
    
    #Specific molecules/mass fragments
    parse.strCheck(currentUserInput.specificMolecules,'specificMolecules')
    parse.strCheck(currentUserInput.specificMassFragments,'specificMassFragments')
    #Chosen Molecules and Mass Fragments are both lists
    currentUserInput.chosenMoleculesNames = parse.listCast(currentUserInput.chosenMoleculesNames)
    currentUserInput.chosenMassFragments = parse.listCast(currentUserInput.chosenMassFragments)
    
    #chosenMoleculesNames should have leading and trailing whitespaces removed.
    currentUserInput.chosenMoleculesNames = parse.stripListOfStrings(currentUserInput.chosenMoleculesNames)
    #currentUserInput.exp_mass_fragment_numbers and currentUserInput.moleculesNames are the molecules and the mass fragments from the referece data and collected data, respectively
    #Populate chosenMassFragmentsForParsing based on user input option to get a list of mass fragments
    if currentUserInput.specificMassFragments == 'yes': #if yes, use the user's chosen mass fragments
        chosenMassFragmentsForParsing = copy.deepcopy(currentUserInput.chosenMassFragments)
        #If using specificMassFragments, make sure all selected fragments are in the collected data
        parse.compareElementsBetweenLists(currentUserInput.chosenMassFragments,currentUserInput.exp_mass_fragment_numbers,'chosenMassFragments','Mass Fragments from Data')
    elif currentUserInput.specificMassFragments == 'no': #Otherwise use all mass fragments
        chosenMassFragmentsForParsing = copy.deepcopy(currentUserInput.exp_mass_fragment_numbers)
    #Populate chosenMolecules based on user input option to get a list of molecules
    if currentUserInput.specificMolecules == 'yes': #if yes, use the user's chosen moleclues
        chosenMoleculesForParsing = copy.deepcopy(currentUserInput.chosenMoleculesNames)
        #If using specificMolecules, make sure all selected molecules are in the reference data
        if currentUserInput.tuningCorrection == 'no': #If not making a mixed reference pattern, then use the regular moleculesNames object for comparison.
            parse.compareElementsBetweenLists(currentUserInput.chosenMoleculesNames,currentUserInput.moleculesNames,'chosenMolecules','Molecules from Reference Data')
        if currentUserInput.tuningCorrection == 'yes':#If using a making a reference pattern, check the extended moleculesNames list.
            currentUserInput.moleculesNamesExtended = parse.stripListOfStrings(currentUserInput.moleculesNamesExtended)
            parse.compareElementsBetweenLists(currentUserInput.chosenMoleculesNames,currentUserInput.moleculesNamesExtended,'chosenMolecules','Molecules from Reference Data')
    elif currentUserInput.specificMolecules == 'no': #Otherwise use all molecules
        if currentUserInput.tuningCorrection == 'no': #If not making a mixed reference pattern, then use the regular moleculesNames object for comparison.
            currentUserInput.moleculesNames = parse.stripListOfStrings(list(currentUserInput.moleculesNames))
            chosenMoleculesForParsing = copy.deepcopy(currentUserInput.moleculesNames)
        if currentUserInput.tuningCorrection == 'yes':#If using a making a reference pattern, check the extended moleculesNames list.
            currentUserInput.moleculesNamesExtended = parse.stripListOfStrings(currentUserInput.moleculesNamesExtended)
            chosenMoleculesForParsing = copy.deepcopy(currentUserInput.moleculesNamesExtended)
    
    #Molecule Likelihoods and Sensitivity Values are lists with the same length as the number of molecules
    currentUserInput.moleculeLikelihoods = parse.listCast(currentUserInput.moleculeLikelihoods)
    currentUserInput.moleculeLikelihoods = parse.parallelVectorize(currentUserInput.moleculeLikelihoods,len(chosenMoleculesForParsing))
    currentUserInput.sensitivityValues = parse.listCast(currentUserInput.sensitivityValues)
    currentUserInput.sensitivityValues = parse.parallelVectorize(currentUserInput.sensitivityValues,len(chosenMoleculesForParsing))
    
    #Linear Baseline Correction Semi-Automatic variables
    parse.strCheck(currentUserInput.linearBaselineCorrectionSemiAutomatic,'linearBaselineCorrectionSemiAutomatic')
    if currentUserInput.linearBaselineCorrectionSemiAutomatic == 'yes': #if using linear baseline correction semi automatic
        currentUserInput.baselineType = parse.listCast(currentUserInput.baselineType) #Baseline type needs to be a list
        currentUserInput.massesToBackgroundCorrect = parse.listCast(currentUserInput.massesToBackgroundCorrect) #Masses to background correct is a list        
        if len(currentUserInput.massesToBackgroundCorrect) == 0: #If massesToBackgroundCorrect is empty
            currentUserInput.massesToBackgroundCorrect = chosenMassFragmentsForParsing #Use the chosenMassFragments
        #Check that all masses in currentUserInput.massesToBackgroundCorrect are in the collected data
        parse.compareElementsBetweenLists(currentUserInput.massesToBackgroundCorrect,chosenMassFragmentsForParsing,"massesToBackgroundCorrect","chosenMassFragments")
        #Early and Late baseline times are lists
        currentUserInput.earlyBaselineTimes = parse.listCast(currentUserInput.earlyBaselineTimes) 
        currentUserInput.lateBaselineTimes = parse.listCast(currentUserInput.lateBaselineTimes)
        #Early and late baseline times are also the same length as masses to background correct
        currentUserInput.earlyBaselineTimes = parse.parallelVectorize(currentUserInput.earlyBaselineTimes,len(currentUserInput.massesToBackgroundCorrect))
        currentUserInput.lateBaselineTimes = parse.parallelVectorize(currentUserInput.lateBaselineTimes,len(currentUserInput.massesToBackgroundCorrect)) 
    
    #Data Solving Restrictions - Marginal Change Restrictor
    parse.strCheck(currentUserInput.interpolateYorN,'interpolateYorN')
    if currentUserInput.interpolateYorN == 'yes':
        #Marginal Change Restriction and Ignorable Delta Y Threshold are both floats
        currentUserInput.marginalChangeRestriction = float(currentUserInput.marginalChangeRestriction)
        currentUserInput.ignorableDeltaYThreshold = float(currentUserInput.ignorableDeltaYThreshold)
    
    #Data Solving Restrictions - Brute Solving Restrictions
    parse.strCheck(currentUserInput.dataRangeSpecifierYorN,'dataRangeSpecifierYorN')
    parse.strCheck(currentUserInput.signalOrConcentrationRange,'signalOrConcentrationRange')
    parse.strCheck(currentUserInput.csvFile,'csvFile')
    parse.strCheck(currentUserInput.csvFileName,'csvFileName')
    #Data  Upper/Lower Bound are both lists
    currentUserInput.dataLowerBound = parse.listCast(currentUserInput.dataLowerBound)
    currentUserInput.dataUpperBound = parse.listCast(currentUserInput.dataUpperBound)
    currentUserInput.bruteIncrements = parse.listCast(currentUserInput.bruteIncrements) #increments is a list
    currentUserInput.moleculesToRestrict = parse.listCast(currentUserInput.moleculesToRestrict) #Molecules range is a list    
    currentUserInput.moleculesToRestrict = parse.stripListOfStrings(currentUserInput.moleculesToRestrict)
    #if using signal range, then data lower/upper bound and increments needs to be the same length as the number of chosenMassFragments
    #if using concentration range, then they need to be the the same length as number of chosenMolecules
    if currentUserInput.signalOrConcentrationRange == 'signal': #So set lenOfParallelVectorizingBruteSolvingRestrictionVars to be the length of chosenMassFragments if using signal
        lenOfParallelVectorizingBruteSolvingRestrictionVars = len(chosenMassFragmentsForParsing)
    elif currentUserInput.signalOrConcentrationRange == 'concentration': #and set it equal to the length of chosenMolecules if using concentration
        lenOfParallelVectorizingBruteSolvingRestrictionVars = len(chosenMoleculesForParsing)
    #paralellVectorize data upper/lower bound and increments to the appropriate length
    currentUserInput.dataLowerBound = parse.parallelVectorize(currentUserInput.dataLowerBound,lenOfParallelVectorizingBruteSolvingRestrictionVars)
    currentUserInput.dataUpperBound = parse.parallelVectorize(currentUserInput.dataUpperBound,lenOfParallelVectorizingBruteSolvingRestrictionVars)
    currentUserInput.bruteIncrements = parse.parallelVectorize(currentUserInput.bruteIncrements,lenOfParallelVectorizingBruteSolvingRestrictionVars)
    
    #Set Scaling Factor
    parse.strCheck(currentUserInput.scaleRawDataOption,'scaleRawDataOption')
    if currentUserInput.scaleRawDataOption == 'manual':
        currentUserInput.scaleRawDataFactor = float(currentUserInput.scaleRawDataFactor) #scaleRawDataFactor is a float
    
    #Reference Correction Changer
    parse.strCheck(currentUserInput.tuningCorrection,'tuningCorrection')
    #The below two variables are no longer strings. They are now lists with two elements, each of which are strings. TODO: Change their names to referenceFileExistingTuningAndForm and referenceFileDesiredTuningAndForm
    #parse.strCheck(currentUserInput.referenceFileExistingTuningAndForm,'referenceFileExistingTuningAndForm')
    #parse.strCheck(currentUserInput.referenceFileDesiredTuningAndForm,'referenceFileDesiredTuningAndForm')
    
    #Reference Pattern Changer
    parse.strCheck(currentUserInput.extractReferencePatternFromDataOption,'extractReferencePatternFromDataOption')
    #If using reference pattern changer, check that all currentUserInput.rpcMoleculesToChange are in the referenceData
    if currentUserInput.extractReferencePatternFromDataOption == 'yes':
        #The molecules to change, their mass fragments, and time ranges are all lists
        currentUserInput.rpcMoleculesToChange = parse.listCast(currentUserInput.rpcMoleculesToChange)
        currentUserInput.rpcMoleculesToChange = parse.stripListOfStrings(currentUserInput.rpcMoleculesToChange)
        currentUserInput.rpcTimeRanges = parse.listCast(currentUserInput.rpcTimeRanges)
        currentUserInput.rpcTimeRanges = parse.parallelVectorize(currentUserInput.rpcTimeRanges,len(currentUserInput.rpcMoleculesToChange)) #rpcTimeRanges needs to have the same number of time ranges as moleculesToChange
        currentUserInput.rpcMoleculesToChangeMF = parse.listCast(currentUserInput.rpcMoleculesToChangeMF) #rpcMoleculesToChangeMF also needs to be of the same length but the mass fragments to change need to be hard coded in the user input so parallel vectorize is not feasible
        parse.compareElementsBetweenLists(currentUserInput.rpcMoleculesToChange,chosenMoleculesForParsing,'rpcMoleculesToChange','chosenMolecules')
    
    #Reference Mass Fragmentation Threshold
    parse.strCheck(currentUserInput.applyReferenceMassFragmentsThresholds,'applyReferenceMassFragmentsThresholds')
    if currentUserInput.applyReferenceMassFragmentsThresholds == 'yes': #If using reference mass fragmentation threshold
        currentUserInput.referenceMassFragmentFilterThreshold = parse.listCast(currentUserInput.referenceMassFragmentFilterThreshold) #reference value threshold is a list
        #The length of the reference value thresholds needs to be the same length as the number of molecules
        currentUserInput.referenceMassFragmentFilterThreshold = parse.parallelVectorize(currentUserInput.referenceMassFragmentFilterThreshold,len(chosenMoleculesForParsing))
                                                                                                                                                                                                                                                                            
        currentUserInput.referenceSignificantFragmentThresholds = parse.parallelVectorize(currentUserInput.referenceSignificantFragmentThresholds,len(chosenMoleculesForParsing))
    
    #Data Threshold Filter
    parse.strCheck(currentUserInput.lowerBoundThresholdChooser,'lowerBoundThresholdChooser')
    if currentUserInput.lowerBoundThresholdChooser == 'yes': #if using lowerBoundThresholdFilter
        #masstes to lower bound threshold filter and lower bound threshold percent/absolute are all three lists
        currentUserInput.massesToLowerBoundThresholdFilter = parse.listCast(currentUserInput.massesToLowerBoundThresholdFilter)
        currentUserInput.lowerBoundThresholdPercentage = parse.listCast(currentUserInput.lowerBoundThresholdPercentage)
        currentUserInput.lowerBoundThresholdAbsolute = parse.listCast(currentUserInput.lowerBoundThresholdAbsolute)        
        if len(currentUserInput.massesToLowerBoundThresholdFilter) == 0: #If currentUserInput.massesToLowerBoundThresholdFilter is empty
            currentUserInput.massesToLowerBoundThresholdFilter = chosenMassFragmentsForParsing #populate it with chosenMassFragments
        #if lowerBoundThresholdPercentage is empty, then user is option to use lowerBoundThresholdAbsolute
        if len(currentUserInput.lowerBoundThresholdPercentage) == 0:
            #and currentUserInput.lowerBoundThresholdAbsolute needs to be the same length as massesToLowerBoundThresholdFilter
            currentUserInput.lowerBoundThresholdAbsolute = parse.parallelVectorize(currentUserInput.lowerBoundThresholdAbsolute,len(currentUserInput.massesToLowerBoundThresholdFilter))
        elif len(currentUserInput.lowerBoundThresholdAbsolute) == 0: #Otherwise lowerBoundThresholdAbsolute is empty and the user has opted to use lowerBoundThresholdPercentage
            currentUserInput.lowerBoundThresholdPercentage = parse.parallelVectorize(currentUserInput.lowerBoundThresholdPercentage,len(currentUserInput.massesToLowerBoundThresholdFilter))
    
    #Data Smoother
    parse.strCheck(currentUserInput.dataSmootherYorN,'dataSmootherYorN')
    parse.strCheck(currentUserInput.dataSmootherChoice,'dataSmootherChoice')
    if currentUserInput.dataSmootherYorN == 'yes': #If using dataSmoother
        #The headers to confine to in data smoother is a list
        currentUserInput.dataSmootherHeadersToConfineTo = parse.listCast(currentUserInput.dataSmootherHeadersToConfineTo)        
        currentUserInput.dataSmootherHeadersToConfineTo = parse.stripListOfStrings(currentUserInput.dataSmootherHeadersToConfineTo)
        #mass fragments in headers to confine to must be included in chosenMassFragments
        parse.compareElementsBetweenLists(currentUserInput.dataSmootherHeadersToConfineTo,chosenMassFragmentsForParsing,'dataSmootherHeadersToConfineTo','chosenMolecules')
    
    #Raw Signal Threshold
    parse.strCheck(currentUserInput.applyRawSignalThresholds,'applyRawSignalThresholds')
    parse.strCheck(currentUserInput.rawSignalThresholdLimit,'rawSignalThresholdLimit')
    if currentUserInput.applyRawSignalThresholds == 'yes': #If using applyRawSignalThresholds
        #raw signal threshold value, sensitivity value, raw signal threshold divider, and raw signal threshold limit percent are all lists
        currentUserInput.rawSignalThresholdValue = parse.listCast(currentUserInput.rawSignalThresholdValue)
        currentUserInput.sensitivityThresholdValue = parse.listCast(currentUserInput.sensitivityThresholdValue)
        currentUserInput.rawSignalThresholdDivider = parse.listCast(currentUserInput.rawSignalThresholdDivider)
        currentUserInput.rawSignalThresholdLimitPercent = parse.listCast(currentUserInput.rawSignalThresholdLimitPercent)
        #sensitivityThreshold parallelVectorized to length of chosenMolecules
        #rawSignalThresholdValue, Divider, and LimitPercent all parallelVectorized to length of chosenMassFragments
        currentUserInput.rawSignalThresholdValue = parse.parallelVectorize(currentUserInput.rawSignalThresholdValue,len(chosenMassFragmentsForParsing))
    #   #TODO Commented out until bug in referenceThreshold is fixed    
    #    currentUserInput.sensitivityThresholdValue = parse.parallelVectorize(currentUserInput.sensitivityThresholdValue,len(chosenMoleculesForParsing))
        currentUserInput.rawSignalThresholdDivider = parse.parallelVectorize(currentUserInput.rawSignalThresholdDivider,len(chosenMassFragmentsForParsing))
        currentUserInput.rawSignalThresholdLimitPercent = parse.parallelVectorize(currentUserInput.rawSignalThresholdLimitPercent,len(chosenMassFragmentsForParsing))

    #Uncertainties
    
    
    #Negative Analyzer
    parse.strCheck(currentUserInput.negativeAnalyzerYorN,'negativeAnalyzerYorN')
    currentUserInput.NegativeAnalyzerTopNContributors = int(currentUserInput.NegativeAnalyzerTopNContributors)
    currentUserInput.NegativeAnalyzerBaseNumberOfGridIntervals = int(currentUserInput.NegativeAnalyzerBaseNumberOfGridIntervals)

    #Data Analysis Methods
    #All must be strings
    parse.strCheck(currentUserInput.solverChoice,'solverChoice')
    parse.strCheck(currentUserInput.uniqueOrCommon,'uniqueOrCommon')
    parse.strCheck(currentUserInput.slsFinish,'slsFinish')
    parse.strCheck(currentUserInput.objectiveFunctionType,'objectiveFunctionType')
    parse.strCheck(currentUserInput.distinguished,'distinguished')
    parse.strCheck(currentUserInput.fullBrute,'fullBrute')
    parse.strCheck(currentUserInput.SLSUniqueExport,'SLSUniqueExport')
    parse.strCheck(currentUserInput.finalOptimization,'finalOptimization')
        
    #Concentration Finder
    parse.strCheck(currentUserInput.concentrationFinder,'concentrationFinder')
    if currentUserInput.concentrationFinder == 'yes':
        #First cast the concentrationFinder variables as lists
        currentUserInput.moleculesTSC_List = parse.listCast(currentUserInput.moleculesTSC_List)
        currentUserInput.moleculesTSC_List = parse.stripListOfStrings(currentUserInput.moleculesTSC_List)
        currentUserInput.moleculeSignalTSC_List = parse.listCast(currentUserInput.moleculeSignalTSC_List)
        currentUserInput.massNumberTSC_List = parse.listCast(currentUserInput.massNumberTSC_List)
        currentUserInput.moleculeConcentrationTSC_List = parse.listCast(currentUserInput.moleculeConcentrationTSC_List)
        #Units needs to be a string, if it is not a string, return an error
        parse.strCheck(currentUserInput.unitsTSC,'unitsTSC')																																        
        
        if currentUserInput.TSC_List_Type == 'MultipleReferencePatterns': #If using multiple reference patterns then the user must input 1 value to use for each reference file or a value for each reference file
            #Then parallelize these variables to have the same length as number of reference patterns
            currentUserInput.moleculesTSC_List = parse.parallelVectorize(currentUserInput.moleculesTSC_List,len(currentUserInput.referencePatternsFileNamesList))
            currentUserInput.moleculeSignalTSC_List = parse.parallelVectorize(currentUserInput.moleculeSignalTSC_List,len(currentUserInput.referencePatternsFileNamesList))
            currentUserInput.massNumberTSC_List = parse.parallelVectorize(currentUserInput.massNumberTSC_List,len(currentUserInput.referencePatternsFileNamesList))
            currentUserInput.moleculeConcentrationTSC_List = parse.parallelVectorize(currentUserInput.moleculeConcentrationTSC_List,len(currentUserInput.referencePatternsFileNamesList))
            #NOTE: vectorizing these lists for 'SeparateMoleculesFactors' occurs in RatioFinder
            
    #Output Files
    #All must be strings
    parse.strCheck(currentUserInput.preProcessedDataOutputName,'preProcessedDataOutputName')
    parse.strCheck(currentUserInput.resolvedScaledConcentrationsOutputName,'resolvedScaledConcentrationsOutputName')
    parse.strCheck(currentUserInput.scaledConcentrationsPercentages,'scaledConcentrationsPercentages')
    parse.strCheck(currentUserInput.concentrationsOutputName,'concentrationsOutputName')
    parse.strCheck(currentUserInput.simulatedSignalsOutputName,'simulatedSignalsOutputName')
        
    #Iterative Analysis
    parse.strCheck(currentUserInput.TotalConcentrationsOutputName,'TotalConcentrationsOutputName')
        

    return None

    
def userInputValidityCheck(UserChoices): #Right now, currentUserInputModule is typically "G"
    #The incompatibilities dictionary is hardcoded.
    incompatibilitiesDict = {}
    incompatibilitiesDict['ReferencePatternChanger']=['ReferencePatternTimeChooser'] #These features are not compatible as of March 12th, 2019.  
    settingsCompatibilityCheck(UserChoices, incompatibilitiesDict)

    #The dependencies dictionary is hardcoded.
    dependenciesDict = {}
    #Note the form below: in right hand side, there are tuples. Foor the user's chocices, the value of the variable in index 0 must match the hardcoded value in index 1, otherwise the original feature has an incompatibility.
    dependenciesDict['SLSUniqueExport']={'yes':[(UserChoices['dataAnalysisMethods']['uniqueOrCommon'],'unique'),(UserChoices['dataAnalysisMethods']['solverChoice'],'sls')]}
    settingsDependenciesCheck(UserChoices, dependenciesDict)

    #Forcing of choices:
    if UserChoices['dataAnalysisMethods']['SLSUniqueExport'] == 'yes':
        if UserChoices['dataAnalysisMethods']['uniqueOrCommon'] != 'unique' or UserChoices['dataAnalysisMethods']['solverChoice'] != 'sls':
            UserChoices['dataAnalysisMethods']['SLSUniqueExport'] = 'no'
            print("Incompatible choice detected: forcing SLSUniqueExport to no.")
 
    if 'implicitSLScorrection' in UserChoices['dataAnalysisMethods']:
        if UserChoices['dataAnalysisMethods']['implicitSLScorrection'] == True:
            if UserChoices['dataAnalysisMethods']['uniqueOrCommon'] != 'unique' or UserChoices['dataAnalysisMethods']['solverChoice'] != 'sls':
                UserChoices['dataAnalysisMethods']['implicitSLScorrection'] = False
                print("Incompatible choice detected: forcing implicitSLScorrection to False.")
   
    #Filling settings variables dictionary so that variables can be populated from it. This is basically a mapping. See user input file for details.
    #The original variable names were single variables. Now, we are using a dictionary type structure (right side of equal signs) so they are being mapped to the single variables (left side of equal sign)
    #TODO: Consider if G.iterativeAnalysis = True or False should be changed to G.IterativeAnalysis_On or something like that, but will break backwards compatibility unless special care is taken.
    #Also to consider if other variables should change to have names like G.specificMolecules_chosenMoleculesNames. Probably not necessary since we have the dictionaries.
    SettingsVDictionary = {}  
    SettingsVDictionary['referencePatternsFileNamesList']   = UserChoices['inputFiles']['referencePatternsFileNamesList']
    SettingsVDictionary['referencePatternsFormsList']   = UserChoices['inputFiles']['referencePatternsFormsList']
    SettingsVDictionary['referencePatternTimeRanges']   = UserChoices['inputFiles']['referencePatternTimeRanges']
    SettingsVDictionary['dataToAnalyzeFileName']   = UserChoices['inputFiles']['dataToAnalyzeFileName']
    SettingsVDictionary['ionizationDataFileName']   = UserChoices['inputFiles']['ionizationDataFileName']
    
    SettingsVDictionary['preProcessing'] = UserChoices['preProcessing']['on'] 
    SettingsVDictionary['dataAnalysis'] = UserChoices['dataAnalysis']['on']
    SettingsVDictionary['dataSimulation'] = UserChoices['dataSimulation']['on']
    
    SettingsVDictionary['grapher'] = UserChoices['grapher']['on']
    SettingsVDictionary['stopAtGraphs'] = UserChoices['grapher']['stopAtGraphs']
    
    SettingsVDictionary['timeRangeLimit']= UserChoices['timeRangeLimit']['on']
    SettingsVDictionary['timeRangeStart'] = UserChoices['timeRangeLimit']['timeRangeStart']
    SettingsVDictionary['timeRangeFinish'] = UserChoices['timeRangeLimit']['timeRangeFinish']

    SettingsVDictionary['iterativeAnalysis'] = UserChoices['iterativeAnalysis']['on'] 
    SettingsVDictionary['TotalConcentrationsOutputName']    = UserChoices['iterativeAnalysis']['TotalConcentrationsOutputName']
    SettingsVDictionary['iterationSuffix']    = UserChoices['iterativeAnalysis']['iterationSuffix']
    SettingsVDictionary['unusedMolecules']    = UserChoices['iterativeAnalysis']['unusedMolecules']
    SettingsVDictionary['oldReferenceFileName']    = UserChoices['iterativeAnalysis']['oldReferenceFileName']
    SettingsVDictionary['oldDataToAnalyzeFileName']    = UserChoices['iterativeAnalysis']['oldDataToAnalyzeFileName']
    SettingsVDictionary['nextRefFileName']    = UserChoices['iterativeAnalysis']['nextRefFileName']
    SettingsVDictionary['nextExpFileName']    = UserChoices['iterativeAnalysis']['nextExpFileName']
    SettingsVDictionary['iterationNumber']    = UserChoices['iterativeAnalysis']['iterationNumber']  
     
    SettingsVDictionary['specificMolecules'] = UserChoices['specificMolecules']['on']
    SettingsVDictionary['chosenMoleculesNames'] = UserChoices['specificMolecules']['chosenMoleculesNames']
    SettingsVDictionary['specificMassFragments'] = UserChoices['specificMassFragments']['on']
    SettingsVDictionary['chosenMassFragments'] = UserChoices['specificMassFragments']['chosenMassFragments']
    
    SettingsVDictionary['moleculeLikelihoods'] = UserChoices['moleculeLikelihoods']['moleculeLikelihoods'] 
    
    SettingsVDictionary['sensitivityValues']  =  UserChoices['sensitivityValues']['sensitivityValues']
    
    SettingsVDictionary['linearBaselineCorrectionSemiAutomatic']    = UserChoices['linearBaselineCorrectionSemiAutomatic']['on']
    SettingsVDictionary['baselineType']    = UserChoices['linearBaselineCorrectionSemiAutomatic']['baselineType']
    SettingsVDictionary['massesToBackgroundCorrect']    = UserChoices['linearBaselineCorrectionSemiAutomatic']['massesToBackgroundCorrect']
    SettingsVDictionary['earlyBaselineTimes']    = UserChoices['linearBaselineCorrectionSemiAutomatic']['earlyBaselineTimes']
    SettingsVDictionary['lateBaselineTimes']    = UserChoices['linearBaselineCorrectionSemiAutomatic']['lateBaselineTimes']
    SettingsVDictionary['backgroundMassFragment']    = UserChoices['linearBaselineCorrectionManual']['backgroundMassFragment']
    SettingsVDictionary['backgroundSlopes']    = UserChoices['linearBaselineCorrectionManual']['backgroundSlopes']
    SettingsVDictionary['backgroundIntercepts']    = UserChoices['linearBaselineCorrectionManual']['backgroundIntercepts']
    
    SettingsVDictionary['interpolateYorN'] = UserChoices['interpolateYorN']['on'] 
    SettingsVDictionary['marginalChangeRestriction'] = UserChoices['interpolateYorN']['marginalChangeRestriction']
    SettingsVDictionary['ignorableDeltaYThreshold'] = UserChoices['interpolateYorN']['ignorableDeltaYThreshold']

    SettingsVDictionary['dataLowerBound']    = UserChoices['bruteSolvingRestrictions']['dataLowerBound']
    SettingsVDictionary['dataUpperBound']    = UserChoices['bruteSolvingRestrictions']['dataUpperBound']
    SettingsVDictionary['dataRangeSpecifierYorN']    = UserChoices['bruteSolvingRestrictions']['dataRangeSpecifierYorN']
    SettingsVDictionary['signalOrConcentrationRange']    = UserChoices['bruteSolvingRestrictions']['signalOrConcentrationRange']
    SettingsVDictionary['csvFile']    = UserChoices['bruteSolvingRestrictions']['csvFile']
    SettingsVDictionary['moleculesToRestrict']    = UserChoices['bruteSolvingRestrictions']['moleculesToRestrict']
    SettingsVDictionary['csvFileName']    = UserChoices['bruteSolvingRestrictions']['csvFileName']
    SettingsVDictionary['bruteIncrements']    = UserChoices['bruteSolvingRestrictions']['bruteIncrements']
    SettingsVDictionary['permutationNum']    = UserChoices['bruteSolvingRestrictions']['permutationNum']
    SettingsVDictionary['maxPermutations']    = UserChoices['bruteSolvingRestrictions']['maxPermutations']
  
    SettingsVDictionary['scaleRawDataYorN']   = UserChoices['scaleRawDataYorN']['on']
    SettingsVDictionary['scaleRawDataOption']   = UserChoices['scaleRawDataYorN']['scaleRawDataOption']
    SettingsVDictionary['scaleRawDataFactor']   = UserChoices['scaleRawDataYorN']['scaleRawDataFactor']

    SettingsVDictionary['tuningCorrection']    = UserChoices['tuningCorrection']['on']
    if 'createMixedTuningPattern' not in UserChoices['tuningCorrection']:
        UserChoices['tuningCorrection']['createMixedTuningPattern'] = True
    SettingsVDictionary['createMixedTuningPattern']  = UserChoices['tuningCorrection']['createMixedTuningPattern']
    SettingsVDictionary['referenceFileExistingTuningAndForm']    = UserChoices['tuningCorrection']['referenceFileExistingTuningAndForm']
    SettingsVDictionary['referenceFileDesiredTuningAndForm']    = UserChoices['tuningCorrection']['referenceFileDesiredTuningAndForm']
    SettingsVDictionary['referenceCorrectionCoefficients']    = UserChoices['tuningCorrection']['referenceCorrectionCoefficients']

    try: #to make sure old unit tests and analyses work.
        #if 'tuningCorrectorGasMixtureMoleculeNames' in UserChoices['tuningCorrection'].keys():
        SettingsVDictionary['tuningCorrectorGasMixtureMoleculeNames'] = UserChoices['tuningCorrection']['tuningCorrectorGasMixtureMoleculeNames']
    except: #to make sure old unit tests work.
        SettingsVDictionary['tuningCorrectorGasMixtureMoleculeNames'] = [] 
        UserChoices['tuningCorrection']['tuningCorrectorGasMixtureMoleculeNames'] = []
    try:
        SettingsVDictionary['referenceCorrectionCoefficients_cov']    = UserChoices['tuningCorrection']['referenceCorrectionCoefficients_cov']
    except:
        SettingsVDictionary['referenceCorrectionCoefficients_cov']    = [0,0,0] #TODO: This is to keep some old unit tests running. Ideally they should be fixed.
    SettingsVDictionary['extractReferencePatternFromDataOption']   = UserChoices['extractReferencePatternFromDataOption']['on']
    SettingsVDictionary['rpcMoleculesToChange']   = UserChoices['extractReferencePatternFromDataOption']['rpcMoleculesToChange']
    SettingsVDictionary['rpcTimeRanges']   = UserChoices['extractReferencePatternFromDataOption']['rpcTimeRanges']
    SettingsVDictionary['rpcMoleculesToChangeMF']    = UserChoices['extractReferencePatternFromDataOption']['rpcMoleculesToChangeMF'] 

    SettingsVDictionary['applyReferenceMassFragmentsThresholds']   = UserChoices['applyReferenceMassFragmentsThresholds']['on']
    SettingsVDictionary['referenceMassFragmentFilterThreshold']   = UserChoices['applyReferenceMassFragmentsThresholds']['referenceMassFragmentFilterThreshold']
    SettingsVDictionary['referenceSignificantFragmentThresholds']   = UserChoices['applyReferenceMassFragmentsThresholds']['referenceSignificantFragmentThresholds']
    
    SettingsVDictionary['lowerBoundThresholdChooser']   = UserChoices['lowerBoundThresholdChooser']['on'] 
    SettingsVDictionary['massesToLowerBoundThresholdFilter']   = UserChoices['lowerBoundThresholdChooser']['massesToLowerBoundThresholdFilter']
    SettingsVDictionary['lowerBoundThresholdPercentage']   = UserChoices['lowerBoundThresholdChooser']['lowerBoundThresholdPercentage']
    SettingsVDictionary['lowerBoundThresholdAbsolute']    = UserChoices['lowerBoundThresholdChooser']['lowerBoundThresholdAbsolute'] 

    SettingsVDictionary['dataSmootherYorN']   = UserChoices['dataSmootherYorN']['on']
    SettingsVDictionary['dataSmootherChoice']   = UserChoices['dataSmootherYorN']['dataSmootherChoice']
    SettingsVDictionary['dataSmootherTimeRadius']   = UserChoices['dataSmootherYorN']['dataSmootherTimeRadius']
    SettingsVDictionary['dataSmootherPointRadius']   = UserChoices['dataSmootherYorN']['dataSmootherPointRadius']
    SettingsVDictionary['dataSmootherHeadersToConfineTo']   = UserChoices['dataSmootherYorN']['dataSmootherHeadersToConfineTo']
    SettingsVDictionary['polynomialOrder']   = UserChoices['dataSmootherYorN']['polynomialOrder']

    SettingsVDictionary['applyRawSignalThresholds']   = UserChoices['applyRawSignalThresholds']['on']
    SettingsVDictionary['rawSignalThresholdValue']   = UserChoices['applyRawSignalThresholds']['rawSignalThresholdValue']
    SettingsVDictionary['sensitivityThresholdValue']   = UserChoices['applyRawSignalThresholds']['sensitivityThresholdValue']
    SettingsVDictionary['rawSignalThresholdDivider']   = UserChoices['applyRawSignalThresholds']['rawSignalThresholdDivider']
    SettingsVDictionary['rawSignalThresholdLimit']   = UserChoices['applyRawSignalThresholds']['rawSignalThresholdLimit']
    SettingsVDictionary['rawSignalThresholdLimitPercent']   = UserChoices['applyRawSignalThresholds']['rawSignalThresholdLimitPercent']
 
    SettingsVDictionary['calculateUncertaintiesInConcentrations'] 	=	    UserChoices['uncertainties']['calculateUncertaintiesInConcentrations'] 
    SettingsVDictionary['referencePatterns_uncertainties'] 	=	    UserChoices['uncertainties']['referencePatterns_uncertainties'] 
    SettingsVDictionary['dataToAnalyze_uncertainties']	=	    UserChoices['uncertainties']['dataToAnalyze_uncertainties']
    SettingsVDictionary['referenceCorrectionCoefficientsUncertainties'] 	=	    UserChoices['uncertainties']['referenceCorrectionCoefficientsUncertainties'] 
    SettingsVDictionary['referenceCorrectionCoefficientsIonizationUncertainties'] 	=	    UserChoices['uncertainties']['referenceCorrectionCoefficientsIonizationUncertainties'] 

 
    SettingsVDictionary['negativeAnalyzerYorN']   =UserChoices['negativeAnalyzerYorN']['on']
    if 'NegativeAnalyzerTopNContributors' in UserChoices['negativeAnalyzerYorN']:
        SettingsVDictionary['NegativeAnalyzerTopNContributors']   = UserChoices['negativeAnalyzerYorN']['NegativeAnalyzerTopNContributors']
    if 'NegativeAnalyzerBaseNumberOfGridIntervals' in UserChoices['negativeAnalyzerYorN']:
        SettingsVDictionary['NegativeAnalyzerBaseNumberOfGridIntervals']   = UserChoices['negativeAnalyzerYorN']['NegativeAnalyzerBaseNumberOfGridIntervals']
    
    SettingsVDictionary['solverChoice']   = UserChoices['dataAnalysisMethods']['solverChoice']
    SettingsVDictionary['uniqueOrCommon']   = UserChoices['dataAnalysisMethods']['uniqueOrCommon']
    SettingsVDictionary['slsWeighting']= UserChoices['dataAnalysisMethods']['slsWeighting']
    SettingsVDictionary['slsFinish']   = UserChoices['dataAnalysisMethods']['slsFinish']
    SettingsVDictionary['slsUniquePositiveConcentrationsOnly']   = UserChoices['dataAnalysisMethods']['slsUniquePositiveConcentrationsOnly']
    SettingsVDictionary['objectiveFunctionType']   = UserChoices['dataAnalysisMethods']['objectiveFunctionType']
    SettingsVDictionary['distinguished']   = UserChoices['dataAnalysisMethods']['distinguished']
    SettingsVDictionary['fullBrute']   = UserChoices['dataAnalysisMethods']['fullBrute']
    SettingsVDictionary['SLSUniqueExport']   = UserChoices['dataAnalysisMethods']['SLSUniqueExport']
    if 'implicitSLScorrection' in UserChoices['dataAnalysisMethods']:
        SettingsVDictionary['implicitSLScorrection'] = UserChoices['dataAnalysisMethods']['implicitSLScorrection']
    else:
        SettingsVDictionary['implicitSLScorrection'] = False #This is maintain backwards compatibility with old unit tests.
    SettingsVDictionary['finalOptimization']   = UserChoices['dataAnalysisMethods']['finalOptimization']

    SettingsVDictionary['concentrationFinder']   = UserChoices['concentrationFinder']['on']
    SettingsVDictionary['TSC_List_Type']   = UserChoices['concentrationFinder']['TSC_List_Type']
    SettingsVDictionary['moleculesTSC_List']   = UserChoices['concentrationFinder']['moleculesTSC_List']
    SettingsVDictionary['massNumberTSC_List']   = UserChoices['concentrationFinder']['massNumberTSC_List']
    SettingsVDictionary['moleculeSignalTSC_List']   = UserChoices['concentrationFinder']['moleculeSignalTSC_List']
    SettingsVDictionary['moleculeConcentrationTSC_List']   = UserChoices['concentrationFinder']['moleculeConcentrationTSC_List']
    SettingsVDictionary['unitsTSC']   = UserChoices['concentrationFinder']['unitsTSC']
    
    SettingsVDictionary['preProcessedDataOutputName']   = UserChoices['outputFiles']['preProcessedDataOutputName']
    SettingsVDictionary['resolvedScaledConcentrationsOutputName']   = UserChoices['outputFiles']['resolvedScaledConcentrationsOutputName']
    SettingsVDictionary['scaledConcentrationsPercentages']   = UserChoices['outputFiles']['scaledConcentrationsPercentages']
    SettingsVDictionary['concentrationsOutputName']   = UserChoices['outputFiles']['concentrationsOutputName']
    SettingsVDictionary['simulatedSignalsOutputName']   = UserChoices['outputFiles']['simulatedSignalsOutputName']
    
    SettingsVDictionary['ExportAtEachStep']= UserChoices['ExportAtEachStep']['on']
    
    SettingsVDictionary['generatePercentages'] = UserChoices['generatePercentages']['on']
    
    SettingsVDictionary['checkpoint']   = UserChoices['checkpoint']['checkpoint']
    SettingsVDictionary['start']   = UserChoices['checkpoint']['start'] 
    SettingsVDictionary['timeSinceLastCheckpoint']   = UserChoices['checkpoint']['timeSinceLastCheckpoint'] 

    return SettingsVDictionary
    
def settingsCompatibilityCheck(UserChoices, incompatibilitiesDict): #Right now, currentUserInputModule is typically "G"
    return None
    
def settingsDependenciesCheck(UserChoices, dependenciesDict, lastKey = None): #Right now, currentUserInputModule is typically "G"
    for key in dependenciesDict:
        if type(dependenciesDict[key])== type({}): #Checking for dictionary type. If it's a dictionary, need to call this function recursively.
            settingsDependenciesCheck(UserChoices, dependenciesDict[key], key)
        elif type(dependenciesDict[key])== type([]): #checking for list type, in which case it's a list of tuple pairs.
            dependenciesTuplesList = dependenciesDict[key]
            for tuplePair in dependenciesTuplesList:
                if tuplePair[0]==tuplePair[1]:
                    pass
                else:
                    print("Warning: Userinput with variable name |"+key+"| under |" + lastKey + "| is incompatible with the user choice of value |"+tuplePair[0] + "| used elsewhere in the userinput.") #Note: this warning needs to be better... may require changing the structure of dependenciesDict.
    return None

#This function populates module variable values from what's in a dictionary with 'key' becoming variable name, and value becoming value.
def populateModuleVariablesFromDictionary(moduleForAddingTo, inputDictionary): 
    for key in inputDictionary:
        setattr(moduleForAddingTo, key, inputDictionary[key])
    return None

#This function populates module variable values from what's in a dictionary with 'key' becoming variable name, and value becoming value. If any of the subitems are a dictionary, they are recursively called.
def populateModuleVariablesFromNestedDictionary(moduleForAddingTo, inputDictionary): 
    for key in inputDictionary:
        if (type(inputDictionary[key]) == type({})) and inputDictionary[key]: #this checks if it's a dictionary. If it is, then we need to call the function again (recursively). We have an exception for inputDictionary[key] being referenceCorrectionCoefficients because that one is supposed to be a dictionary.
            populateModuleVariablesFromNestedDictionary(moduleForAddingTo, inputDictionary[key]) 
        else:
            setattr(moduleForAddingTo, key, inputDictionary[key])
    return None