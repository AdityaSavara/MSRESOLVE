import os
import sys 
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
if os.path.basename(__file__) != "DefaultUserInput.py":
    from DefaultUserInput import *
thisModuleObject = sys.modules[__name__]
if str(__name__) == "DefaultUserInput":
    print("***This is in DefaultUserInput.py. DefaultUserInput gets loaded before UserInput to make sure all variables are initially populated with defaults.***")
    #Everything below this point is the variables' new dictionary format
    UserChoices = {} #initialize a dictionary to store all user choices
if str(__name__) != "DefaultUserInput":
    print("***This is in the actual UserInput.py.***")

#//Input Files//
UserChoices['inputFiles'] = {} #initialize the inputFiles container
UserChoices['inputFiles']['referenceFileNamesList'] = ['AcetaldehydeNISTRefMixed2.csv'] #enter the file name of the file containing reference information
UserChoices['inputFiles']['referenceFormsList'] = 'xyyy' #form is either 'xyyy' or 'xyxy' (if using reference pattern time chooser enter as list with forms for each individual reference file ['xyyy','xyyy','xyyy'])
UserChoices['inputFiles']['referencePatternTimeRanges'] = [] #Leave empty if not using reference pattern time chooser []
UserChoices['inputFiles']['collectedFileName'] = '2-CrotAcetExp#2.csv'	#enter the file name with raw mass spectrometer data

UserChoices['inputFiles']['ionizationDataFileName'] = '181017ProvidedIonizationData.csv' #the name of the file containing the ionization data

#Iterative Analysis
UserChoices['iterativeAnalysis'] = {} #initialize the iterativeAnalysis container
#Options are True, False, or '<name of iteration>'
UserChoices['iterativeAnalysis']['on'] = False
#the chosenMoleculesNames argument is used for iterative analysis, so make sure that input is accurate
#the chosenMassFragments argument is also used for iterative analysis, so make sure that input is accurate as well
#Below are Only used in iterative analysis. Most are just initializing and should not be touched.
UserChoices['iterativeAnalysis']['TotalConcentrationsOutputName'] = 'TotalConcentrations.csv'
UserChoices['iterativeAnalysis']['iterationSuffix'] = ''
UserChoices['iterativeAnalysis']['unusedMolecules'] =''
UserChoices['iterativeAnalysis']['oldReferenceFileName'] = []
UserChoices['iterativeAnalysis']['oldCollectedFileName'] ='' 
UserChoices['iterativeAnalysis']['nextRefFileName'] = []
UserChoices['iterativeAnalysis']['nextExpFileName'] = ''
UserChoices['iterativeAnalysis']['iterationNumber'] = None #just initializing.


#do you wish for the program to institute preproccessing and/or Data analysis?
#note that preproccesing must be done at least once before being bypassed 
#options for preProcessing and dataAnalysis are 'yes', 'load', or 'skip'
UserChoices['preProcessing'] = {} #initialize the preProcessing container
UserChoices['preProcessing']['on'] = 'yes'
UserChoices['dataAnalysis'] = {} #initialize the dataAnalysis container
UserChoices['dataAnalysis']['on'] = 'yes'
#options for dataSimulation are 'yes' or 'no' 
UserChoices['dataSimulation'] = {} #initialzie the dataSimulation container
UserChoices['dataSimulation']['on'] = 'yes'

#//Graphing//
UserChoices['grapher'] = {} #initialize grapher container
#option allowing you to view a graph of determined concentrations
UserChoices['grapher']['on'] = 'yes' #yes will graph function no will not
UserChoices['grapher']['stopAtGraphs'] = True #True will cause stopping at graphs.

#//Time Range//
UserChoices['timeRangeLimit'] = {} #initialize the timeRangeLimit container
#This function limits the data analyzed and proccessed to a certain subset of the total data
UserChoices['timeRangeLimit']['on'] = 'yes'	#if you wish to enable this function enter 'yes' otherwise 'no'
UserChoices['timeRangeLimit']['timeRangeStart'] = 176  #start time (-int)
UserChoices['timeRangeLimit']['timeRangeFinish'] = 900	#finish time (-int)

#//Chosen Molecules
UserChoices['specificMolecules'] = {} #initialize the specificMolecules container
#To choose only specific molecules to solve, input in a list of strings  below
UserChoices['specificMolecules']['on'] = 'no'
UserChoices['specificMolecules']['chosenMoleculesNames'] = ['Crotyl Alcohol']

#//Chosen Mass Fragments//
UserChoices['specificMassFragments'] = {} #initialize the specificMassFragments container
#To choose only specific mass fragments from collected data, input below:
UserChoices['specificMassFragments']['on'] = 'no'	#if you wish to enable this function enter 'yes' otherwise 'no'
UserChoices['specificMassFragments']['chosenMassFragments'] = [57] #enter the mass fragments you wish to include in calculations in the format [x,y,z...]

#//Molecule Likelihoods//
UserChoices['moleculeLikelihoods'] = {} #initialize the moleculeLikelihoods container
#To specify the percentage chance of detecting a particular molecule. This must be the same length as the number of molecules in the reference file, or have no value.
UserChoices['moleculeLikelihoods']['moleculeLikelihoods'] = [] #This should be like this [], or like this: [0.8, 1.0, 0.01,... 1.0] where the decimals are the user's guess of the likelihood of each molecule being present.
#//Sensivity Values//

UserChoices['sensitivityValues'] = {} #initialize the  container
#Sensitivity values allow the user the specify the threshold of each molecule individually, or apply one threshold to all molecules 
UserChoices['sensitivityValues']['sensitivityValues'] = []

#TODO 2/3/18: 
# Change so that late baseline times are omitted with a blank list for that mass (or all masses) rather than with zeros, 
# since this is not a good way of doing things.  Furthermore, after looking at the code, it does not even look like the code 
# is programmed to expect 0s, it looks like the code expects a blank list in order to skip the late baseline times.  
# I think maybe  the idea was that by putting [0,0] as the range, it cannot have any points in that range (not a finite range) 
# and therefore gets excluded from being fit.  Even if that's true, it needs to be checked that it's working.
#//Baseline Correction - Linear Semiautomatic//
UserChoices['linearBaselineCorrectionSemiAutomatic'] = {} #initialize the linearBaselineCorrectionSemiAutomatic container
#To change delete background average/slope based on time-interval, input below information below. 'yes' enables 'no' disables
UserChoices['linearBaselineCorrectionSemiAutomatic']['on'] = 'no'  #selection can be 'yes' or 'no'
UserChoices['linearBaselineCorrectionSemiAutomatic']['baselineType'] = ['linear'] 	#baselineType may be either 'flat' or 'linear'	
# if you would like to apply this correction to all fragments, leave as []
UserChoices['linearBaselineCorrectionSemiAutomatic']['massesToBackgroundCorrect'] = [2, 18, 27, 28, 31, 39, 41, 44, 57, 70]			#mflist: enter mass list delimited With commas [m1, m2, m3]
# to apply a uniform time range to all fragments, only insert one time range as such [[x,y]]
UserChoices['linearBaselineCorrectionSemiAutomatic']['earlyBaselineTimes'] = [[177.0, 177.1]]	# to apply different times for each fragment enter time pairs as such [[x,y],[z,w]..]
# to apply a uniform time range to all fragments, only insert one time range [[x,y]]
UserChoices['linearBaselineCorrectionSemiAutomatic']['lateBaselineTimes'] = [[901.0,902.0],[496.0,503.0],[901.0,902.0],[496.0,503.0],[901.0,902.0],[901.0,902.0],[901.0,902.0],[496.0,503.0],[901.0,902.0],[901.0,902.0]]	#if you do not wish to enter a second time list enter 0's [[0,0],[0,0]...] or [[0,0]]

#//Baseline Correction - Linear  Manual//
UserChoices['linearBaselineCorrectionManual'] = {} #initialize the linearBaselineCorrectionManual container
#To manually eliminate background slope/average, input below. If you do not wish to change any fragment leave function empty (like so [])
UserChoices['linearBaselineCorrectionManual']['backgroundMassFragment'] = []
UserChoices['linearBaselineCorrectionManual']['backgroundSlopes'] = []
UserChoices['linearBaselineCorrectionManual']['backgroundIntercepts'] = []

#//Data Solving Restrictions - Marginal Change Restrictor//
UserChoices['interpolateYorN'] = {} #initialize the interpolateYorN container
#To input data ranges for certain molecules, input below
UserChoices['interpolateYorN']['on'] = 'no'
#These factors are used to determine the search parameters for the brute force analysis method
# However, it is also needed for the interpolater which allows the data ranges to be small enough to be solvable by brute force
UserChoices['interpolateYorN']['marginalChangeRestriction'] = 2.0
UserChoices['interpolateYorN']['ignorableDeltaYThreshold'] = 0.01

#//Data Solving Restrictions - Brute Solving Restrictions 
UserChoices['bruteSolvingRestrictions'] = {} #initialize the bruteSolvingRestrictions container
#dataLowerBound and dataUpperbound put absolute bounds on the values brute searches for across all time points.
#I believe (if I am not mistaken) that dataRangeSpecifier puts bounds for each time point.
UserChoices['bruteSolvingRestrictions']['dataLowerBound'] = []
UserChoices['bruteSolvingRestrictions']['dataUpperBound'] = []
#NOTE: This feature (dataRangeSpecifierYorN) may not be compatible with the 
# "Load" pre-processing feature.
UserChoices['bruteSolvingRestrictions']['dataRangeSpecifierYorN'] = 'no' 
UserChoices['bruteSolvingRestrictions']['signalOrConcentrationRange'] = 'signal'	#'signal' or 'concentration'
UserChoices['bruteSolvingRestrictions']['csvFile'] = 'yes'	#'yes' or 'no'
UserChoices['bruteSolvingRestrictions']['moleculesToRestrict'] = []
UserChoices['bruteSolvingRestrictions']['csvFileName'] = 'rangestemplate.csv'
#NOTE: The increment choice of the user is then possibly overridden based on 
# the values of maxPermutations (the number of molecules and bruteIncrements might 
# cause too large of a number of permutations, in which case larger bruteIncrements 
# may be used).
#bruteIncrements sets the size of the increments for Brute (e.g., if we said  0.01 bar, it would make the 
# separation between points 0.01 bar in the grid, for that axis). 
UserChoices['bruteSolvingRestrictions']['bruteIncrements'] = []
UserChoices['bruteSolvingRestrictions']['permutationNum'] = 1000
UserChoices['bruteSolvingRestrictions']['maxPermutations'] = 100001

#// Set Scaling Factor?
UserChoices['scaleRawDataYorN'] = {} #initialize the scaleRawDataYorN container
UserChoices['scaleRawDataYorN']['on'] = 'no' #This variable is currently unused, but later choosing "no" will set things to "manual" and change "scaleRawDataFactor" to 1.
UserChoices['scaleRawDataYorN']['scaleRawDataOption'] = 'manual' #Choices are 'manual' or 'auto'
#'auto' automatically scales the data so that the smallest value is equal to 1
#If manual is chosen above, this option allows the user to scale all raw data by a factor of their choosing 
UserChoices['scaleRawDataYorN']['scaleRawDataFactor'] = 1
#Note that 1 is the default and will make no alteration to the data

#//Tuning Corrector - Reference Correction Coefficients//
UserChoices['measuredReferenceYorN'] = {} #initialize the measuredReferenceYorN container
#TODO Reference Correction Coefficients feature should be upgraded to enable separate coefficients for each molecule to allow mixing and matching of reference patterns
#TODO This can be tested by looking at the exported reference file and comparing it to the existing reference file
#To change reference data based on mass dependent 2nd degree polynomial fit, input polynomial below. If you do not wish to use this function, simply leave as default
UserChoices['measuredReferenceYorN']['on'] ='no'
UserChoices['measuredReferenceYorN']['referenceFileExistingTuningAndForm'] =['ReferenceCollected.csv','xyyy']
UserChoices['measuredReferenceYorN']['referenceFileDesiredTuningAndForm'] =['ReferenceLiterature.csv','xyyy']
#The reference correction coefficients are always used.  If measuredReferenceYorN is 'yes' then the coefficients are overwritten (Ashi thinks it will generate a new reference pattern)
UserChoices['measuredReferenceYorN']['referenceCorrectionCoefficients'] = {'A': 0.0, 'B': 0.0, 'C': 1.0}	
                            #default is 'A': 0.0, 'B': 0.0, 'C': 1.0.   Used as.... Factor = A*X^2 + B*X + C, so A=0,B=0,C=1.0 means the final factor is 1.0 and independent of molecular weight.
UserChoices['measuredReferenceYorN']['referenceCorrectionCoefficients_cov'] = [0,0,0] #Covariance for reference correction coefficients for tuning corrector. Default is 0,0,0. Can be 9 x 9 covariance.


#//Reference Pattern Changer // (rpc)
UserChoices['extractReferencePatternFromDataOption'] = {} #initialize the extractReferencePatternFromDataOption container
#To change reference data based on collected data at a certain time, enter mass fragments for the molecule and times below
UserChoices['extractReferencePatternFromDataOption']['on'] = 'no'
UserChoices['extractReferencePatternFromDataOption']['rpcMoleculesToChange'] = ['Crotyl Alcohol']
#rpcTimeRanges and rpcMoleculesToChangeMF are nested lists.  Each nested list corresponds to a molecule in rpcMoleculesToChange
#To make this easier to visualize, each nested list is placed on its own line so the first line refers to the first molecule, second line refers to the second molecule and so on
UserChoices['extractReferencePatternFromDataOption']['rpcTimeRanges'] = [
                                                                         [300,500], #For each molecule to be changed, a pair of times is required.
                                                                         ]
#The first mass fragment is the base fragment and it will not be changed.  The fragments following the first one are all altered based on the signal of the first fragment from the collected data
UserChoices['extractReferencePatternFromDataOption']['rpcMoleculesToChangeMF'] = [
                                                                                  [70,57], #For each molecule for using the rpc on, make a new line with a list of masses (length of each should be greater than 1).
                                                                                  ]

#//Reference Mass Fragmentation Threshold//
UserChoices['minimalReferenceValue'] = {} #initialize the minimalReferenceValue container
# if you want to exclude tiny fragmentation peaks
UserChoices['minimalReferenceValue']['on'] = 'no'
UserChoices['minimalReferenceValue']['referenceValueThreshold'] = [2.0]
UserChoices['minimalReferenceValue']['referenceSignificantFragmentThresholds'] = [6.0]

#//Data Threshold Filter//
UserChoices['lowerBoundThresholdChooser'] = {} #initialize the lowerBoundThresholdChooser container
#To change the lower bound below which data is eliminated, change below; lowerBoundThresholdChooser ='yes' or 'no'
#The idea is that below an absolute (or percent based) threshold the intensity will be set to 0.
UserChoices['lowerBoundThresholdChooser']['on'] = 'no' 
UserChoices['lowerBoundThresholdChooser']['massesToLowerBoundThresholdFilter'] = [] # leave as [ ] to apply identical percentages or absolute thresholds to all masses.
UserChoices['lowerBoundThresholdChooser']['lowerBoundThresholdPercentage'] = [0.25]  # 1.0 is max value. leave as [ ] to only use the absolute threshold. Always include a decimal. 
UserChoices['lowerBoundThresholdChooser']['lowerBoundThresholdAbsolute'] = []  # leave as [ ] to only use the percentage threshold. Always include a decimal.

#TODO change the name option from point/timerange to 
# abscissaPointradius / abscissaDistanceRadius
#//Data Smoothing//
UserChoices['dataSmootherYorN'] = {} #initialize the dataSmootherYorN container
#This section is for the data smoother function which, by default, is enabled. 
#Data smoothing can be conducted by a time basis or by a data point basis
UserChoices['dataSmootherYorN']['on'] = 'yes'
UserChoices['dataSmootherYorN']['dataSmootherChoice'] = 'timerange'	#options are 'pointrange' or 'timerange'
# abscissaPointRadius and absc
UserChoices['dataSmootherYorN']['dataSmootherTimeRadius'] = 7
UserChoices['dataSmootherYorN']['dataSmootherPointRadius'] = 5
UserChoices['dataSmootherYorN']['dataSmootherHeadersToConfineTo'] = [] #Masses on which to perform smoothing.
                           #Should be a subset of 'choosenMassFragments' above.
                           # leave array empty [], to smooth all masses
UserChoices['dataSmootherYorN']['polynomialOrder'] = 1  #During the local smoothing, a linear fit (or polynomial fit) is applied.

#//Raw Signal Threshold//
UserChoices['rawSignalThresholdMethod'] = {} #initialize the rawSignalThresholdMethod container
#To change the threshold at which raw signals are not longer relevant, change below (similar to above function, but for rows instead of columns)
#We think the reference to the 'above function' in the previous line is referring to Data Threshold Filter
#These signals get converted into 0.
#WARNING: This function is highly complex and should be considered a work in progress. It cannot be confirmed to work properly (as of 7/18/17).
UserChoices['rawSignalThresholdMethod']['on'] = 'no'
UserChoices['rawSignalThresholdMethod']['rawSignalThresholdValue'] = [.0000001]
UserChoices['rawSignalThresholdMethod']['sensitivityThresholdValue'] = [1] #this is the number in the Reference given the relative intensity of the signal of the mass fragment
UserChoices['rawSignalThresholdMethod']['rawSignalThresholdDivider'] = []
#Part of previous entry function, but this function enables the user to change the sum of raw signals, allowing molecules with very high concentrations not to affect previous funciton
UserChoices['rawSignalThresholdMethod']['rawSignalThresholdLimit'] = 'no'
UserChoices['rawSignalThresholdMethod']['rawSignalThresholdLimitPercent'] = []

#//Uncertainties for Calculating Uncertainties in Concentrations//
UserChoices['uncertainties'] = {}
UserChoices['uncertainties']['calculateUncertaintiesInConcentrations'] = False
UserChoices['uncertainties']['referenceFileUncertainties'] = 2 #which can be a float/integer for absolute uncertainties or the value True (or the value 'File'. Will expect same file name as reference file with _uncertainties after that) or the value None (False will also be set to None) . For example, the value 2 would mean a 2% uncertainty for the value 100, but a 50% uncertainty for the value of 4.
UserChoices['uncertainties']['collectedFileUncertainties'] = 'Auto' # or 'File' or 'None'. Will expect same file name as collected file with _uncertainties after that).
UserChoices['uncertainties']['referenceCorrectionCoefficientsUncertainties'] = None #Else a dictionary of uncertainties for 'A', 'B', 'C'. Not yet implemented.
UserChoices['uncertainties']['referenceCorrectionCoefficientsIonizationUncertainties'] = None #Not yet implemented.

#//Negative Analyzer//
UserChoices['negativeAnalyzerYorN'] = {} #initialize the negativeAnalyzerYorN container
#if enabled ('yes') Negative Analyzer will prevernt negative valued concentrations from being compututed.  
UserChoices['negativeAnalyzerYorN']['on'] = 'no'
UserChoices['negativeAnalyzerYorN']['NegativeAnalyzerTopNContributors'] = 5
UserChoices['negativeAnalyzerYorN']['NegativeAnalyzerBaseNumberOfGridIntervals'] = 5

#//Data Analysis Methods
UserChoices['dataAnalysisMethods'] = {} #initialize the dataAnalysisMethods container
#Below the path for the analysis of the data; sls or inverse
UserChoices['dataAnalysisMethods']['answer'] = 'inverse'	#'inverse' or 'sls'; sls is suggested
UserChoices['dataAnalysisMethods']['uniqueOrCommon'] = 'common'	#'unique' or 'common'; common is suggested
UserChoices['dataAnalysisMethods']['slsWeighting'] = [1,0,0,0] #The first uses uncertainties weighting. The second solves for largest concentrations first. The third uses reference peak height. The fourth uses the signal intensity.  All can be on at the same time. 
UserChoices['dataAnalysisMethods']['slsFinish'] = 'brute'	#'brute' or 'inverse'; brute is suggested
UserChoices['dataAnalysisMethods']['slsUniquePositiveConcentrationsOnly'] = False #Can be true or false. This is faster but less accurate than NegativeAnalyzer
UserChoices['dataAnalysisMethods']['objectiveFunctionType'] = 'ssr'	#objectiveFunctionType = 'ssr', 'sar', 'weightedSAR' or 'weightedSSR' 
UserChoices['dataAnalysisMethods']['distinguished'] = 'yes'
UserChoices['dataAnalysisMethods']['fullBrute'] = 'yes'
UserChoices['dataAnalysisMethods']['SLSUniqueExport'] = 'yes'
UserChoices['dataAnalysisMethods']['implicitSLScorrection'] = False
UserChoices['dataAnalysisMethods']['finalOptimization'] = 'None' #options are 'None' or... 'Nelder-Mead','Powell','CG','BFGS','Newton-CG','L-BFGS-B','TNC','COBYLA','SLSQP','dogleg','trust-ncg','trust-xact','trust-krylov' from scipy.optimize.minimze

#//Concentration Finder//
UserChoices['concentrationFinder'] = {} #initialize the concentrationFinder container
#this last set of inputs is where you enter your conversion factors from raw signal to concentration, unlike most rows, do not leave brackets around chosen numbers
#here you put in a known raw signal intensity and the known concentration it corresponds to. 
#TODO Note: Concentration Finder is not compatible with simultaneous use of multiple reference files AND separate molecules' factors as of 181022. Currently, it a user may use either one or the other.
UserChoices['concentrationFinder']['on'] = 'no'
UserChoices['concentrationFinder']['TSC_List_Type'] = 'MultipleReferencePatterns' #Options are 'MultipleReferencePatterns' or 'SeparateMoleculesFactors'
UserChoices['concentrationFinder']['moleculesTSC_List'] = 'Acetaldehyde'
UserChoices['concentrationFinder']['massNumberTSC_List'] = 29
UserChoices['concentrationFinder']['moleculeSignalTSC_List'] = 1.66945 #This is the list of intensity values that correspond to a known concentration to scale with for the same mass fragments
UserChoices['concentrationFinder']['moleculeConcentrationTSC_List'] = 0.05	#this is the concentration/pressure associated with the signal.
UserChoices['concentrationFinder']['unitsTSC'] = 'bar'	#this string is the unit for the concentration. The unit will not be used in calculations so any units may be used

#//Output Files//
UserChoices['outputFiles'] = {} #initialize the outputFiles container
#the last section designates the various output files, all are suppose to be csv values
#If files do not exist they will be generated
UserChoices['outputFiles']['preProcessedDataOutputName'] = 'PreProcessedData.csv'
UserChoices['outputFiles']['resolvedScaledConcentrationsOutputName'] = 'ScaledConcentrations.csv'
UserChoices['outputFiles']['scaledConcentrationsPercentages'] = 'ScaledConcentrationPercentages.csv'
UserChoices['outputFiles']['concentrationsOutputName'] = 'ResolvedConcentrations.csv'
UserChoices['outputFiles']['simulatedSignalsOutputName'] = 'SimulatedRawSignals.csv'

UserChoices['ExportAtEachStep'] = {} #initialize the ExportAtEachStep container
UserChoices['ExportAtEachStep']['on'] = 'yes'
UserChoices['generatePercentages'] = {} #initialize the generatePercentages container
UserChoices['generatePercentages']['on'] = 'no'

UserChoices['checkpoint'] = {} #initialize the  container
UserChoices['checkpoint']['checkpoint'] = '' #just initializing, not really necessary to have here.
UserChoices['checkpoint']['start'] = '' #just initializing, not really necessary to have here.
UserChoices['checkpoint']['timeSinceLastCheckpoint'] = '' #just initializing, not really necessary to have here.



from userInputValidityFunctions import userInputValidityCheck
from userInputValidityFunctions import populateModuleVariablesFromDictionary
SettingsVDictionary = userInputValidityCheck(UserChoices)
populateModuleVariablesFromDictionary(thisModuleObject, SettingsVDictionary)
####End of temporary code####
__var_list__ = ['referenceFileNamesList','referenceFormsList','collectedFileName','referencePatternTimeRanges','ionizationDataFileName','iterativeAnalysis','iterationNumber','iterationSuffix','unusedMolecules','oldReferenceFileName', 'oldCollectedFileName', 'nextRefFileName', 'nextExpFileName','preProcessing','dataAnalysis','dataSimulation','grapher','stopAtGraphs','timeRangeLimit','timeRangeStart','timeRangeFinish','specificMolecules','chosenMoleculesNames','specificMassFragments','chosenMassFragments','moleculeLikelihoods','sensitivityValues','linearBaselineCorrectionSemiAutomatic','baselineType','massesToBackgroundCorrect','earlyBaselineTimes','lateBaselineTimes','backgroundMassFragment','backgroundSlopes','backgroundIntercepts','interpolateYorN','marginalChangeRestriction','ignorableDeltaYThreshold','dataLowerBound','dataUpperBound','dataRangeSpecifierYorN','signalOrConcentrationRange','csvFile','moleculesToRestrict','csvFileName','bruteIncrements','permutationNum','maxPermutations','scaleRawDataOption','scaleRawDataFactor','measuredReferenceYorN','referenceFileExistingTuningAndForm','referenceFileDesiredTuningAndForm','referenceCorrectionCoefficients','referenceCorrectionCoefficients_cov','extractReferencePatternFromDataOption','rpcMoleculesToChange','rpcMoleculesToChangeMF','rpcTimeRanges','minimalReferenceValue','referenceValueThreshold','referenceSignificantFragmentThresholds','lowerBoundThresholdChooser','massesToLowerBoundThresholdFilter','lowerBoundThresholdPercentage','lowerBoundThresholdAbsolute','dataSmootherYorN','dataSmootherChoice','dataSmootherTimeRadius','dataSmootherPointRadius','dataSmootherHeadersToConfineTo','polynomialOrder','rawSignalThresholdMethod','rawSignalThresholdValue','sensitivityThresholdValue','rawSignalThresholdDivider','rawSignalThresholdLimit','rawSignalThresholdLimitPercent','negativeAnalyzerYorN','NegativeAnalyzerTopNContributors','NegativeAnalyzerBaseNumberOfGridIntervals','calculateUncertaintiesInConcentrations' , 'referenceFileUncertainties' ,'collectedFileUncertainties','referenceCorrectionCoefficientsUncertainties', 'referenceCorrectionCoefficientsIonizationUncertainties' ,'answer','uniqueOrCommon','slsWeighting','slsFinish','slsUniquePositiveConcentrationsOnly','objectiveFunctionType','distinguished','fullBrute','SLSUniqueExport', 'implicitSLScorrection', 'finalOptimization', 'concentrationFinder','moleculesTSC_List','TSC_List_Type','moleculeSignalTSC_List','massNumberTSC_List','moleculeConcentrationTSC_List','unitsTSC','preProcessedDataOutputName','resolvedScaledConcentrationsOutputName','scaledConcentrationsPercentages','concentrationsOutputName','simulatedSignalsOutputName','TotalConcentrationsOutputName','ExportAtEachStep','generatePercentages','checkpoint','start','timeSinceLastCheckpoint', 'iterationNumber']