import os
import sys 
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
if os.path.basename(__file__) != "DefaultUserInput.py":
    from DefaultUserInput import *

#USER INPUT FILE
#//Input Files//
referencePatternsFileNamesList = ['AcetaldehydeNISTRefMixed2.tsv','AcetaldehydeNISTRefMixed2Edit.csv'] #enter the file name of the file containing reference information
referencePatternsFormsList = 'xyyy'	#form is either 'xyyy' or 'xyxy'
referencePatternTimeRanges = [[176,200],[250,400]] #Leave empty if not using reference pattern time chooser []
dataToAnalyzeFileName = '2-CrotAcetExp#2Truncated.csv'	#enter the file name with raw mass spectrometer data

ionizationDataFileName = '_ProvidedIonizationData.csv'

#Iterative Analysis
#Options are True, False, or '<name of iteration>'
iterativeAnalysis = False
#the chosenMolecules argument is used for iterative analysis, so make sure that input is accurate
#the chosenMassFragments argument is also used for iterative analysis, so make sure that input is accurate as well
iterationNumber = None

#do you wish for the program to institute preproccessing and/or Data analysis?
#note that preproccesing must be done at least once before being bypassed 
#options for preProcessing and dataAnalysis are 'yes', 'load', or 'skip'
preProcessing= 'yes'
dataAnalysis='yes'
#options for dataSimulation are 'yes' or 'no' 
dataSimulation='yes'

#//Graphing//
#option allowing you to view a graph of determined concentrations
grapher = 'no' #yes will graph function no will not


#//Time Range//
#This function limits the data analyzed and proccessed to a certain subset of the total data
timeRangeLimit = 'yes'	#if you wish to enable this function enter 'yes' otherwise 'no'
timeRangeStart = 176  #start time (-int)
timeRangeFinish = 200	#finish time (-int)

#//Chosen Molecules
#To choose only specific molecules to solve, input in a list of strings  below
specificMolecules = 'no'
chosenMoleculesNames = ['CO', 'CO2', 'H2O']

#//Chosen Mass Fragments//
#To choose only specific mass fragments from collected data, input below:
specificMassFragments = 'yes'	#if you wish to enable this function enter 'yes' otherwise 'no'
chosenMassFragments = [2, 18, 27, 28, 31, 39, 41, 44, 57, 70] #enter the mass fragments you wish to include in calculations in the format [x,y,z...]

#//Molecule Likelihoods//
#To specify the percentage chance of detecting a particular molecule. This must be the same length as the number of molecules in the reference file, or have no value.
moleculeLikelihoods = [] #This should be like this [], or like this: [0.8, 1.0, 0.01,... 1.0] where the decimals are the user's guess of the likelihood of each molecule being present.
#//Sensivity Values//
#Sensitivity values allow the user the specify the threshold of each molecule individually, or apply one threshold to all molecules 
sensitivityValues = []

#TODO 2/3/18: 
# Change so that late baseline times are omitted with a blank list for that mass (or all masses) rather than with zeros, 
# since this is not a good way of doing things.  Furthermore, after looking at the code, it does not even look like the code 
# is programmed to expect 0s, it looks like the code expects a blank list in order to skip the late baseline times.  
# I think maybe  the idea was that by putting [0,0] as the range, it cannot have any points in that range (not a finite range) 
# and therefore gets excluded from being fit.  Even if that's true, it needs to be checked that it's working.
#//Baseline Correction - Linear Semiautomatic//
#To change delete background average/slope based on time-interval, input below information below. 'yes' enables 'no' disables
linearBaselineCorrectionSemiAutomatic = 'no'  #selection can be 'yes' or 'no'
baselineType = ['linear'] 	#baselineType may be either 'flat' or 'linear'	
# if you would like to apply this correction to all fragments, leave as []
massesToBackgroundCorrect = [2, 18, 27, 28, 31, 39, 41, 44, 57, 70]			#mflist: enter mass list delimited With commas [m1, m2, m3]
# to apply a uniform time range to all fragments, only insert one time range as such [[x,y]]
earlyBaselineTimes = [[177.0, 177.1]]	# to apply different times for each fragment enter time pairs as such [[x,y],[z,w]..]
# to apply a uniform time range to all fragments, only insert one time range [[x,y]]
lateBaselineTimes = [[901.0,902.0],[496.0,503.0],[901.0,902.0],[496.0,503.0],[901.0,902.0],[901.0,902.0],[901.0,902.0],[496.0,503.0],[901.0,902.0],[901.0,902.0]]	#if you do not wish to enter a second time list enter 0's [[0,0],[0,0]...] or [[0,0]]


#//Baseline Correction - Linear  Manual//
#To manually eliminate background slope/average, input below. If you do not wish to change any fragment leave function empty (like so [])
backgroundMassFragment = []
backgroundSlopes = []
backgroundIntercepts = []

#//Data Solving Restrictions - Marginal Change Restrictor//
#To input data ranges for certain molecules, input below
interpolateYorN = 'no'
#These factors are used to determine the search parameters for the brute force analysis method
# However, it is also needed for the interpolater which allows the data ranges to be small enough to be solvable by brute force
marginalChangeRestriction = 2.0
ignorableDeltaYThreshold = 0.01
dataToAnalyze_uncertainties_radiusType = 'pointrange'
referenceCorrectionCoefficientsUncertainties = None
referenceCorrectionCoefficientsIonizationUncertainties = None
#//Data Solving Restrictions - Brute Solving Restrictions 
#dataLowerBound and dataUpperbound put absolute bounds on the values brute searches for across all time points.
#I believe (if I am not mistaken) that dataRangeSpecifier puts bounds for each time point.
dataLowerBound = []
dataUpperBound = []
#NOTE: This feature (dataRangeSpecifierYorN) may not be compatible with the 
# "Load" pre-processing feature.
dataRangeSpecifierYorN = 'no' 
signalOrConcentrationRange = 'signal'	#'signal' or 'concentration'
csvFile = 'yes'	#'yes' or 'no'
moleculesToRestrict = []
csvFileName = 'rangestemplate.csv'
#NOTE: The increment choice of the user is then possibly overridden based on 
# the values of maxPermutations (the number of molecules and bruteIncrements might 
# cause too large of a number of permutations, in which case larger bruteIncrements 
# may be used).
#bruteIncrements sets the size of the increments for Brute (e.g., if we said  0.01 bar, it would make the 
# separation between points 0.01 bar in the grid, for that axis). 
bruteIncrements = []
permutationNum = 1000
maxPermutations = 100001



#// Set Scaling Factor?
scaleRawDataOption = 'manual' #Choices are 'manual' or 'auto'
#'auto' automatically scales the data so that the smallest value is equal to 1
#If manual is chosen above, this option allows the user to scale all raw data by a factor of their choosing 
scaleRawDataFactor = 1
#Note that 1 is the default and will make no alteration to the data

#//Reference Correction Coefficients//
#To change reference data based on mass dependent 2nd degree polynomial fit, input polynomial below. If you do not wish to use this function, simply leave as default
tuningCorrection='no'
referenceFileExistingTuningAndForm='AcetaldehydeMeasuredRef.csv'
referenceFileDesiredTuningAndForm ='AcetaldehydeOnlyNISTRef.csv'
referenceCorrectionCoefficients = {'A': 0.0, 'B': 0.0, 'C': 1.0}	
                            #default is 'A': 0.0, 'B': 0.0, 'C': 1.0
tuningCorrectorGasMixtureConcentrations = []
tuningCorrectorGasMixtureSignals = []

#//Reference Pattern Changer // (rpc)
#To change reference data based on collected data at a certain time, enter mass fragments for the molecule and times below
extractReferencePatternFromDataOption = 'no'
rpcMoleculesToChange = ['Crotyl Alcohol']
#rpcTimeRanges and rpcMoleculesToChangeMF are nested lists.  Each nested list corresponds to a molecule in rpcMoleculesToChange
#To make this easier to visualize, each nested list is placed on its own line so the first line refers to the first molecule, second line refers to the second molecule and so on
rpcTimeRanges = [
                 [300,500], #For each molecule to be changed, a pair of times is required.
                ]
#The first mass fragment is the base fragment and it will not be changed.  The fragments following the first one are all altered based on the signal of the first fragment from the collected data
rpcMoleculesToChangeMF = [
                          [70,57], #For each molecule for using the rpc on, make a new line with a list of masses (length of each should be greater than 1).
                         ]



#//Reference Mass Fragmentation Threshold//
# if you want to exclude tiny fragmentation peaks
applyReferenceMassFragmentsThresholds= 'no'
referenceMassFragmentFilterThreshold = [6.0]


#//Data Threshold Filter//
#To change the lower bound below which data is eliminated, change below; lowerBoundThresholdChooser ='yes' or 'no'
#The idea is that below an absolute (or percent based) threshold the intensity will be set to 0.
lowerBoundThresholdChooser = 'no' 
massesToLowerBoundThresholdFilter = [] # leave as [ ] to apply identical percentages or absolute thresholds to all masses.
lowerBoundThresholdPercentage = [0.25]  # 1.0 is max value. leave as [ ] to only use the absolute threshold. Always include a decimal. 
lowerBoundThresholdAbsolute = []  # leave as [ ] to only use the percentage threshold. Always include a decimal.


#TODO change the name option from point/timerange to 
# abscissaPointradius / abscissaDistanceRadius
#//Data Smoothing//
#This section is for the data smoother function which, by default, is enabled. 
#Data smoothing can be conducted by a time basis or by a data point basis
dataSmootherYorN = 'yes'
dataSmootherChoice = 'timerange'	#options are 'pointrange' or 'timerange'
# abscissaPointRadius and absc
dataSmootherTimeRadius = 7
dataSmootherPointRadius = 5
dataSmootherHeadersToConfineTo = [] #Masses on which to perform smoothing.
                           #Should be a subset of 'choosenMassFragments' above.
                           # leave array empty [], to smooth all masses
polynomialOrder = 1  #During the local smoothing, a linear fit (or polynomial fit) is applied.


#//Raw Signal Threshold//
#To change the threshold at which raw signals are not longer relevant, change below (similar to above function, but for rows instead of columns)
#These signals get converted into 0.
#WARNING: This function is highly complex and should be considered a work in progress. It cannot be confirmed to work properly (as of 7/18/17).
applyRawSignalThresholds = 'no'
rawSignalThresholdValue = [.0000001]
sensitivityThresholdValue = [1] #this is the number in the Reference given, the relative intensity of the signal of the mass fragment
rawSignalThresholdDivider = []
#Part of previous entry function, but this function enables the user to change the sum of raw signals, allowing molecules with very high concentrations not to affect previous funciton
rawSignalThresholdLimit = 'no'
rawSignalThresholdLimitPercent = []


#//Negative Analyzer//
#if enabled ('yes') Negative Analyzer will prevernt negative valued concentrations from being compututed.  
negativeAnalyzerYorN = 'no'


#//Data Analysis Methods
#Below the path for the analysis of the data; sls or inverse
solverChoice = 'sls'	#'inverse' or 'sls'; sls is suggested
uniqueOrCommon = 'common'	#'unique' or 'common'; common is suggested
slsWeighting = [0,0,1,0]
slsFinish = 'brute'	#'brute' or 'inverse'; brute is suggested
objectiveFunctionType = 'ssr'	#objectiveFunctionType = 'ssr', 'sar', 'weightedSAR' or 'weightedSSR' 
distinguished = 'yes'
fullBrute = 'yes'
SLSUniqueExport = 'yes'


#//Concentration Finder//
#this last set of inputs is where you enter your conversion factors from raw signal to concentration, unlike most rows, do not leave brackets around chosen numbers
#here you put in a known raw signal intensity and the known concentration it corresponds to. 
concentrationFinder = 'no'
moleculesTSC_List = 'Acetaldehyde'
moleculeSignalTSC_List = 1.66945
massNumberTSC_List = 29
moleculeConcentrationTSC_List = 0.05	#pressure can also be used in subsitute
unitsTSC = 'bar'	#the units will not be used in calculations so any units may be used

calculateUncertaintiesInConcentrations = False

#//Output Files//
#the last section designates the various output files, all are suppose to be csv values
#If files do not exist they will be generated
preProcessedDataOutputName= 'PreProcessedData.csv'
resolvedScaledConcentrationsOutputName = 'ScaledConcentrations.csv'
scaledConcentrationsPercentages = 'ScaledConcentrationPercentages.csv'
concentrationsOutputName= 'ResolvedConcentrations.csv'
simulatedSignalsOutputName= 'SimulatedRawSignals.csv'

#Only used in iterative analysis
TotalConcentrationsOutputName = 'TotalConcentrations.csv'


ExportAtEachStep = 'yes'
generatePercentages = 'no'

#Below is where __var_list__ used to be.