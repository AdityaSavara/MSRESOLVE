This directory test the Best Mass Fragment Chooser. 

NOTE: Now, by default, it uses ExtentOfSLSUniqueSolvable.py and does not do on the fly SLS, nor does it do the significance checks etc.

The Best Mass Fragment Chooser selects mass fragment combinaitons that have unique fragments and can therefore pass the SLS Unique Fragments solution method. It ranks the mass fragment combination in terms the sum of the significance factors across the reference data for the selected masss fragments and returns the fragment combinations having the largest significance sum. It also returns a reference file containing refernce data for only the specified mass fragments and molecules.

The user specified inputs are:
moleculesToMonitor, moleculesLikelihood, numberOfMassFragsToMonitor, referenceFileName, referenceForm, referenceIntensityThreshold=5, onTheFlySLS=False, keep_N_ValuesInRoughUniquenessCheck=1000, keep_N_ValuesInSignificanceFactorCheck=1000, finalNumberOfCombinationsToKeep=10

Pre-Checks
    keep_N_ValuesInRoughUniquenessCheck can be set as either an integer or False if the feature is not desired. It will preform a rough uniqueness check that will keep and N number of mass fragment combinations contain the largest number of zeros. The more zeros, the higher the chance of unique fragments.

    keep_N_ValuesInSignificanceFactorCheck can be set as either an integer or False if the feature is not desired. It will keep only an N number of mass fragment combinations contianing the largest sum of significance values.

onTheFlySLS runs SLS while iterating through all possible combinations of mass fragments. When either of the pre-checks are used it only preforms SLS when the mass fragment combinaiton is stored. onTheFlySLS is used when the pre-checks are not used.

There are 4 unit tests:
    1) Uses the default settings of the Best Mass Fragment Chooser to evaluate 3 mass fragments for Ethylene(Ethene) and Ethanol.
    2) 4 mass fragments are selected for Ethylene(Ethene), Ethanol, and Crotyl Alcohol. The pre-checks are turned off and onTheFlySLS used.
    3) 4 mass fragments are selected for Ethylene(Ethene), Ethanol, and Crotyl Alcohol. onTheFlySLS is used with pre-check. This is expected to give the same results as test 2.
	4) moleculesToMonitor=['Ethylene (Ethene)', 'Ethanol', 'Crotyl Alcohol' , 'Acetaldehyde' ] and 6 fragments. With pre-check.
	5) moleculesToMonitor=['Ethylene (Ethene)', 'Ethanol', 'Crotyl Alcohol' , 'Acetaldehyde' ] and 6 fragments. Without pre-check.
	
	The unit tests above "10" have useExtentOfSLSUniqueSolvable = True, and add progressively more molecules to check.

	