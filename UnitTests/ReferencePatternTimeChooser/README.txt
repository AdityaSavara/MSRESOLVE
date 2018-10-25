ReferencePatternTimeChooser Unit Tester ReadMe

This feature allows the user to use different reference fragmentation patterns during data analysis in MS-RESOLV-SG

The feature consists of three functions:SelectReferencePattern, InterpolateReferencePatterns, PrepareReferenceObjectsAndCorrectionValues, 

SelectReferencePattern checks the current time in data analysis and determines whether or not to continue using the current reference file, use the next reference file, or interpolate the two reference files
It returns an index that tells the program which file to use, and the current reference data

InterpolateReferencePattern will interpolate two reference files if the user has a gap in their defined time range
It returns the interpolated reference object

PrepareReferenceObjectsAndCorrectionValues preprocesses reference objects and matches their correction values.  It is called for the interpolated reference object in SelectReferencePattern


test_1.py tests the use of two reference files by using three test input files and running MSRESOLVE.main() three times and reading the exported scaledConcentrations.csv files.
The experiment data used is truncated only including time ranges of 300 to 309
test_1_inputA.py uses AcetaldehydeNISTRefMixed.csv
test_1_inputB.py uses AcetaldehydeNISTRefMixedEdit.csv
test_1_inputC.py uses AcetaldehydeNISTRefMixed.csv from 300 to 303 and AcetaldehydeNISTRefMixedEdit.csv from 306 to 309.
The scaled concentrations are stored as scaledConcentrationAArray, scaledConcentrationBArray, and scaledConcentrationCArray
The scaled concentrations resulting from test_1_inputC from 300 to 303 are compared to the scaled concentrations resulting from test_1_inputA from 300 to 303
The scaled concentrations resulting from test_1_inputC from 306 to 309 are compared to the scaled concentrations resulting from test_1_inputB from 306 to 309

test_2.py tests the use of the interpolater
The same experimental data from test_1 is used.
Here MSRESOLVE.main() was run once using AcetaldehydeNISTRefMixed.csv and the standardized reference signals were exported
Then MSRESOLVE.main() was run again using AcetaldehydeNISTRefMixedEdit.csv and the standardized reference signals were exported
HandInterpolatedReferenceData.csv contains the interpolation of these files at time 304.8878 (the last time before 305)
Data Analysis is run using the time range 300 to 305 so at 305, currentReferenceData (a global variable) can be accessed in test_2.py
Comparing the standardized_reference_intensities of the two, there were three small rounding errors that was close enough to overwrite the expected_results

test_3.py tests simulated raw data signals when using two reference files
The test is run exactly the same as test_1 but instead of comparing scaledConcentrations, the test compares simulatedRawSignals