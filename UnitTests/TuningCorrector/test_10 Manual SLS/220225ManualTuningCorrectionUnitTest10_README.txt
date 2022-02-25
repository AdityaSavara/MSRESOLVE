Unit Test: test_10.py

Test_10 Quick Summary:              (Note: It's better to read this "README" file then to read the ProcDoc)

1)	Unit test_5.py was modified to include SLS. Once modified, test_5.py was renamed to test_10.py.

2)	Our original reference pattern "AcetaldehydeMeasured.csv" was scaled to 100 for each molecule, all intensities below 4 manually filtered out. This can be shown in "AcetaldehydeMeasured_Scaled_to_100.xlsx"

3)	Using the "180613MadixKoMSCorrectionsFactorExample_NO_HIDING.xlsx" file, we input the scaled molecule patterns (scaled to 100 for each molecule) from the our scaled original pattern file "AcetaldehydeMeasured_Scaled_to_100.xlsx" into MadixKo to get our correction factors (Non-Tuning Corrected). This is shown in the "180613MadixKo_Test_10.xlsx" file.

4)  	We must run test_10.py to know which significant Mass is used & SLS molecule order. This will tell us which correction factor we need to use for each molecule. You will need the measured csv file "2-CrotAcetExp#2". The original reference file "AcetaldehydeMeasured.csv" and the existing tuning reference file "ReferenceLiterature.csv". Once ran, we will get our SLS molecule order and chosen mass fragment in the export file "ExportedSLSUniqueMoleculesAndChosenMassFragments.csv". We will then rename this file to "Masses-Molecules_SLS_UnitTest_Test_10.xlsx".

5)	We acquire a Nist refernce pattern for the molecules used. This is the "ReferenceLiterature.csv" that was used in the last step.  These patterns are then scaled to 100 for each molecule in the file "ReferenceLiterature_Standardized_100". Zeros where then removed from the patterns. This can be found in the "NistRef_Zeros_Removed.csv"
	We must also acquire a measured reference file with zeros removed. Using the "ExportedReferencePatternDesiredOriginal.csv" file (This file will be acquired after runing test_10.py) we will then remove the zeros and unneeded molecules, Kept: (Acetaldehyde, (E) 2-Butenal, Ethylene, Ethanol, Crotyl Alcohol). This can be found in the "MeasureRef_Zeros_Removed.csv" file.

6)	Using the "TuningCorrectorGasMixtureHypotheticalReferenceMeasuredVsSimulated.xlsx"  excel file as a reference, we are able to create a polynomial fit equation between the measured reference pattern (MeasureRef_Zeros_Removed.csv) and the Nist reference pattern (NistRef_Zeros_Removed.csv). This equation can be used to create ratios for each molecular mass that will then be multiplied to the Nist reference pattern of "1Butanal" to get a Nist-TuningCorrected pattern for "1Butanal". The file with the created polynomial fit can be found in the "Created_Polynomial_Fit_Test_10.xlsx" file.  This pattern will also be standardized to 100 and have intensities less than 4 removed.

7)	We now create a Mixed pattern, this reference pattern will include the original reference patterns (from AcetaldehydeMeasured_Scaled_to_100.xlsx) for molecules Acetaldehyde, CO, Ethylene (Ethene), Ethanol, Crotyl Alcohol, H2, H2O, and the Nist-TuningCorrected pattern for "1Butanal" ("1Butanal" tuning corrected made in the "Created_Polynomial_Fit_Test_10.xlsx") . This will give us our Mixed reference pattern ("Mixed_Ref.csv").

8)	Using the "180613MadixKoMSCorrectionsFactorExample_NO_HIDING.xlsx" excel file, we input the molecule patterns from our Mixed reference pattern ("Mixed_Ref.csv") to get our new correction factors. This can be seen in the file "180613MadixKo_Test_10_TC_1butanal".

9)	We now perform SLS manually by picking a time from our collected signals file "2-CrotAcetExp#2", 176.848 for this example, and subtracting each molecule signal and their contributions (contributions are the intensites of each molecule as provided in the "Mixed_Ref.csv". Order of molecule subtraction and chosen mass fragment in "Masses-Molecules_SLS_UnitTest_Test_10.xlsx") for each mass using our correction factors found in "180613MadixKo_Test_10_TC_1butanal". Using the Madix & Ko formula, we can obtain our solved concentration for each molecule. Manual SLS can be viewed in the "Test_10_SLS.xlsx" file. 
