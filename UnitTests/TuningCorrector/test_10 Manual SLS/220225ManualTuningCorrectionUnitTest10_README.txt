Unit Test: test_10.py

Test_10 Quick Summary:              (Note: It's better to read this "README" file then to read the ProcDoc)

1)	Unit test_5 was modified to include SLS

2)	Our original reference pattern was scaled to 100 for each molecule, all intensities below 4 manually filtered out

3)	Using the "180613MadixKoMSCorrectionsFactor.." excel file, we input the scaled molecule patterns (scaled to 100 for each molecule) from the our scaled original pattern into MadixKo to get our correction factors (Non-Tuning Corrected). We must run MSRESOLVE to know which significant Mass is used. This will tell us which correction factor we need to use for each molecule.

4)	We acquire a Nist refernce pattern for the molecules used. These patterns are then scaled to 100 for each molecule.

5)	Using the "TuningCorrectorGasMixtureHypotheticalReferenceMeasuredVsSimulated"  excel file as a reference, we are able to create a polynomial fit equation between the original reference pattern and the Nist reference pattern. This equation can be used to create ratios for each molecular mass that will then be multiplied to the Nist reference pattern of "1Butanal" to get a Nist-TuningCorrected pattern for "1Butanal". This pattern will also be standardized to 100 and have intensities less than 4 removed.

6)	We now create a Mixed pattern, this reference pattern will include the original reference patterns for Acetaldehyde, CO, Ethylene (Ethene), Ethanol, Crotyl Alcohol, H2, H2O, and the Nist-TuningCorrected pattern for "1Butanal". This will give us our Mixed reference pattern.

7)	Using the "180613MadixKoMSCorrectionsFactor.." excel file, we input the molecule patterns from our Mixed reference pattern to get our new correction factors. We must run MSRESOLVE to know which significant Mass is used. This will tell us which correction factor we need to use for each molecule.

8)	We now perform SLS manually by picking a time (176.848 for this example), and subtracting each molecule signal and their contributions for each mass using our correction factor. We will obtain our solved concentration for each molecule.
