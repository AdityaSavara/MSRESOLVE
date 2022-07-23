There is a file called "ExtractedReferencePattern.tsv" which came from actual calibration measurements on a mass spectrometer.
There is a file called "LiteratureReference.tsv" which came from downloading the JDXConverter and NIST spectra for JDX.

For this example, we have created a fictional case of a gas mixture by multiplying the patterns in the ratio of ethane:ethene:ethyne of 10:3:1.
This is in the file TuningCorrectorGasMixtureHypotheticalReferenceMeasuredSignals.tsv



Proc Doc type things:
 -- First added at line 5351 
     if len(G.UserChoices['tuningCorrection']['tuningCorrectorGasMixtureMoleculeNames']) > 0:
 -- Now make a reference object from LiteratureReference.tsv so that simulated signals can be made.
 -- GenerateReferenceDataList(referencePatternsFileNamesList,referencePatternsFormsList,AllMID_ObjectsDict={}):
 -- Need to be careful because "ReferenceInputPreProcessing" is being used.
  
     
     
 -- made  Collected_Data.csv the transpose of TuningCorrectorGasMixtureHypotheticalReferenceMeasuredSignals.tsv
 -- The reason Collected_Data was edited is because we need to compare the masses of the measured data and the reference data to remove unnecessary masses when creating a reference object.

###NEEDED TO MAKE THE NIST REFERENCE PATTERN THE LITERATURE PATTERN BECAUSE OTHERWISE THE EXPERIMENTAL DATA MASSES DID NOT GET TRIMMED.###



Creation of the fictional gas mixture was performed in FictionalGasMixtureCreationFile, Because only the ratio matters for TuningCorrectorGasMixture.


***

Test_2.py adds in the tuning correction intensity feature by use of the Standard Reference pattern (but has uncertainties off)
Test_3.py adds in the tuning correction intensity feature by use of the Standard Reference pattern (and has uncertainties on)

Test_2.py and Test_3.py outputs have identical concentrations for ethene, ethane, and 1butanal.
Test_2.py and Test_3.py outputs have differing concentrations for ethyne.

If we look at the SLSUniqueMoleculesAndChosenMassFragments, we see that the reason is that test_3.py switches to using m24 for ethyne rather than m26.
Test_4.py is a copy of test_3.py, but the SLS solving has been changed to focus on largest reference fragment. Accordingly, test_4.py has solving and output that matches test_2.py

***
To evaluate whether the tuning corrector intensity feature is working correctly, we note:

In 1, 2, and 4 the concentrations are solved based on:

30	ethane
15	1butanal
28	ethene
26	ethyne

From the file TuningCorrectorGasMixtureHypotheticalReferenceMeasuredVsSimulated.xlsx, we see that at mass 15 the measured ratio is much higher than the literature.  That means that the tuning correction for the pattern will push the tuning corrected 1butanal m15 to be higher.
Comparing ExportedReferencePatternMixed.tsv and ExportedReferencePatternStandardForCorrectionValuesMixedStandardTuning.tsv, we see that m15 relative intensity has gone from 5.4 in Standard Tuning to 7.32 in the "ExportedReferencePatternMixed" where it is tuning corrected.  The tuning corrector intensity feature should move in the opposite direction. That means it should bring the ScaledConcentration of 1butanal down in test_2.py relative to test_1.py. Test_4.py has the same scaled concentrations as test_2.py.  That is the trend observed, but test_4.py and test_2.py have negative concentrations for 1butanal, which is a bit surprising.  To look at a more clear example, we can look at the first molecule solved by SLS which is ethane, solved by m30. Based on what is in  TuningCorrectorGasMixtureHypotheticalReferenceMeasuredVsSimulated.xlsx, the pattern should go down at m30 from standard tuning to tuning corrected (and it does, in ExportedReferencePatternStandardForCorrectionValuesMixedStandardTuning.tsv and ExportedReferencePatternExternalTuningCorrected.tsv) : Thus we would expect the concentration to go up for ethane, when going from test_1.py to test_2.py. It does. Test_4.py has the same scaled concentrations as test_2.py.