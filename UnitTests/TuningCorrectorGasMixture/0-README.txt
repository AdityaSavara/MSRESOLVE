There is a file called "ExtractedReferencePattern.csv" which came from actual calibration measurements on a mass spectrometer.
There is a file called "LiteratureReference.csv" which came from downloading the JDXConverter and NIST spectra for JDX.

For this example, we have created a fictional case of a gas mixture by multiplying the patterns in the ratio of ethane:ethene:ethyne of 10:3:1.
This is in the file TuningCorrectorGasMixtureHypotheticalReferenceMeasuredSignals.csv



Proc Doc type things:
 -- First added at line 5351 
     if len(G.UserChoices['measuredReferenceYorN']['tuningCorrectorGasMixtureMoleculeNames']) > 0:
 -- Now make a reference object from LiteratureReference.csv so that simulated signals can be made.
 -- GenerateReferenceDataList(referenceFileNamesList,referenceFormsList,AllMID_ObjectsDict={}):
 -- Need to be careful because "ReferenceInputPreProcessing" is being used.
  
     
     
 -- made  Collected_Data.csv the transpose of TuningCorrectorGasMixtureHypotheticalReferenceMeasuredSignals.csv
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