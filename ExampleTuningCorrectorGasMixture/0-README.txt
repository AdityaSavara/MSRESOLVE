There is a file called "ExtractedReferencePattern.csv" which came from actual calibration measurements on a mass spectrometer.
There is a file called "LiteratureReference.csv" which came from downloading the JDXConverter and NIST spectra for JDX.

For this example, we have created a fictional case of a gas mixture by multiplying the patterns in the ratio of ethane:ethene:ethyne of 10:3:1.
This is in the file TuningCorrectorGasMixtureHypotheticalReferenceMeasuredSignals.csv


The reference pattern listed for analysis has internally collected spectra which are retained in the final pattern (ethyne, ethene, and ethane). This example takes an external 1butanal pattern (from NIST) and corrects it to match the collected spectra's spectrometer.

Proc Doc type things:
 -- First added at line 5351 
     if len(G.UserChoices['tuningCorrection']['tuningCorrectorGasMixtureMoleculeNames']) > 0:
 -- Now make a reference object from LiteratureReference.csv so that simulated signals can be made.
 -- GenerateReferenceDataList(referencePatternsFileNamesList,referencePatternsFormsList,AllMID_ObjectsDict={}):
 -- Need to be careful because "ReferenceInputPreProcessing" is being used.
  
     
     
 -- made  Collected_Data.csv the transpose of TuningCorrectorGasMixtureHypotheticalReferenceMeasuredSignals.csv
 -- The reason Collected_Data was edited is because we need to compare the masses of the measured data and the reference data to remove unnecessary masses when creating a reference object.

###NEEDED TO MAKE THE NIST REFERENCE PATTERN THE LITERATURE PATTERN BECAUSE OTHERWISE THE EXPERIMENTAL DATA MASSES DID NOT GET TRIMMED.###



Creation of the fictional gas mixture was performed in FictionalGasMixtureCreationFile, Because only the ratio matters for TuningCorrectorGasMixture.