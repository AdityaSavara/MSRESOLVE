TuningCorrector Unit Tester ReadMe

This feature is for correcting reference patterns between different mass spectrometer tunings. There are actually two aspects of this feature:

a)TuningCorrectorPattern which changes the pattern of External reference spectra to match how it would appear on one's own mass spectrometer. This is for solving which molecules are present and subtracting the signals.  In practice, it is for externally adjusted reference patterns (like NIST) to better align with one's own data (converting external spectra to match your spectrometer's tuning). However, it can also be used for adjusting spectra from one's own spectrometer to match tuning from external ones (such as NIST's). The way the feature works is it looks at reference patterns collected from two mass spectrometres, then it uses a fit to make polynomial based tuning correction so that it can make the reference pattern from one mass spectrometer look like it was collected on the other one.


b)TuningCorrectorIntensity which changes the sensitivityCorrectionValues as a function of mass fragment so that an accurate relative concentration of species can be determined. This involves two steps, each using a polynomial correction for the direction of converting one's own mass spectrometer patterns into the direction of the standard corrections. First, convert the patterns into the standard tunining version (that way the fragmentation pattern would be correct for summing up the fragments), that way standard tuning correction factors are calculated. The first step is necessary because otherwise the summation term and calculated ratios of the fragments will be incorrect.  Then, as a second step, we applying a tuning correction (actually the same factor) to the final correction factor, which is the equivalent of adding a tuning correction to the observed intensity. So ultimately, we have CorrectionFactor_standardTuning * Signal_standardTuningEquivalent = Concentration_standardTuning.

The typical usage is shown in test_4.py

#apparently needed to have dataAnalysis on to use this feature (as checked Sep 2019, ideally should not need to).

Note: Before approximately Sept 23rd-30th 2021, the regular reference file was listed in the unit tests as the existing tuning. After Sept 30th 2021, during standard usage, the 'existing' tuning is the external (NIST/Literature) tuning.  The following files now get exported to try to bring clarity: 'ExportedReferencePatternOriginal" "ExportedReferencePatternExisting" "ExportedReferencePatternExternal" "ExportedReferencePatternExternalTuningCorrected"  "ExportedReferencePatternMixed"

The file AcetaldehydeMeasured.csv is actually a mixed reference pattern. However, the columns source fields have been renamed to say "Measured". The tuning correction tests here are just a 'toy' model reference pattern to check the feature. In real life, the measured reference pattern will have more difference than the literature reference.

test_1.py has AcetaldehydeMeasured.csv and ReferenceCollected.csv with the same tuning. It applies a TuningCorrection to AcetaldehydeMeasured.csv (and also to ReferenceCollected.csv).  The Desired tuning file is ReferenceLiterature.csv. As of Sept 30th 2021, this is not the typical usage (this is essentially reverse of the current typical usage).
ReferenceLiterature.csv does not have as many molecules as ReferenceCollected.csv and the higher masses have lower intensity in ReferenceLiterature.csv (for one of the molecules, crotyl alcohol). So the feature uses a polynomial function and applies it to *all* molecules in ReferenceCollected.csv to make it look more like ReferenceLiterature.csv (lowers the intensity).

In the test_1.py, the reference threshold filter is off.
In  test_2.py, the reference threshold filter is on.

In test_3.py, the two files are referenceFileExistingTuning = ['ReferenceCollected.csv','xyyy'] and referenceFileDesiredTuning =['ReferenceLiterature.csv','xyyy'].  The original reference file is set as ReferenceLiterature. As of Sept 30th 2021, this is not the typical usage (this is essentially reverse of the current typical usage).

test_4.py is a copy of test_1.py, only now the returnMixedPattern feature is set to true, so ReferenceLiterature is tuned to match ReferenceCollected. This file **is** the typical usage.
test_5.py is a copy of test_4.py, only now the desired pattern is set as blank, which should give the same output as test_4. 
test_6.py is a copy of test_5.py, only now the existing pattern is set as blank, but the new ReferencePatternStandard is populated, so that the existing pattern will be populated from that one. This also means that there is a tuningCorrectionIntensity feature usage. This test had output that matched test_5.py exactly before the tuningCorrectionIntensity feature was implemented. In test_6.py, the effects of tuningCorrectionIntensity are quite small. So test_7 and test_8 were created for checking the feature properly.

test_7.py and test_8.py are copies of test_6.py but uses referenceThreshold filtering , specific mass fragments, and also SLSUniqueExport so that mass 70 ends up being used. We see that ExportedSLSUniqueMoleculesAndChosenMassFragments indicates the below, so that we know Crotonaldehyde is solved with mass 70.

2	H2
18	H2O
45	Ethanol
70	(E) 2-Butenal (Crotonaldehyde
31	Crotyl Alcohol

In test_7.py, only the ExternalTuningCorrection feature is used.
In test_8.py, the tuningCorrectionIntensity feature is used by populating the variable referenceFileStandardTuning.

In both test_7 and test_8, the butenal is solved for using mass 70.
test_7 and test_8 should have the same reference pattern at the end: they do.
test_7 and test_8 should have different concentrations at the end: they do.
The test_8 scaled concentrations has a lower butenal concentration.

By embedding a print statement inside MSRESOLVE, this 'manual' printing showed that the intensity correction factor for mass 70 is 0.4085935090718888 (because The Reference Collected file is a little higher at mass 70 than the Reference Literature). When only the second step of the tuningCorrectionIntensity feature was implemented, we could see that test_8.py had Butenal of 0.017563494273429126 in comparison to the test_7.py butenal of 0.042985250336757, which is consistent with the same 0.408 factor. This example is not a realistic real world example because normally there will be many molecules showing the same trend in tuning. Here, just one molecule was adjusted between the patterns such that the "average" tuning correction does not really work well on average. This is just a unit test, it is not really an example case.  After the implementation of tuningCorrectionIntensity the value at mass 70 changes to 0.015542924457281676.  more significantly, the CO solved value as well as some of the other small molecules increase in concentration quite a bit. 

#FIXME: test_9.py is with createMixedTuningPattern turned off and also uses a shortened measured pattern (this test is not realistic and was made just to check that MSRESOLVE would lengthen arrays as needed). Test_9.py was not checked to see what the outcomes were on the resolved concentrations, just that a mixed pattern was not made and that the analysis ran to completion.  However, inspection of the scaled concentrations that come out for test_9.py show only one molecule, which is not the expected behavior.

***
In many real applications of this feature, what is desired is to predict from an external reference what the fragmentation pattern would be on one's own spectrometer.  In that situation, the "ReferenceCollected.csv" is the desired pattern one and the "ReferenceLiterature.csv" is the existing pattern to be adjusted. These names may become further adjusted to "PatternToMatch" and "PatternToAdjust" or something like that.