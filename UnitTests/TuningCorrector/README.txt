TuningCorrector Unit Tester ReadMe

This feature is for correcting reference patterns between different mass spectrometer tunings. In practice, it is for externally adjusted reference patterns (like NIST) to better align with one's own data (converting external spectra to match your spectrometer's tuning). However, it can also be used for adjusting spectra from one's own spectrometer to match tuning from external ones (such as NIST's).


The way the feature works is it looks at reference patterns collected from two mass spectrometres, then it uses a fit to make polynomial based tuning correction so that it can make the reference pattern from one mass spectrometer look like it was collected on the other one.

#apparently need to have dataAnalysis on to use this feature (ideally should not need to, but as of Sep 2019, do need to).

NOTE: ON SEPT 23rd 2021, THE REGULAR INPUT FILE IS CONSIDERED 'EXISTING' TUNING. THIS WILL CHANGE WITH THE MIXED REFERENCE PATTERN UPDATE.
AFTER THAT, THE REGULAR INPUT FILE WILL BE CONSIDERED THE 'DESIRED' TUNING. IT MAY BE A GOOD TIME TO RENAME THE VARIABLES TO 'EXTERNAL' AND 'EXISTING' TUNING OR SOMETHING LIKE THAT. THIS MAY MAKE IT MORE CLEAR THAT THE EXTERNAL TUNING WILL BE CHANGED TO MATCH THE EXISTING ONE.
WILL EXPORT FILES AS 'ExportedReferencePatternOriginal" "ExportedReferencePatternExisting" "ExportedReferencePatternExternal" "ExportedReferencePatternExternalTuningCorrected"  "ExportedReferencePatternMixed"

test_1.py has AcetaldehydeNISTRefMixed2.csv and ReferenceCollected.csv with the same tuning. It applies a TuningCorrection to AcetaldehydeNISTRefMixed2.csv (and also to ReferenceCollected.csv).  The Desired tuning file is ReferenceLiterature.csv.
ReferenceLiterature.csv does not have as many molecules as ReferenceCollected.csv and the higher masses have lower intensity in ReferenceLiterature.csv (for one of the molecules, crotyl alcohol). So the feature uses a polynomial function and applies it to *all* molecules in ReferenceCollected.csv to make it look more like ReferenceLiterature.csv (lowers the intensity).

In the test_1.py, the reference threshold filter is off.
In  test_2.py, the reference threshold filter is on.

In test_3.py, the two files are referenceFileExistingTuning = ['ReferenceCollected.csv','xyyy'] and referenceFileDesiredTuning =['ReferenceLiterature.csv','xyyy'].  The original reference file is set as ReferenceLiterature.

test_4.py is a copy of test_1.py, only now the returnMixedPattern feature is set to true, so ReferenceLiterature is tuned to match ReferenceCollected.
test_5.py is a copy of test_4.py, only now the desired pattern is set as blank, which should give the same output as test_4.
test_6.py is a copy of test_5.py, only now the existing pattern is set as blank, but the new ReferencePatternStandard is populated, so that the existing pattern will be populated from that one. This also means that there is a tuningCorrectionIntensity feature usage. This test had output that matched test_5.py exactly before the tuningCorrectionIntensity feature was implemented.

test_7.py is a copy of test_6.py but uses referenceThreshold filtering , specific mass fragments, and also SLSUniqueExport so that mass 70 ends up being used. We see that ExportedSLSUniqueMoleculesAndChosenMassFragments indicates the below, so that we know Crotonaldehyde is solved with mass 70.

2	H2
18	H2O
45	Ethanol
70	(E) 2-Butenal (Crotonaldehyde
31	Crotyl Alcohol



In many real applications of this feature, what is desired is to predict from an external reference what the fragmentation pattern would be on one's own spectrometer.  In that situation, the "ReferenceCollected.csv" is the desired pattern one and the "ReferenceLiterature.csv" is the existing pattern to be adjusted. These names may become further adjusted to "PatternToMatch" and "PatternToAdjust" or something like that.