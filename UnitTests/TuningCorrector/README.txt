TuningCorrector Unit Tester ReadMe

This feature is for correcting reference patterns between different mass spectrometer tunings. In practice, it is for 'correcting' ones own mass spectrometer's tuning with the assumption that the NIST tuning is 'correct'.


The way the feature works is it looks at reference patterns collected from two mass spectrometres, then it uses a fit to make polynomial based tuning correction so that it can make the reference pattern from one mass spectrometer look like it was collected on the other one.

#apparently need to have dataAnalysis on to use this feature (should not need to, but as of Sep 2019, do need to.

test_1.py takes the ReferenceCollected.csv and ReferenceLiterature.csv

ReferenceLiterature.csv does not have as many molecules as ReferenceCollected.csv and the higher masses have lower intensity in ReferenceLiterature.csv.
So the feature uses a polynomial function and applies it to *all* molecules in ReferenceCollected.csv to make it look more like ReferenceLiterature.csv.

In the test_1.py, the reference threshold filter is off.
In  test_2.py, the reference threshold filter is on.