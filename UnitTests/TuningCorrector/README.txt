TuningCorrector Unit Tester ReadMe

This feature is for correcting reference patterns between different mass spectrometer tunings. In practice, it is for externally adjusted reference patterns (like NIST) to better align with one's own data (converting external spectra to match your spectrometer's tuning). However, it can also be used for adjusting spectra from one's own spectrometer to match tuning from external ones (such as NIST's).


The way the feature works is it looks at reference patterns collected from two mass spectrometres, then it uses a fit to make polynomial based tuning correction so that it can make the reference pattern from one mass spectrometer look like it was collected on the other one.

#apparently need to have dataAnalysis on to use this feature (ideally should not need to, but as of Sep 2019, do need to).

test_1.py takes the ReferenceCollected.csv and ReferenceLiterature.csv

ReferenceCollected.csv does not have as many molecules as ReferenceLiterature.csv and the higher masses have lower intensity in ReferenceCollected.csv
So the feature uses a polynomial function and applies it to *all* molecules in ReferenceLiterature.csv to make it look more like ReferenceCollected.csv.

In the test_1.py, the reference threshold filter is off.
In  test_2.py, the reference threshold filter is on.

In test_3.py, the two files are referenceFileExistingTuning = ['ReferenceLiterature.csv','xyyy'] and referenceFileDesiredTuning =['ReferenceCollected.csv','xyyy'].  The original reference file is set as ReferenceCollected.



In many real applications of this feature, what is desired is to predict from an external reference what the fragmentation pattern would be on one's own spectrometer.  In that situation, the "ReferenceCollected.csv" is the desired pattern one and the "ReferenceLiterature.csv" is the existing pattern to be adjusted. These names may become further adjusted to "PatternToMatch" and "PatternToAdjust" or something like that.