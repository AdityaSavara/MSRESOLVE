This is the unit test for iterative analysis.
test_1.py is the test file.

Basically, a non_iterative analysis is done using test_1_initial_input_noniterative.py

Then, an iterative analysis is done using test_1_initial_input_iterative.py.

This test was carefully crafted, there are some comments inside test_1.py
This test file tests the iterative analysis feature's compatibility with the use of multiple reference patterns at different time ranges where there is a gap in between time ranges requiring interpolation.
The reference file AcetaldehydeNISTRefMixed2.tsv was used from time 176 to time 200.  The reference file AcetaldehydeNISTRefMixed2Edit.csv was used from time 250 to 400.  So interpolation of reference patterns is required from 200 to 250.
The only difference between the data in these reference files is the signal of Acetaldehyde at m28.  AcetaldehydeNISTRefMixed2.tsv has a signal of 6103 at m28 while AcetaldehydeNISTRefMixed2Edit.csv has a signal of 6000 at m28.

When using iterative analysis, correction factors at interpolated points are calculated by interpolating correction values from the two reference files directly.  When not using iterative analysis, the signals from two reference patterns are interpolated and the correction values are determined from the interpolated signals.
So there is a slight difference in the data in scaledConcentrations.csv (output used to check non-iterative results) and TotalConcentrations.csv (output used to check iterative results).  Adding a relative tolerance of 1E-4 and absolute tolerance of 1E-8 takes care of any rounding.

