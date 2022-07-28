This is the unit test for iterative analysis.
test_1.py is the test file.

Basically, a non_iterative analysis is done using test_1_initial_input_noniterative.py

Then, an iterative analysis is done using test_1_initial_input_iterative.py.

This test was carefully crafted, there are some comments inside test_1.py
This test file tests the iterative analysis feature's compatibility with the use of multiple reference patterns at different time ranges.
The reference file AcetaldehydeNISTRefMixed2.tsv was used from time 176 to time 250.  The reference file AcetaldehydeNISTRefMixed2Edit.csv was used from time 250 to 400.
The only difference between the data in these reference files is the signal of Acetaldehyde at m28.  AcetaldehydeNISTRefMixed2.tsv has a signal of 6103 at m28 while AcetaldehydeNISTRefMixed2Edit.csv has a signal of 6000 at m28.
