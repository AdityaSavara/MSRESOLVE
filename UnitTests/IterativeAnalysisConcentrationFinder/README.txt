This is the unit test for iterative analysis.
test_1.py is the test file.

Basically, a non_iterative analysis is done using test_1_initial_input_noniterative.py

Then, an iterative analysis is done using test_1_initial_input_iterative.py.

This test was carefully crafted, there are some comments inside test_1.py
This test file tests iterative analysis's compatibility with concentration finder where the known concentrations at certain signals pertain to separate molecules.
The ResolvedConcentrations.csv file (concentrations from non-iterative run) is compared to TotalConcentrations.csv (concentrations from iterative runs)