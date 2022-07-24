Test 1 and 2 use SLSUnique with no "collected data" uncertainties. There are different SLS weightings for these.
Test 3 uses inverse and has no "collected data" uncertainties.
Note: Test 3 has two possible outcomes: one if the uncertainties modules is present, and one if it is not. no "collected data" uncertainties.
Test 4 uses inverse and has collected data uncertainties as well as the absolute uncertainties feature for the reference file, based on the 'auto' uncertainties feature. This data is particularly not noisy, so the error is not that much increased.
Test 5 uses the collected data uncertainties and the reference data uncertainties for the SLS.
test 6 inverse
test 7
test 8 sls chosen molecules, chosen fragments
test 9
test 10 inverse which is a copy of test 3 except that it uses a file for the uncertainties as input
test 11 inverse chosen molecules
test 12 sls
test 13 sls, chosen molecules, chosen fragments. --> ends up using inverse after solving one molecule.
test 14 inverse with  tuning corrector and tuning corrector uncertainties that are propagated (unfortunately, values are arbitrary)
test 15 inverse with same tuning corrector as 14 but no tuning correctors to be propagated.  Comparing the relative uncertainties in the concentrations of of 14 and 15 shows that 14 does in fact have higher uncertainties.
test 16: copy of test_14 but with Mixed Reference pattern as True but with no external pattern provided. Output matches 14 as expected.  However, for bizarre reasons that are not clear, test_16.py fails with pytest during datasmoothing (while test_14.py does not) if using the base package without uncertainties. Accordingly, test_16 has been changed to a manual test.
test 17 dataForAnalysis file with _absolute_uncertainties based on test_6.py.  The _absolute_uncertainties are actually the same uncertainties as in test_6.py: they were exported manually (by inserting a savetxt statement into MSRESOLVE), and are stored in 0-test_6-rawsignals_absolute_uncertainties.csv.  That data were used as to create  20190817A-integer-data_absolute_uncertainties.csv, so that the same outputs should come from test_17 and test_6. That is what is observed, within numerical error.
test 18 dataForAnalysis file with _absolute_uncertainties set as 50% of the original signals, based on test_17.  The file 20190817A-integer-data18.csv is the same as the file without the "18" in the name, but the file 20190817A-integer-data18_absolute_uncertainties.csv has 50% of the original signal size as the uncertainties. The same scaled concentrations as test_6.py should arise, and the uncertainties should be at least half the size of the final concentrations. And that is indeed what happens: test_18 has the same concentration outputs as test_6.py, but has ~50% relative uncertainties (slightly above 50%) for each concentration.