This is the ReadMe for the concentrationFinder unit test.

ConcentrationFinder uses signals of known molecules' concentrations to determine conversion factors that are used to convert scaled to CO concentrations to concentrations of a user determined unit set.

The unit test uses a reference pattern containing Acetaldehyde and Acetaldehy_easy_to_ionize that share the exact same data with the exception of their ionization factors and signal intensities at m29 and m29.2
Acetaldehyde has an ionization factor of 1, a signal intensity of 9999 at m29, and a signal intensity of 0 at m29.2
Acetaldehyde_easy_to_ionize has an ionization factor of 2, a signal intensity of 0 at m29, and a signal intensity of 9999 at m29.2

The collected data is truncated with only m29 and m29.2 present and their signals are 1 for each data point.  Since m29's only contribution is from Acetaldehyde and m29.2's only contribution is from Acetaldehyde_easy_to_ionize, we would expect the scaled concentrations to be the same for both molecules.
However, the ratio of ionization factors as 2 so we expect the scaled concentrations to differ by a factor of 2.

We convert to concentration based on our "known" concentration of Acetaldehyde (0.05 bar when the signal is 1.66945) and expect it to also differ by a factor of 2.

Running test_1.py runs MSRESOLVE.main() and then gets the concentrationsarray from the global resultsObjects dictionary.
Taking the ratio of the resolve Acetaldehye to Acetaldehyde_easy_to_ionize, we get 1.9982073753308842 which is close to 2 as expected.


Test_2.py uses AcetaldehydeNISTRefMix2_test_2.csv and 2-CrotAcetExp#2Truncated2.csv.  The reference file contains only Acetaldehyde with known ionization factor of 1.  The collected data file contains one mass fragment (m29) with a signal of 1 at each data point.
The same reference file is used twice but the 'known' concentration of Acetaldehyde differs at different times.
From times 1 to 4, the 'known' concentration of Acetaldehyde is 0.05 bar at a signal of 1.66945.
From times 5 to 8, the 'known' concentration of Acetaldehyde is 0.1 bar at a signal of 1.66945.
Since the 'known' concentration differs by 2, we expect resolved concentrations to also differ by 2 since the collected data has a uniform signal of 1 and the two reference patterns are identical.

Test_3.py uses AcetaldehdyeNISTRefMix2_test_1.csv and 2-CrotAcetExp#2Truncated.csv (the same reference file and collected data file as test_1.py.
The main difference between this test and test_1.py is the conversion factor used for Acetaldehyde_Easy_To_Ionize is overwritten based on values input by the user.
In this case we say that we know the concentration of Acetaldehyde to be 0.05 bar at a m29 signal of 1.66945 and Acetaldehyde_Easy_To_Ionize to be 0.15 bar at a m29.2 signal of 1.66945.
Since the ratio of 0.05 to 0.15 is a factor of 3, we expect the ratio of resolved concentrations of Acetaldehyde_Easy_To_Ionize to Acetaldehye to be 3.
