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

