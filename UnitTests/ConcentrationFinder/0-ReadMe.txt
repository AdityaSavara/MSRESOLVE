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

Test_2 is testing the concentrationFinder feature with numerous reference patterns.  Different concentration scaling factors can be used at different times.
Test_2.py uses AcetaldehydeNISTRefMix2_test_2.csv and 2-CrotAcetExp#2Truncated2.csv.  The reference file contains only Acetaldehyde with known ionization factor of 1.  The collected data file contains one mass fragment (m29) with a signal of 1 at each data point.
The same reference file is used twice but the 'known' concentration of Acetaldehyde differs at different times.
From times 1 to 4, the 'known' concentration of Acetaldehyde is 0.05 bar at a signal of 1.66945.
From times 5 to 8, the 'known' concentration of Acetaldehyde is 0.1 bar at a signal of 1.66945.
Since the 'known' concentration differs by 2, we expect resolved concentrations to also differ by 2 since the collected data has a uniform signal of 1 and the two reference patterns are identical.

Test_3 is testing the concentrationFinder feature with different scaling factors for different molecules.  This test also demonstrates any unlisted molecule (in the TSC list from the user input file) will use a concentration factor based on the first listed molecule's known concentration.
Test_3.py uses AcetaldehdyeNISTRefMix2_test_3.csv and 2-CrotAcetExp#2Truncated3.csv. The purpose of this test is for having separate concentration factors for separate molecules. This is the same reference file and collected data file as test_2.py with Acetaldehyde_copy added to the reference data and m29.3 added to the collected data.  Acetaldehyde_copy is a copy of Acetaldehyde but the signal at m29 is 0 and at m29.3 is 9999.  m29.3, like m29 and m29.2, has a signal of 1 at each time point.
In this case we use the test input file to indicate that we "know" the concentration of Acetaldehyde to be 0.05 bar at a m29 signal of 1.66945 and Acetaldehyde_Easy_To_Ionize to be 0.15 bar at a m29.2 signal of 1.66945.
Since the ratio of 0.05 to 0.15 is a factor of 3, we expect the ratio of resolved concentrations of Acetaldehyde_Easy_To_Ionize to Acetaldehye to be 3.
In this test, Acetaldehyde_copy will use the same conversion factor as Acetaldehyde, because all unspecified molecules use the first molecule listed's scaling factor, so the resolved concentrations of the two molecules should be the same, or in other words the ratio of the two will be 1.

Test_4 is testing the concentrationFinder feature with numerous reference patterns that require interpolating due to a gap in between time ranges.  
Test_4.py is set up almost indentically to test_2.py where the same reference file is used for two different time ranges and the same collected data file is used.
The first time range is from 1 to 1 and the last time range is from 8 to 8.  By doing this, every molecule that is not in the first or last time point will be solved by using a conversion factor that was obtained via interpolation.
Like test_2, the concentration at point 8 will be twice the concentration at point 1.
To check that our interpolation is accurate, we can interpolate the concentration of our molecule from the first time point and the last time point at each individual time point (i.e. at time point 3, interpolate the concentration of time 3 between the concentrations from time points 1 and 8)
 