Started with Test8 from the uncertainties tests, then worked on it to make an SLS which has a few chosen molecules and can test implicit feature.

The particular choices made here during solving are unusual. It is a fictitious example to that creates a more noticeable result than most scenarios.


During teh solving, we see from file ExportedSLSUniqueMoleculesAndChosenMassFragments that these massses are used: 
70	2butenalE(crotonaldehyde)
44	2buten1ol(crotyl alcohol)
54	13butadiene

From looking at ConvertedSpectra we see there is significant overlap between these molecules.

Exported8ReferenceThresholdFilter.csv shows:
Mass	2buten1ol(crotyl alcohol)	2butenalE(crotonaldehyde)	13butadiene
44	10.78107811	0	0
54	9.130913091	0	95.19951995
70	0	82.00820082	0
71	6.96069607	8.00080008	0

Showing that there is a clear solving order of Butenal, then Butenol, then Butadiene.
The mass 54 of 2butenalE(crotonaldehyde) has been made zero by the filtering. As a consequence, if 2butenalE(crotonaldehyde) Concentration is positive, then the 1,3-butadiene has been overestimated with the current settings.  

The mass 70 of 2buten1ol(crotyl alcohol) has also been made zero by filtering, which means that 2butenalE(crotonaldehyde) is over-estimated, also.

In test_1, the implicitSLS feature is set to False.
In test_2, the implicitSLS feature is set to True.  <-- the outcome is a bit strange with butenal (the first molecule solved) going from positive concentration in test_1 to negative in test_2. In looking into why, we find that the decimal ratio of correction for 2butenalE(crotonaldehyde) from 2buten1ol(crotyl alcohol) has changed due to mass 70 filtering out. Initially, the big diffference for butenal at the end is shocking, but looking at the concentrations that come out of test_1, we see that the 2buten1ol(crotyl alcohol) is on the order of 1000 times more in concentration for this example, which is enough that more than compensates for any supposed  2butenalE(crotonaldehyde).  Based on what's written above, we expect the 2butenalE(crotonaldehyde) to go down, and the 13butadiene to go down also. In this specific case, 2butenalE(crotonaldehyde) actually goes from slightly positive in concentration to very negative.  We also see that (importantly) the uncertainties for butenal are much higher in test 2.

#NOTE: There are two cases to consider, one where we should use the original concentration (currently implemented) and one where should use the revised concentration during recursion. What is best dependso on subtraction order etc. The current implementation works even with "circular" cases and that was one of the major reasons it was chosen.