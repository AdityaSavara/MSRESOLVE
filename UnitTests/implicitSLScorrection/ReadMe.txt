Started with Test8 from the uncertainties tests, then worked on it to make an SLS which has a few chosen molecules and can test implicit feature.

During teh solving, we see from file ExportedSLSUniqueMoleculesAndChosenMassFragments that these massses are used: 2buten1ol(crotyl alcohol)	2butenalE(crotonaldehyde) have 70 and 71, 13butadiene has 54
From looking at ConvertedSpectra we see there is significant overlap between these molecules.
The order solved is: 2butenalE(crotonaldehyde) 2buten1ol(crotyl alcohol)		13butadiene

Exported8ReferenceThresholdFilter.csv shows:
Mass	2buten1ol(crotyl alcohol)	2butenalE(crotonaldehyde)	13butadiene
54	9.130913091	0	95.19951995
70	0	82.00820082	0
71	6.96069607	8.00080008	0

Showing that there is a clear solving order of Butenal, then Butenol, then Butadiene.
The mass 54 of 2butenalE(crotonaldehyde) has been made zero by the filtering. As a consequence, if 2butenalE(crotonaldehyde) Concentration is positive, then the 1,3-butadiene has been overestimated with the current settings.  

The mass 70 of 2buten1ol(crotyl alcohol) has also been made zero by filtering, which means that 2butenalE(crotonaldehyde) is over-estimated, also.

In test_1, the implicitSLS feature is set to False.
In test_2, the implicitSLS feature is set to True.  <-- the outcome is a bit strange, because we find that the decimal ratio of correction for 2butenalE(crotonaldehyde) from 2buten1ol(crotyl alcohol) due to mass 70 filtering out. Initially, this is shocking, but looking at the concentrations that come out of test_1, we see that the 2buten1ol(crotyl alcohol) is on the order of 1000 times more in concentration for this example, which is enough that more than compensates for any supposed  2butenalE(crotonaldehyde).  Based on what's written above, we expect the 2butenalE(crotonaldehyde) to go down, and the 13butadiene to go down also. In this specific case, 2butenalE(crotonaldehyde) actually goes from slightly positive in concentration to very negative.  Though this is just a fictitious example.


#NOTE: I think there are two cases to consider, one where we should use the original concentration (currently implemented) and one where should use the revised concentration during recursion. Need to think about that.