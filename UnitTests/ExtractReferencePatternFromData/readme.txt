ExtractReferencePatternFromData Unit Test

This function has the following inputs:
ExperimentData -- Class Object containing the collected data
ReferenceData -- Class Object containing the reference pattern
rpcMoleculesToChange -- a list of strings containing the name of the molecule to change
rpcMoleculesToChangeMF -- a list of lists containing the molecules' mass fragments that the user wants to change based on collected data
rpcTimeRanges -- the time range within the collected data the user wants to extract data

The output is an edited copy of reference data.

What this function does is it takes the average intensity of two (or more) user specified mass fragments during a specified time range.  The first mass fragment selected is the "base" fragment.
The reference signal of a specified mass fragment for a specified molecule is multiplied by the ratio of the average intensity of a specified mass fragment to the average intensity of the base mass fragment.

For the Unit Tests:
Test 1 uses a truncated piece of 2-CrotAcetExp#2.csv and the signals for m57 and m70 have been edited to where m70 is twice m57.
The variables input are:
rpcMoleculesToChange = ['Crotyl Alcohol']
rpcMoleculesToChangeMF = [[57,70]]
rpcTimeRanges = [[569,577]]
So here m70 for Crotyl Alcohol should be twice what it was before since the ratio of the average of m70 to the average of m57 is approximately 2.

Test 2 uses the same input as before but now
rpcMoleculesToChangeMF = [[70,57]]
So now m57 for Crotyl Alcohol should half what it was since the ratio of the average of m57 to the average of m70 is approximately 0.5

Test 3 uses the full 2-CrotAcetExp#2.csv and the original signals
Variable inputs are:
rpcMoleculesToChange = ['CO2']
rpcMoleculesToChangeMF = [[28,44]]
rpcTimeRanges = [[300,600]]
The reference signal for mass fragement 44 for CO2 changed from 10456 to 3551 since the ratio of the average signal of m44 to the average signal of m28 from 300 to 600 is about 3

Test 4 uses the full 2-CrotAcetExp#2.csv and the original signals
Variable inputs are:
rpcMoleculesToChange = ['Crotyl Alcohol','CO2']
rpcMoleculesToChangeMF = [[57,70],[28,44]]
rpcTimeRanges = [[300,600],[300,600]]
So here the reference signal for Crotyl Alcohol at mass fragment 70 will be replaced with the product of itstelf and the ratio of the average intensity of m70 from time 300 to 600 to the average intensity of m57 from time 300 to 600
Also the reference signal for CO2 at mass fragment 28 will be replaced with the procut of itself and the ratio of the average intensity of m44 from time 300 to 600 to the average intensity of m28 from time 300 to 600