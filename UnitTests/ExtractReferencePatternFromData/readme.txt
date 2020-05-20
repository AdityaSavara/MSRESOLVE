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
So here m70 will be altered so that it is twice m57 for crotyl alcohol

Test 2 uses the same input as before but now
rpcMoleculesToChangeMF = [[70,57]]
So now m57 will be altered to be half of m57 for crotyl alcohol

Test 3 uses the full 2-CrotAcetExp#2.csv and the original signals
Variable inputs are:
rpcMoleculesToChange = ['CO2']
rpcMoleculesToChangeMF = [[28,44]]
rpcTimeRanges = [[300,600]]
In CO2, the reference signal for m44 should be altered to be about 3 times larger than the signal for m28

Test 4 uses the full 2-CrotAcetExp#2.csv and the original signals
Variable inputs are:
rpcMoleculesToChange = ['Crotyl Alcohol','CO2']
rpcMoleculesToChangeMF = [[57,70],[28,44]]
rpcTimeRanges = [[300,600],[300,600]]
m70 for crotyl alcohol is replaced with the product of the ratio of m70/m57 from 300 to 600 and the signal at m57
Likewise, m44 for CO2 is replaced with the product of the ratio of m44/m28 from 300 to 600 and the signal at m28

Test 5 uses the full 2-CrotAcetExp#2.csv and the original signals
Variable inputs are:
rpcMoleculesToChange = ['Crotyl Alcohol']
rpcMoleculesToChangeMF = [[57,70,44]]  
rpcTimeRanges = [[300,600]]
m70 for crotyl alcohol is replaced with the product of the ratio of m70/m57 from 300 to 600 and the signal at m57. ~1.89 NOTE: The correct number is the ratio of the averages (not the average ratio) see 2-CrotAcetExp#2_300_to_600.xlsx.
m44 for crotyl alcohol is replaced with the product of the ratio of m44/m57 from 300 to 600 and the signal at m57 (results in a large number of ~176)

Test 6 uses the full 2-CrotAcetExp#2.csv and the original signals
Variable inputs are:
rpcMoleculesToChange = ['Crotyl Alcohol','CO2']
rpcMoleculesToChangeMF = [[57,70,44],[28,44]]
rpcTimeRanges = [[300,600],[300,600]]
m70 for crotyl alcohol is replaced with the product of the ratio of m70/m57 from 300 to 600 and the signal at m57
m44 for crotyl alcohol is replaced with the product of the ratio of m44/m57 from 300 to 600 and the signal at m57 (results in a large number)
m44 for CO2 is replaced with the product of the ratio of m44/m28 from 300 to 600 and the signal at m28

