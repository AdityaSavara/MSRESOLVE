This is the unit test for iterative analysis.
test_1.py is the test file.

Basically, a non_iterative analysis is done using test_1_initial_input_noniterative.py

Then, an iterative analysis is done using test_1_initial_input_iterative.py.

This test was carefully crafted, there are some comments inside test_1.py
The test_1 file also tests iterative analysis's compatibility with concentration finder where the known concentrations at certain signals pertain to separate molecules.

In test_1_initial_input_iterative file, the following lines were included to use the concentration finder feature:
#//Concentration Finder//
#this last set of inputs is where you enter your conversion factors from raw signal to concentration, unlike most rows, do not leave brackets around chosen numbers
#here you put in a known raw signal intensity and the known concentration it corresponds to. 
concentrationFinder = 'yes'
TSC_List_Type = 'SeparateMolecularFactors'
moleculesTSC_List = ['Acetaldehyde','CO','Ethanol']
moleculeSignalTSC_List = [1.66945,1.2,1.1]
massNumberTSC_List = [29,28,44]
moleculeConcentrationTSC_List = [0.05,0.01,0.03]	#pressure can also be used in subsitute
unitsTSC = 'bar'	#the units will not be used in calculations so any units may be used

The ResolvedConcentrations.csv file (concentrations from non-iterative run) is compared to the concentrations from iterative runs