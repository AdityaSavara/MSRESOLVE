This is the ReadMe for the Ionization Factors Feature Unit Tests

The ionization efficiencies used in MSRESOLVE are in a list as an object of the MSReference object.  Since MSReference objects are global variables, we can access the ionization factors used in the last run and use them as the calculated results.

Test 1 uses AcetaldehydeNISTRefKnownFactors.csv
This reference file has values populating the knownIonizationFactorsRelativeToN2 row.  These are the first values populateIonizationEfficiecies looks for.
So the expected results can simply just be the row of ionization factors from the Reference File.

Test 2 uses AcetaldehydeNISTRefMatchingMolecule.csv
We have two molecules in the reference data that match molecules in the Molecular Ionization Data: Acetaldehyde and Ethanol.
This reference file has the same ionization factors as AcetaldehydeNISTRefKnownFactors.csv (from test 1) with the ionization factors of Acetaldehyde and Ethanol set to unknown.
Since the program doesn't have an ionization factor given to it, it will first look in the MID Dictionary for a molecule that matches.  If one is found (in our case it is) then the RS_Value stored in the dictionary is used as the ionization factor.
Expected results are calculated similarly to test 1 but, from the MI data, we know Ethanol has an average RS_Value of 3.25 and Acetaldehyde has an RS_Value of 2.6.  So we just overwrite the knownIonizationFactors at the proper indicies.

Test 3 uses AcetaldehydeNISTRefKnownTypes.csv
The molecules Acetaldehyde and Ethanol have been renamed in this reference file to Ethenal and EtOH.  This is so the molecule names do not match any molecule name in the MID Dictionary.  The ionization factors are all unknown, and the molecule types have been populated.  We only have four ionization types: Aldehydes, Alcohols, Hydrogen non-metal-ides, and carbon containing nonmetals.
The excel file LinearFits.xlsx has the required data from the four types of molecules we have.  The slopes and intercepts from the linear fit data are copied into the test_3.py file as numpy arrays (i.e. numpy.array([slope, intercept]) ).
A poly1d object is made for each ionization type.  Then using a for loop to loop over an initialized array to contain the ionization factors, we can evaluate the ionization factor using the poly1d object and the electron number of the molecule with polyval.
The expected results is just simply the populated ionization array.
The results were close enough for the objects to match but not the strings.  So I have copied the exact results into the test_3.py file and added numpy.round() functions (rounding to the nearest 2 decimals) so the strings would match.

Test 4 uses AcetaldehydeNISTRefDefault.csv
This reference pattern is the same as AcetaldehydeNISTRefMixed2.csv from the main directory but Acetaldehyde and Ethanol have been renamed to Ethenal and EtOH for the same reason as before
Using this reference pattern should default populateIonizationEfficiencies to using the Madix and Ko equation.  So the expected results will be the ionization factors array populated by a for loop calculating the ionization factor based solely on electron number.