This is the ReadMe for the Ionization Factors Feature Unit Tests

The ionization efficiencies used in MSRESOLVE are in a list as an object of the MSReference object.  Since MSReference objects are global variables, we can access the ionization factors used after MSRESOLVE.main() has been executed and use them as the calculated results.

Test 1 uses AcetaldehydeNISTRefKnownFactors.csv
This reference file has values populating the knownIonizationFactorsRelativeToN2 row.  These are the first values populateIonizationEfficiecies looks for.
So the expected results can simply just be the row of ionization factors from the Reference File.

Test 2 uses AcetaldehydeNISTRefMatchingMolecule.csv
We have two molecules in the reference data that match molecules in the Molecular Ionization Data: Acetaldehyde and Ethanol.
This reference file has the same knownIonizationFactorsRelativeToN2 as AcetaldehydeNISTRefKnownFactors.csv (from test 1) with the exception of Acetaldehyde and Ethanol which are set to 'unknown'.
Since the program doesn't have an ionization factor given to it, it will first look in the MID Dictionary for a molecule that matches.  If one is found (in our case we have two) then the average of the RS_Values stored in the dictionary is used as the ionization factor.
Expected results are set as the knownIonizationFactorsRelativeToN2 (same as test 1) but, from the MI data, we know Ethanol has an average RS_Value of 3.25 and Acetaldehyde has an RS_Value of 2.6.  So we just overwrite the knownIonizationFactors with these values at the proper indicies.

Test 3 uses AcetaldehydeNISTRefKnownTypes.csv
The molecules Acetaldehyde and Ethanol have been renamed in this reference file to Ethenal and EtOH.  This is so the molecule names do not match any molecule name in the MID Dictionary.  The knownIonizationFactorsRelativeToN2 are all set to 'unknown', and the molecule types have been populated.  We only have four ionization types: Aldehydes, Alcohols, Hydrogen non-metal-ides, and carbon containing nonmetals.
The excel file LinearFits.xlsx has the required data for the four types of molecules we have.  The slopes and intercepts from the linear fit data are copied into the test_3.py file as numpy arrays (i.e. typeCoefficients = numpy.array([slope, intercept]) ).
A poly1d object is made for each ionization type's polynomial coefficients.  An array the same length as the number of electron numbers is initialized as a row of zeros.  This array will store the ionization factors.  Using a for loop to iterate over this array, we can evaluate the ionization factor using the poly1d object and the electron number of the molecule with polyval.
The expected results is just simply the populated ionization array.  The objects match within the defined tolerances but the strings do not due to rounding.

Test 4 uses AcetaldehydeNISTRefDefault.csv
This reference pattern is the same as AcetaldehydeNISTRefMixed2.csv from the main directory but Acetaldehyde and Ethanol have been renamed to Ethenal and EtOH so they are not overwritten with the average RS_Value from the MID Dictionary.
Using this reference pattern should default populateIonizationEfficiencies to using the Madix and Ko equation.  An array of zeros having the same length as ElectronNumbers is initialized.  It is populated using a for loop and calculating the ionization factor using the Madix and Ko equation.