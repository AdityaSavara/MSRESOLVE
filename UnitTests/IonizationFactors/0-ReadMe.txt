This is the ReadMe for the Ionization Factors Feature Unit Tests

The ionization efficiencies used in MSRESOLVE are in a list as an object of the MSReference object.  Since MSReference objects are global variables, we can access the ionization factors used after MSRESOLVE.main() has been executed and use them as the calculated results.

Test 1 uses AcetaldehydeNISTRefKnownFactors.csv
This reference file has values populating the relativeIonizationEfficiencies row.  These are the first values populateIonizationEfficiecies looks for.
So the expected results can simply just be the row of ionization factors from the Reference File.

Test 2 uses AcetaldehydeNISTRefMatchingMolecule.csv
We have two molecules in the reference data that match molecules in the Molecular Ionization Data: Carbon Monoxide and Carbon Dioxide.
Since the program doesn't have an ionization factor given to it, it will look in the MID Dictionary for a molecule that matches next.  If one is found (in our case we have two) then the average of the RS_Values stored in the MID dictionary is used as the ionization factor.
From the MI data, we know Carbon Monoxide has an average RS_Value of 1.05 and Carbon Dioxide has an average RS_Value of 1.4.  So we just overwrite the ionizationFactorsRN2 with these values at the proper indicies.
Expected results is a sliced array containing only the ionizatoin factors for CO and CO2.  The calculated results is the sliced array of the MSReference variable ionizationEfficienciesList containing only the ionization efficiencies for CO and CO2.

Test 3 uses AcetaldehydeNISTRefKnownTypes.csv
The relativeIonizationEfficiencies are all set to 'unknown', and the molecule types for CO and CO2 have been populated.  The only type used in the reference file is Main Group.
The excel file LinearFits.xlsx has the required data to calculate slopes and intercepts for Main Group molecules.  The slopes and intercepts from the linear fit data are copied into the test_3.py file as lists (i.e. typeCoefficients = [slope, intercept] ).
A poly1d object is made for each ionization type's polynomial coefficients.  An array the same length as the number of electron numbers is initialized as a row of zeros.  This array will store the ionization factors.  Using a for loop to iterate over this array, we can evaluate the ionization factor using the poly1d object and the electron number of the molecule with polyval.
The expected results is just simply the populated ionization array.  The objects match within the defined tolerances but the strings do not due to rounding.

Test 4 uses AcetaldehydeNISTRefDefault.csv
This reference pattern is the same as AcetaldehydeNISTRefMixed2.tsv from the main directory.
Using this reference pattern should default populateIonizationEfficiencies to using the Madix and Ko equation.  An array of zeros having the same length as ElectronNumbers is initialized.  It is populated using a for loop and calculating the ionization factor using the Madix and Ko equation.