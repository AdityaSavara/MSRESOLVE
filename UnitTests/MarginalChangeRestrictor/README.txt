This directory contains the unit tests of the marginal change restrictor function from MSRESOLVE.

In order to be Pytest compatible and UnitTesterSG compatible, the user input file and xyyy data functions also are included.  
To run the unit tests separately, execute them out of the Anaconda Prompt, or run from Spyder directly.

There are 4 unit tests designed to test different aspects of the Marginal Change Restrictor.

	Test 1 tests an abscissa of 2 values and a 2 by 3 working_data matrix with a zero in the first row. This test is testing the insertion of a new row following a zero containing half of the ignorable delta Y threshold and the insertion of new rows if they exceed the marginal change restriction.
	Test 2 tests an abscissa of 2 values and a 2 by 3 working_data matrix with two sign changes.This test is testing the insertion of two new rows following a sign change containing a 0, the insertion of new rows surrounding a zero containing half of the ignorable delta Y threshold and the insertion of new rows if they exceed the marginal change restriction.
	Test 3 tests the use of the Superfluous Row Cleanup in the data points created in test 2 since test 2 inserted 2 unnecessary rows.
	Test 4 tests the Interpolate Accompanying Arrays function with the abscissa created in test 2 and a 2 x 7 example concentrationBounds matrix.

The expected results are saved as pickled files and text files. Retaining these files is necessary as reference outputs for the test when unit testing in future versions.
