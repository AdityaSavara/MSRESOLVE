README

This directory contains the unit tests of the roughUniquenessCheck function from MSRESOLVE.

In order to be Pytest compatible and UnitTesterSG compatible, the user input file and xyyy data functions also are included.  
To run the unit tests separately, execute them out of the Anaconda Prompt, or run from Spyder directly.

The roughUniquenessCheck is a funciton used by the bestMassFragChooser to rank a mass fragment combination by the number of zeros in the reference data array. The lower the the sum of the rowSumsList, the less zeros are present in the data array. It will return a list containing the mass fragment combination with the smallest number of zeros along with the sum of the rowSumsList. The last return value is a boolean that represents if the last mass fragment iterated through was stored in the list.

Two test cases are represented in 4 test files
	1)This is the expected use of the function that stores 3 mass fragment combinations.
	2)All of the row sums are the same. The funciton should only store the 1st 4 mass fragment combinations.
	3)Uses set expected results. The expected result is set as the mass fragment combination with the largest number of zeros. This is [3,4,5,6,7].
	4)Uses set expected results. The expected result if the first mass fragment combinaiton in the list, [1,2,3,4,5] since all of the reference arrays are the same.

The expected results are saved as pickled files and text files. Retaining these files is necessary as reference outputs for the test when unit testing in future versions.