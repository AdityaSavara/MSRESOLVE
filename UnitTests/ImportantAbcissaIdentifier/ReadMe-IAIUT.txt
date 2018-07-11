This directory contains the unit tests of the important abscissa identifier function from MSRESOLVESG. 

In order to be Pytest compatible and UnitTesterSG compatible, the user input file and xyyy data functions also are included.  
To run the unit tests separately, execute them out of the Anaconda Prompt, or run from Spyder directly. 

There are 5 unit tests. All of them have an identical format, but they test different aspects of the function.

  Test 1 tests a simple 3 by 4 matrix of reference values, all of which are integers. This test is testing the basic algorithm of the function. 
  Test 2 tests the same matrix, but has a moleculesLikelihood array and a sensitivityValues array. This test is of a subfunction, ListLengthChecker. The two arrays are intentionally of the wrong length, so a written warning should be printed, but no error should be thrown. 
  Test 3 tests a single column of reference data, and should return the largest value in it. 
  Test 4 tests both a negative and a zero in the reference data. Both of these values should be replaced with zeros during the run of the program. 
  Test 5 tests a full array of moleculesLikelihood, so that it will be applied within the function (without alteration to the reference data). 
  
The expected results are saved as pickled files and text files. Retaining these files is necessary as reference outputs for the test when unit testing in future versions.
