README

This directory contains the unit tests of the significanceFactorCheck function from MSRESOLVE.

In order to be Pytest compatible and UnitTesterSG compatible, the user input file and xyyy data functions also are included.  
To run the unit tests separately, execute them out of the Anaconda Prompt, or run from Spyder directly.

The significanceFactorCheck funciton calculates the significance factors for all elements of a reference array. It then sums these elements and stores the largest significance factors in a list of a user specifecd size along with the mass fragment combination the sum belongs to. The last value in the return is a boolean representing if the last iteration was stored in the list.

There are 2 test cases:
	1) Calculates and stores all of the significance factors when a value of 60 is continuously inserted in the second row.
	2) Only stores the top 3 significance factors when 54 is inserted in a diagonal.