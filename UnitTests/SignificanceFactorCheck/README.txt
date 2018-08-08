README

This directory contains the unit tests of the significanceFactorCheck function from MSRESOLVE.

In order to be Pytest compatible and UnitTesterSG compatible, the user input file and xyyy data functions also are included.  
To run the unit tests separately, execute them out of the Anaconda Prompt, or run from Spyder directly.

The significanceFactorCheck funciton calculates the significance factors for all elements of a reference array. It then sums the significance factors across the reference array. The sum is multiplied by a negative so it can be added to a sorted list containing the greatest magnitude of significance sums. An N number of largest magnitude significance sums and stored in order of descending magnitude. The corresponding N mass fragment combinations are stored in a parallel list according to how their significance sums are inserted in list of magnitudes. The last value in the funciton return is a boolean representing if the last iteration was stored in the list.

There are 2 test cases and 4 test files:
	1) Calculates and stores all of the significance factors when a value of 60 is continuously inserted in the second row.
	2) Only stores the top 3 significance factors when 54 is inserted in a diagonal.
	3) Test case one, but uses set results and only compares the mass fragment combinaiton at the top of the best mass fragments list. In this case the fragment is [9,8,7,6,5]
	3) Test case two, but uses set results and only compares the mass fragment combinaiton at the top of the best mass fragments list. In this case the fragment is [9,8,7,6,5]
