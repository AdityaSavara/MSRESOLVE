This directory contains unit tests that are works in progress. They are not called by pytestDriver in their current form. 

180807 Charles

I have added 2 tests in progress:testKeepOnlyYYYYCol.py and TestRemoveSignals.py. 
testKeepOnlyYYYYCol.py tests the KeepOnlySelectedYYYYColumns function from XYYYDataFunctionsSG.
TestRemoveSignals.py tests the RemoveSignals function from XYYYDataFunctionsSG.

Both of these tests currently are written to print the results to the console and do not yet use the UnitTester module. 

In order to make them operational, add functions to call UnitTester on the results, save the output, and rename them to be accesible by pytest.
