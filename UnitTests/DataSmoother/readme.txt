KeepOnlySelectedYYYYColumns Unit Tests

This function takes in 3 inputs (2 arrays and a list) and outputs 2 arrays.  The YYYYdata variable just contains an array of values.
The headerValues variable is a one-dimensional array of numbers (as strings) that are the "titles" to the values in data.
headerValuesToKeep is the list containing the titles of the headers (as strings) that you want to keep.

Test1 is rather easy to visualize.  The YYYYdata is a numpy array of 14 columns and 5 rows.  The headersValue variable is an array of 14 columns and 1 row.  
These values are the headers of the columns in YYYYdata that have the same index. (i.e. '2' maps to the first column, '18' maps to the second column, etc.)
headerValuesToKeep is just ['2','18'] which tells KeepOnlySelectedYYYYColumns to keep the columns with the header 2 and 18.
In this case that happens to be the first two columns in YYYYData so the output will be two arrays: one with the first two columns, and one with the headers to keep.

Test2 uses the same YYYYdata as test1.  The headerValues variable is an array with 14 columns and 1 row. The headerValuesToKeep variable is a list of 10 values.
We expect KeepOnlySelectedYYYYColumns to output an array of YYYYdata that has 10 columns and 5 row where each column is a column from YYYYdata that has a value in both headerValues and headerValuesToKeep.
The other array that KeepOnlySelectedYYYYColumns to output should be a one dimensional array of 10 values where each value is the header to the columns in the new YYYYdata.

Both tests gave expected outputs.