KeepOnlySelectedYYYYColumns Unit Tests

This function takes in 3 inputs (2 arrays and a list) and outputs 2 arrays.  The YYYYdata variable just contains an array of values.
The headerValues variable is an array of numbers (as strings) that are the "titles" to the values in data.
headerValuesToKeep is the list containing the titles of the headers (as strings) that you want to keep.

Test1 uses a YYYYdata as a numpy array of 14 columns and 5 rows.  The headerValues variable is an array with 14 columns and 1 row. The headerValuesToKeep variable is a list of 10 values.
We expect KeepOnlySelectedYYYYColumns to output an array of 10 columns and 5 rows.  The 10 columns will be the 10 columns in YYYYdata that correspond to the same index in headerValues for those values in headerValuesToKeep

Test2 is a little easier to visualize.  The same YYYYdata and headerValues is used from test1.  headerValuesToKeep has been changed to ['2','18'] which tells KeepOnlySelectedYYYYColumns to keep the columns with the header 2 and 18.
In this case that happens to be the first two columns in YYYYData so the output will be two arrays: one with the first two columns, and one with the headers to keep.

Both tests gave expected outputs.