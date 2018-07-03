Data Smoother Unit Tests

DataSmoother takes in 8 variables and outputs smootedData.
The 8 variables consists of:
data which is just an array of data
abscissa which is a one-dimensional array of values of time or temperature corresponding to the values in the same index in data
headers which is a list containing the names of the data in the columns of the same index
dataSmootherChoice is a string (either 'timerange' or 'pointrange') that lets the user decide which method to use.
dataSmootherTimeRadius is an int defining the increment of time that dataSmoother uses to determine which points to use in fitting data points across
dataSmootherPointRadius is an int defining the number of points before and after the current point that are used to fit data.
headersToConfineTo is a list containing the header names of the columns the user wants to change.  If left empty all the columns will change.
polynomialOrder is just the order the user wants to use to fit the data.

*NOTE* dataSmootherPointRadius is not actually used in the function.  Rather dataSmootherTimeRadius is used in both scenarios.  So if you are using pointrange, the radius defined by the user should be placed in the dataSmootherTimeRadius variable.

For the test1 and test2, the data is an array of two columns and seven rows, and abscissa being a one-dimensional array of seven values.
Test1 uses timerange and time radius set to 2 with polynomial order of 2.  So we expect dataSmoother to fit points in data to a 2nd degree polynomial using the two times before and the two times after the data point being changed.

Test2 sets headersToConfineTo to [34] so we expect only the first column to change.

