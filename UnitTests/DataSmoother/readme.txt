Data Smoother Unit Tests

DataSmoother takes in 8 variables and outputs smootedData.
The 8 variables consists of:
data which is just an array of data
abscissa which is an array of values of time or temperature corresponding to the values in the same index in data
headers which is a list containing the names of the data in the columns of the same index
dataSmootherChoice is a string (either 'timerange' or 'pointrange') that lets the user decide which method to use.
dataSmootherTimeRadius is an int defining the increment of time that dataSmoother uses to collect and fit data points across
dataSmootherPointRadius is an int defining the number of points before and after the current point that are used to fit data.
*NOTE* dataSmootherPointRadius is not actually used in the function, dataSmootherTimeRadius is so the user should just use dataSmootherTimeRadius for both
headersToConfineTo is a list containing the header names of the columns the user wants to change.  If left empty all the columns will change.
polynomialOrder is just the order the user wants to use to fit the data.

All these tests use headers = [34,35], data as an array of two columns and seven rows, and abscissa being a one-dimensional array of seven values.
Test1 uses timerange and time radius set to 2 with polynomial order of 2.  So we expect dataSmoother to fit points in data to a 2nd degree polynomial using the two times before and the two times after the data point being changed.

Test2 sets headersToConfineTo to [34] so we expect only the first column to change.

Test3 sets dataSmootherChoice to 'pointrange' and uses point radius (time radius) to be 2.  headersToConfineTo is left blank.  Abscissa is changed to np.array([1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]) to test pointrange
Here we expect dataSmoother to use the two points before and the two points after the point being changed to be fit with a 2nd order polynomial.



