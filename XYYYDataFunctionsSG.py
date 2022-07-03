"""
Functions for use with MS data. All of these functions assume that the list of Y data lists 
have been turned into numpy arrays and transposed. Then the data series for each mass is a column 
in the new numpy array and the number of columns is the same as the number of elements of 
'choosenMassFragments' in UserInput.py.
"""

import numpy
import math
import copy
import pandas 


#Function to retrieve Y values and Y Value uncertainties given  A, B, C coefficients and covMatrixOfParameters
#list of Parameters is "A,B,C"
def returnPolyvalEstimatesAndUncertainties(x_values, abcCoefficients, abcCoefficients_covMat):
    n = len(abcCoefficients) - 1
    x_values = numpy.array(x_values)
    TT = numpy.vstack([x_values**(n-i) for i in range(n+1)]).T
    y_predicted = numpy.dot(TT, abcCoefficients)  # matrix multiplication calculates the polynomial values
    Cov_y = numpy.dot(TT, numpy.dot(abcCoefficients_covMat, TT.T)) # Cov_y = TT*C_z*TT.T
    y_predicted_uncertainties = numpy.sqrt(numpy.diag(Cov_y))  # Standard deviations are sqrt of diagonal
    return y_predicted, y_predicted_uncertainties


#Function to retrieve Y values uncertainties given  A, B, C coefficients and covMatrixOfParameters
#list of Parameters is "A,B,C"
def returnPolyvalEstimatedUncertainties(arrayOfAbscissaValues, listOfParameters, covMatrixOfParameters):
    n = len(listOfParameters) - 1
    TT = numpy.vstack([arrayOfAbscissaValues**(n-i) for i in range(n+1)]).T
    yi = numpy.dot(TT, listOfParameters)  # matrix multiplication calculates the polynomial values
    Cov_y = numpy.dot(TT, numpy.dot(covMatrixOfParameters, TT.T)) # Cov_y = TT*C_z*TT.T
    uncertaintiesArray = numpy.sqrt(numpy.diag(Cov_y))  # Standard deviations are sqrt of diagonal
    return uncertaintiesArray

    #Takes a 1D array or list and returns a comma separated string.
def arrayLikeToCSVstring(inputArray):
    #First check if the objects are strings. If they are, we will have to remove single quotes from our final string.
    if (type(inputArray[0]) == type("string") ) or (type(inputArray[0]) == type(numpy.str_("string"))) :
        stringObjects = True
    else:
        stringObjects = False
    listString = str(list(inputArray))
    CSVstring = listString[1:-1]
    if stringObjects == True:
        CSVstring= CSVstring.replace("'","")
    return CSVstring
def AppendColumnsToCSV(CSVName, YYYYData, columnheaders, rowIndex = [], rowIndexHeader = []):
    #rowIndex and rowIndexHeader are only used if there is not already a file with the specified name. 
   
    #create pandas data frame of new data
    newColumns = pandas.DataFrame(data = YYYYData, columns = columnheaders)
    
    #make sure that the solved concentrations file exists
    try:
        #read in the previously solved data    
        if '.csv' in CSVName:
            totalColumns = pandas.read_csv(CSVName)
        if '.tsv' in CSVName:
            try:
                totalColumns = pandas.read_csv(CSVName, sep = '\t', encoding='utf8')
            except:
                totalColumns = pandas.read_csv(CSVName, sep = '\t', encoding='utf16')
    #If the file doesn't exist
    except FileNotFoundError:
        #then the times need to be included in this writing
        index = pandas.DataFrame(data = rowIndex, columns = rowIndexHeader)
        totalColumns = pandas.concat((index,newColumns), axis = 1)
    #Otherwise the new data is simply added to the old
    else:
        totalColumns = pandas.concat((totalColumns,newColumns), axis = 1)
    #all the data is rewritten to the csv file
    if '.csv' in CSVName:
        totalColumns.to_csv(CSVName, index = False, sep =',') 
    if '.tsv' in CSVName:
        totalColumns.to_csv(CSVName, index = False, sep ='\t')     
    return None

def TrimReferenceFileByMolecules(moleculesToSave, referenceFileName):
    
    #it is useful to trim whitespace from each string. This prevents future errors in comparison. 
    for moleculeIndex in range(len(moleculesToSave)):
        moleculesToSave[moleculeIndex] = moleculesToSave[moleculeIndex].strip()
        
    #generate the dataframe
    referenceFile = pandas.read_csv(referenceFileName, header = None)
    #create an array to store indicies of previously used molecules 
    moleculesToDelete = []
    #for each column that isn't in the mass fragment abcsissa
    for moleculeCount in range(len(referenceFile.columns)-1):
        #if the name of that molecule is one to be saved
        if not referenceFile.iloc[1,1+moleculeCount].strip() in moleculesToSave:
            #store the index of that column
            moleculesToDelete.append(1+moleculeCount)
    #delete all the unspecified molecules 
    referenceFile = referenceFile.drop(moleculesToDelete, axis = 1)
#    if unusedReferenceFileName != None:
#        referenceFile.to_csv(unusedReferenceFileName, header = False, index = False)      
#    else: return referenceFile
    return referenceFile

'''
MSDataWriterXYYY() replaces ExportXYYYData() for writing
the export data from the MSData class member function ExportMSData().

filename: type string, its the desired filename
data: 2-d numpy array, shape (# samples in series/len(times), # masses investigated)
abscissa: list of times/temps/etc. when data was taken
dataHeader: list of the masses (e.g. [m2, m34, m44, ...])
'''
def MSDataWriterXYYY(filename, data, abscissa, dataHeader, abscissaHeader):

    # the dataHeader is a list of column names for data
    # it needs to include an entry for abscissa as well
    dataHeader.insert(0,abscissaHeader)

    # for numpy.savetxt() we need the header as a string
    # create a string of the column headers seperated by commas
    dataHeaderAsString = ','.join(map(str,dataHeader))

    # add the abscissa column into the data array
    data = numpy.column_stack((abscissa, data))
    
    # save data to the file
    numpy.savetxt(filename, data, delimiter=',', header=dataHeaderAsString, comments='')
    


'''
The DataSmootherPolynomialSmoothing function takes the currentWindow 
(2-d numpy array) and timelist (1-d list) that are generated with 1 of the 4 options
in DataSmoother (using ExtractWindow...Radius() function) and centered about 'currentTime'. 
Then for a single currentTime or for every element of 'timeslist' it
takes the relevant data and  performs the smoothing and then returns the smoothedData.
'''
def DataSmootherPolynomialSmoothing(timeslist, currentWindowAsArray,  currentTime, polynomialOrder=1, returnAllTimes = False):
    polyFitObject = numpy.polyfit(timeslist,currentWindowAsArray,polynomialOrder)
    if returnAllTimes == False:
        smoothedPoints = numpy.polyval(polyFitObject, currentTime)
    if returnAllTimes == True:
        smoothedPoints = currentWindowAsArray*0.0 #just getting the same size as the data window.
        for timePointIndex in range(len(timeslist)):
            smoothedPoints[timePointIndex] = numpy.polyval(polyFitObject, timeslist[timePointIndex])                
    return  smoothedPoints


'''
GetTimeRadiusIndexLimits() function:
Here we will search for the appropriate abscissa index limits such that the
the given TimeRadius is satisfied. We start the search at the current time and
move outwards in both directions simultaneously
returns tuple (limitIndexLimit, lowerIndexLimit)
the index limits are inclusive
(e.g. for abscissa[currentTimeIndex]=2.5 and dataSmootherTimeRadius=1
if time 1.5 exists in abscissa then its index will be returned, but no lower)
'''
def GetTimeRadiusIndexLimits(abscissa, dataSmootherTimeRadius, currentTimeIndex):


    upperIndexLimit = currentTimeIndex
    lowerIndexLimit = currentTimeIndex

    currentTime = abscissa[currentTimeIndex]
    upperTimeLimit = currentTime + dataSmootherTimeRadius
    lowerTimeLimit = currentTime - dataSmootherTimeRadius

    upperLimitFound = False
    lowerLimitFound = False

    # first make sure that the TimeLimits are not beyond the range of 'abscissa'
    # if so then set the indexLimits to the max/min index of the abscissa
    if (currentTime + dataSmootherTimeRadius >= abscissa[-1]):
        upperLimitFound = True
        upperIndexLimit = len(abscissa) - 1 # highest abscissa index

    if (currentTime - dataSmootherTimeRadius <= abscissa[0]):
        lowerLimitFound = True
        lowerIndexLimit = 0 # lowest abscissa index

    # normally the time range/radius will fall within the abscissa range
    while ((not upperLimitFound) or (not lowerLimitFound)):
        # check/increment  upper
        if (not upperLimitFound):
            if (abscissa[upperIndexLimit + 1] > upperTimeLimit):
                # if true then the time has finally exceeded the limit -> stop inrementing
                upperLimitFound = True
            elif(abscissa[upperIndexLimit + 1] <= upperTimeLimit):
                # still under the limit, increment the index, go to later time
                upperIndexLimit += 1

        # check/decrement  lower
        if (not lowerLimitFound):
            if (abscissa[lowerIndexLimit - 1] < lowerTimeLimit):
                # if true then the time has finally fallen below the limit -> stop decrementing
                # but abscissa[lowerIndexLimit] is still <= lowerTimeLimit
                lowerLimitFound = True
            elif(abscissa[lowerIndexLimit - 1] >= lowerTimeLimit):
                # still above the limit, dencrement the index, go to earlier time
                lowerIndexLimit -= 1

    return (lowerIndexLimit, upperIndexLimit)

'''
Returns 'dataWindowsAsTuples' a list of tuples of the appropriate time and
paired data windows given a choice of radius and radius type (pointrange or timerange).
The tuples in the 'dataWindowsAsTuples' are of the form (timeslist, currentWindowAsArray), there is one
tuple (i.e. element of 'dataWindowsAsTuples') for every time in abscissa. 
 -timeslist is a list of times that is the appropriate subset of abscissa given the radius type and radius.
 -currentWindowAsArray is a 2d numpy array of dimensions [len(abscissa), data.shape[1]], i.e.
  there is a column for every mass and each column is of the same length as the 
  (each data point in the columns is paired to a time)
'''
def GetDataWindows(data, abscissa, radius, dataSmootherChoice):
    dataWindowsXvaluesInArrays = []
    dataWindowsYYYYvaluesInArrays = []
    
    # for every time/element in abscissa populate 'dataWindowsAsTuples'
    for timecounter,timeValue in enumerate(abscissa):
        # if we use a timeradius we need to find the appropriate indices
        # from abscissa that relate to the time window
        if dataSmootherChoice == 'timerange':
            indexLimits = GetTimeRadiusIndexLimits(abscissa, radius, timecounter)

        # if we use pointradius its easier to get the index limits
        elif dataSmootherChoice == 'pointrange':
            # just move radius from the current time index
            # make sure not to go beyond 0 or len(abscissa + 1)
            indexLimits = (max(timecounter - radius, 0) ,
                           min(timecounter + radius,len(abscissa) + 1))

        # now just slice data and abscissa arrays
        # appropriately to form the windows for this timecounter
        timeslist = abscissa[indexLimits[0]:indexLimits[1]]
        currentWindowAsArray= data[indexLimits[0]:indexLimits[1], :]

        # Add to the list
        dataWindowsXvaluesInArrays.append(timeslist)
        dataWindowsYYYYvaluesInArrays.append(currentWindowAsArray)

    return (dataWindowsXvaluesInArrays, dataWindowsYYYYvaluesInArrays)

'''
GetRelevantData() is used when headersToConfineTo != [] (only needed for options 3 and 4). 
In other words when we only want to smooth the data associated with certain masses.
'headers' is a list of all masses that are being analyzed (should have the same number 
of elements as there are Y in XYYY data). headersToConfineTo
is a list of masses that we want to smooth, it is a subset of headers. This function 
takes 'data' and returns a numpy array of the same form 'tempData'
that has only the columns of data series to be smoothed.
'''
def GetRelevantData(data,abscissa, headers, headersToConfineTo):
    # Determine which columns of 'data' we will be smoothing
    # by populating 'extractedColumnTracker' '0' for an irrelevant column
    # '1' for relevant columns
    extractedColumnTracker = numpy.zeros(len(headers))
    for relevantHeaderIndex,relevantHeader in enumerate(headersToConfineTo):
        if relevantHeader in headers: # just make sure it exists
            extractedColumnTracker[list(headers).index(relevantHeader)] = 1.0 # this assumes no
                                                                        # duplicates in 'headers' !!
        else:
            print("You choose to smooth mass data which doesnt exist in "+
                  "'choosenMassFragments'")

    # create extractedData, same form as data but loses the columns that we dont want to smooth
    extractedData = numpy.zeros( [len(data[:,0]), len(headersToConfineTo)])
    tempIndex = 0
    for massIndex,relevance in enumerate(extractedColumnTracker):
        if relevance:
            extractedData[:,tempIndex] = data[:,massIndex]
            tempIndex += 1

    return (extractedData, extractedColumnTracker)



"""
ReconstructSmoothedData is used in DataSmoother when not all mass 
data series are to be smoothed.It takes the 
'smoothedExtractedData' and 'extractedColumnTracker' and uses them
to reconstruct the full 'smoothedData' array. smoothedExtractedData is 
a data array that only has mass data columns corresponding
to masses in 'headersToConfineTo'. 'extractedColumnTracker' is a 
list of 1 and 0 the same length as headers, when there is a 1 it means 
that corresponding mass data should be smoothed, 0 => that the data 
for that mass should not be smoothed. 'unSmoothedData' is a copy of 
the 'data' array that is passed into the main DataSmoother() function.
'ReconstructSmoothedData' is basically
the inverse of GetRelevantData() and allows reconstruction  of a full 
smoothedData array once the relevant mass data series (indicated 
by 'headersToConfineTo' and 'extractedColumnTracker') have been smoothed.
"""
def ReconstructSmoothedData(unsmoothedData, smoothedExtractedData, extractedColumnTracker):
    # Now construct smoothedData from smoothedExtractedData, data and extractedColumnTracker
    # i.e. combine the columns we just smoothed with the columns that didnt need smoothing
    # in the proper order
    tempIndex = 0
    for columnIndex,relevance in enumerate(extractedColumnTracker):
        if relevance: # should be equal number of relevant columns ('1' in relevantColumnts)
                      # and number of columns to replace
            unsmoothedData[:,columnIndex] = smoothedExtractedData[:,tempIndex]
            tempIndex += 1

    # Now that it is smoothed copy over to 'smoothedData'
    smoothedData = copy.deepcopy(unsmoothedData)

    return smoothedData


'''
KeepOnlySelectedYYYYColumns() compares the values of DataHeaders and HeaderValuesToKeep.
Values(numbers) that are in DataHeaders and not in HeaderValuesToKeep are deleted
from DataHeaders. Further the columns of YYYYData correspond to the DataHeaders array,
i.e. shape(YYYData)[1] = len(DataHeaders). When a value is removed from 
DataHeaders the corresponding column is removed from YYYYData. 
Parameters:
YYYYData- A 2-d numpy array, shape(YYYYData) = (*, len(DataHeaders)), i.e. the columns
          of YYYYData correspond to the entries in DataHeaders
DataHeaders- List of floats/strings/etc. (or 1-d numpy array)
HeaderValuesToKeep- List of floats/strings/etc. NOTE: HeaderValuesToKeep is not necessarily 
                      smaller than DataHeaders, it may contain values not included
                      in DataHeaders.
header_dtype_casting is the type the elements will be cast into during the comparison (int, float, or str). If None, no casting occurs
'''
def KeepOnlySelectedYYYYColumns(YYYYData, DataHeaders, HeaderValuesToKeep, Array1D = False, header_dtype_casting=None):
    import ParsingFunctions as parse
    #Initialize a copy of data headers and header values to keep as empty lists
    copyOfDataHeaders = []
    copyOfHeaderValuesToKeep = []    
    if header_dtype_casting != None: #user is using header_dtype_casting argument
        for elem in DataHeaders: #Loop through data headers
            castElem = header_dtype_casting(elem) #cast elem into the proper type
            if header_dtype_casting == str: #If the preferred type is a string, standardize it
                castElem = parse.standardizeString(castElem)
            copyOfDataHeaders.append(castElem) #append to copied list with the casted elem
        for elem in HeaderValuesToKeep: #Loop through header values to keep
            castElem = header_dtype_casting(elem) #cast elem into the proper type
            if header_dtype_casting == str: #If the preferred type is a string, standardize it
                castElem = parse.standardizeString(castElem)
            copyOfHeaderValuesToKeep.append(castElem) #append to copied list with the casted elem
    else: #If left as default, make copies a direct copy of DataHeaders and HeaderValuesToKeep
        copyOfDataHeaders = copy.copy(DataHeaders)
        copyOfHeaderValuesToKeep = copy.copy(HeaderValuesToKeep)
        
    # list to track indices that should be deleted
    deletion_indices = []

    # loop through DataHeaders, record the indices
    # of values not also found in AbscissaValuesToKeep
    # for deletion
    for (valueIndex,value) in enumerate(copyOfDataHeaders):
        if value not in copyOfHeaderValuesToKeep:
            deletion_indices.append(valueIndex)

    # Now remove the unwanted values from DataHeaders
    # and the corresponding unwanted columns from YYYYData
    DataHeaders = numpy.delete(DataHeaders,deletion_indices)
    axis = 1 #This is the appropriate value for a 2D array of YYYYData
    if Array1D:
        axis = 0 #This is the appropriate value for a 1D array of YData
    #remove the data columns
    if type(YYYYData)!= type(None): #if somehow a Nonetype got in, we will do nothing. (Though callling function is still useful to trim headers.)
        YYYYData = numpy.delete(YYYYData,deletion_indices, axis)

    return (YYYYData, DataHeaders)

#TODO: make a function KeepOnlyYYYYRows() that is very similar to this one.
# It may be useful and could replace some of the functionality of ArrayBuilder()
# and UnecessaryMoleculesDeleter()
    

'''
UncertaintiesFromLocalWindows uses a particular timeRadius to determine the uncertainties around that datapoint.
This Function borrows code flow and sub-functions that are used in the DataSmoother method.
headersToConfineTo is not intended to be used, but is kept to have a parallel arguments list to the DataSmoother function.
At present it is assumed that polynomialOrder will always be kept at 1, but it's conceivable somebody might want to do differently.
We keep this function separate from DataSmoother for several reasons: a) this can run when that doesn't, b) this one will always do all masses, c) this one runs before interpolator on purpose.
'''
def UncertaintiesFromLocalWindows(data,abscissa,headers,UncertaintiesWindowsChoice='pointrange',UncertaintiesWindowsTimeRadius=0,UncertaintiesWindowsPointRadius=5,headersToConfineTo=[],polynomialOrder = 1, UncertaintiesType="standardError"):    
    #UncertaintiesWindowsChoice can be 'pointrange' or 'timerange'. The corresponding radius is used.
    #uncertaintiesType can be standardError or aggregateError or standardDeviation
    UncertaintiesFromData = copy.deepcopy(data) #First make a copy of the data to have the right shape. THis does not have the times, just the intensities.  
    AverageResidualsFromData = UncertaintiesFromData*0.0 #just making an array of the same size.
    if UncertaintiesWindowsChoice == 'timerange':
        UncertaintiesWindowsRadius = UncertaintiesWindowsTimeRadius
    if UncertaintiesWindowsChoice == 'pointrange':    
        UncertaintiesWindowsRadius = UncertaintiesWindowsPointRadius
    
    dataSmootherChoice = UncertaintiesWindowsChoice
    smoothedData = copy.deepcopy(data) #Still need to do smoothing in order to accomplish this.
    if UncertaintiesWindowsChoice == 'timerange':
        dataSmootherRadius = UncertaintiesWindowsTimeRadius
    if UncertaintiesWindowsChoice == 'pointrange':    
        dataSmootherRadius = UncertaintiesWindowsPointRadius

    if headersToConfineTo == []:
        # Get a list of time and data windows for each time in the abscissa
        # these windows are used to perform the smoothing
        (dataWindowsXvaluesInArrays, dataWindowsYYYYvaluesInArrays) = GetDataWindows(data,abscissa, dataSmootherRadius, dataSmootherChoice)
        # replace data points with their smoothed counterparts creating 'smoothedData'
        for timecounter, timeslist in enumerate(dataWindowsXvaluesInArrays):            
            currentWindowAsArray = dataWindowsYYYYvaluesInArrays[timecounter]
            # Smooth the data
            smoothedData[timecounter,:] = DataSmootherPolynomialSmoothing(timeslist, currentWindowAsArray, abscissa[timecounter], polynomialOrder)
            smoothedDataFromWindowOfTimes = DataSmootherPolynomialSmoothing(timeslist, currentWindowAsArray, abscissa[timecounter], polynomialOrder, returnAllTimes = True)
            #Below is Calculating the standardError.
            subtractedDataForTimeCounter =  dataWindowsYYYYvaluesInArrays[timecounter] - smoothedDataFromWindowOfTimes #The first index (rows) are time, second are masses. We want to find the standard deviation across masses.
            averageResidualForTimeCounter = numpy.mean(abs(subtractedDataForTimeCounter), axis=0) #averaging across all times for each mass.
            standardDeviationForTimeCounter = numpy.std((subtractedDataForTimeCounter), axis=0) #standard deviation across all times in the window for each mass.
            standardErrorOfTheMeanForTimeCounter = standardDeviationForTimeCounter/(len(subtractedDataForTimeCounter)**0.5)
            AverageResidualsFromData[timecounter] = averageResidualForTimeCounter
            #Below is Calculating the aggregateError.
            #First get the residuals for this timeCounter.
            numOfValues = len(subtractedDataForTimeCounter[:])
            #now, because unfortunately we are working with windows of data, we need to assess if the number of points is odd or even.
            if (numOfValues % 2) == 0: #then it is even. This is the expected case, and if there are 10 points we need index 5 (just past halfway mark).
                indexToExtract = int(numOfValues/2)
            if (numOfValues % 2) != 0: #then it is odd. This is not the expected case, and if there are 9 points we need index 4 (which is the middle point)
                indexToExtract = numOfValues/2.0
                from math import floor
                indexToExtract = int(floor(indexToExtract))
            residualForThisTimeCounter = abs(subtractedDataForTimeCounter[indexToExtract])            
            aggregateUncertaintyForThisTimeCounter = (residualForThisTimeCounter**2 + averageResidualForTimeCounter**2 + standardErrorOfTheMeanForTimeCounter**2)**0.5
            if UncertaintiesType == "standardError":
                UncertaintiesFromData[timecounter] = standardErrorOfTheMeanForTimeCounter
            if UncertaintiesType == "aggregateError": 
                UncertaintiesFromData[timecounter] = aggregateUncertaintyForThisTimeCounter
            if UncertaintiesType == "standardDeviation": 
                UncertaintiesFromData[timecounter] = standardDeviationForTimeCounter
    return UncertaintiesFromData, AverageResidualsFromData        

'''
The DataSmoothing function 'smooths' data over a certain time or datapoint ranges:
it goes through each mass fragment at a certain time and determines a polynomial that modeles the datapoints around
that mass fragment (within point or time radius). Then applies this determined polynomial to recalculate the value of the datapoint.
after all datapoint of a certain time are analyze the function then resets on the next datapoint. 
NOTE: The comments in this function were written in context with mass spectrometry so they are not very general
'''
def DataSmoother(data,abscissa,headers,dataSmootherChoice,dataSmootherTimeRadius,dataSmootherPointRadius,headersToConfineTo,polynomialOrder = 1):
    smoothedData = copy.deepcopy(data)
    if dataSmootherChoice == 'timerange':
        dataSmootherRadius = dataSmootherTimeRadius
    if dataSmootherChoice == 'pointrange':    
        dataSmootherRadius = dataSmootherPointRadius

    ## Option # 1
    #This if statement is for the first two possibilities- if the user does not only want a specific point
    #moved, but rather, all of the abscissa points changed, i.e. all mass series in data smoothed for XYY data
    if headersToConfineTo == []:

        # Get a list of time and data windows for each time in the abscissa
        # these windows are used to perform the smoothing

        (dataWindowsXvaluesInArrays, dataWindowsYYYYvaluesInArrays) = GetDataWindows(data,abscissa, dataSmootherRadius, dataSmootherChoice)


        # replace data points with their smoothed counterparts creating 'smoothedData'
        for timecounter, timeslist in enumerate(dataWindowsXvaluesInArrays):
            
            currentWindowAsArray = dataWindowsYYYYvaluesInArrays[timecounter]
            # Smooth the data
            smoothedData[timecounter,:] = DataSmootherPolynomialSmoothing(timeslist, currentWindowAsArray, abscissa[timecounter], polynomialOrder)   

    ## Option #2
    # if only specific masses in the data should be smoothed
    # the masses to be smoothed are listed in headersToConfineTo
    elif headersToConfineTo != []:
        # 'extractedData' is 'data' with the columns that will not be smoothed removed
        # 'extractedColumnTracker' is a list of 1 and 0 the same length as 'headers'
        # there is a 1 where the mass data is to be smoothed and a 0 when the mass data should not be smoothed
        (extractedData,extractedColumnTracker) = GetRelevantData(data, abscissa, headers, headersToConfineTo)

        # Now everything proceeds exactly as in Option #1, but use 'extractedData' rather
        # than 'data' and after we find 'smoothedExtractedDataData'
        # we will use it to replace the unsmoothed sections of 'data' where required
        # thus forming 'smoothedData'

        # Get a list of time and data windows for each time in the abscissa
        # these windows are used to perform the smoothing
        (dataWindowsXvaluesInArrays, dataWindowsYYYYvaluesInArrays) = GetDataWindows(extractedData,abscissa, dataSmootherRadius, dataSmootherChoice)


        # replace extractedData points with their smoothed counterparts creating 'smoothedData'
        smoothedExtractedData = numpy.zeros(extractedData.shape)
        for timecounter, timeslist in enumerate(dataWindowsXvaluesInArrays):
            
            currentWindowAsArray = dataWindowsYYYYvaluesInArrays[timecounter]
            # Smooth the data
            smoothedExtractedData[timecounter,:] = DataSmootherPolynomialSmoothing(timeslist, currentWindowAsArray, abscissa[timecounter], polynomialOrder)
    
        # Now construct smoothedData from tempSmootedData, data and extractedColumnTracker
        # i.e. combine the columns we just smoothed with the columns that didnt need smoothing
        # and do it in the proper/original order
        smoothedData = ReconstructSmoothedData(smoothedData, smoothedExtractedData, extractedColumnTracker)
           
    return smoothedData

''' 
The below funcitons are used for the marginalChangeRestrictor. The marginalChangeRestrictor prevents two points from differing by more than the 
factor marginalChangeRestriction by inserting an interpolated abscissa value and the corresponding interpolated YYYdata.
NOTE: marginalChangeRestrictor was originally intended for mass spec data so the comments below may correspond with the use with mass spec.
'''

#The analyticalLinearInterpolator will interpolate a single value or new row based on a value to be inserted between two points. It finds the slope
#of a theoretical line between any point in ordinateValues 1 and the corrsponding point in ordinateValues2 for all of the values across the row. It
#then uses the y=mx+b form to find the vector of y-intercepts for the slopes mentioned previously. Using the same y=mx+b form, it creates a vector 
#of the interpolated points across the row corresponding to the x value of the abscissaValueToInterpolateTo 
def analyticalLinearInterpolator (ordinateValues1, ordinateValues2, abscissaValueToInterpolateTo, abscissaValue1, abscissaValue2):
        slopesVector=(ordinateValues2-ordinateValues1)/(abscissaValue2-abscissaValue1)
        interceptVector=ordinateValues1-slopesVector*abscissaValue1
        interpolatedOrdinateValues=slopesVector*abscissaValueToInterpolateTo+interceptVector
        return interpolatedOrdinateValues
    
#This is a helper function made to find the sign of a certain element and return a multiplication factor containing the sign of the element.
#It returns -1 if the value is negative and 1 if the value is positive.
def signOfElement(element):
    sign=1.0
    if (element<0):
        sign=-1.0
    return sign

#For any zero or near-zero value, this funciton inserts rows above and below the row containing the zero/near-zero contianing half of the ignorableDeltaYThreshold directly above
#or below the value in question to prevent the interpolator from interpolating between the two points.
def adjustForNearZerosAndZeros(data, abscissa, IgnorableDeltaYThreshold):
    
    #Step 3.a
    #The nearZeroIndicesArray is created to recod the location of any zeros in the original data. A one represents a non-zero greater than the ignorableDeltaYThreshold while a 
    #zero represents a zero/nearZero This is used to find the locations where rows need to be added in this pre-processing step. It is not used outside of this function.
    nearZeroIndicesArray=(abs(data)>=IgnorableDeltaYThreshold)/1.0
    
    
    #The following statements convert the arrays into lists that would be more effecient for inserting rows
    nearZeroIndicesList=list(nearZeroIndicesArray)
    dataList=list(data)
    abscissaList=list(abscissa)
    
    #Since for loops do not revaluate the len(dataList) everytime it iterates, a number of rows added has to be included that ensures the commands contianed within
    #the loop only operate on the original rows in the dataList.
    rowsAdded=0
    
    #Step 3.b.1
    #This for loop iterates throught the rows of the original dataList and adds rows above and below if a zero or near-zero value is included in the data. The locations
    #of the zeros are recorded in the zeroIndicesArray, while the locations of zeros and near-zeros are kept in the nearZeroIndicesArray. The rows that are isnerted
    #have half of the threshold value inserted above and below the zero/nearZero with the other values in the row being linearly interpolated.
    for rowCounter in range(len(dataList)):
        
        #Keeps track of the current_row and abscissa for the row that is being evaluated
        current_row=dataList[rowCounter+rowsAdded]
        current_row_abscissa=abscissaList[rowCounter+rowsAdded]
        
        #Lists of the interpolated abscissa values that need to be inserted above or below a row based on the inserted half-threshold value.
        abscissaToInsertBelow=[]
        abscissaToInsertAbove=[]
        
        #Array of locations of zeros/nearZeros in the current row that need rows inserted above and/or below
        nearZeroIndicesRow=nearZeroIndicesList[rowCounter]
        
        #This section inserts the rows above if it is not the first row.
        if rowCounter !=0:
            
            previous_row=dataList[rowCounter+rowsAdded-1]
            previous_row_abscissa=abscissaList[rowCounter+rowsAdded-1]
            
            #Step 3.b.i
            #Loops through the values in the nearZeroIndicesRow to determine if any zeros/ are present in the row
            for columnCounter, nearZeroIndexValue in enumerate (nearZeroIndicesRow):
            
                absDifferenceAbove=abs(current_row[columnCounter]-previous_row[columnCounter])
            
                #If a zero/nearZero is present AND the difference between the current row and the above row is significant(greater than the IgnorableDeltaYThreshold)
                #a new abscissa will be calculated based on the halfMinThreshold to be inserted.
                if nearZeroIndexValue==0 and absDifferenceAbove>IgnorableDeltaYThreshold:
                
                    #Makes the sign of the halfMinThreshold the same as the one above.
                    halfMinThresholdAbove=(IgnorableDeltaYThreshold/2)*signOfElement(previous_row[columnCounter])
                    #Interpolates the new abscissa value
                    interpolatedAbscissaValue=analyticalLinearInterpolator(previous_row_abscissa,current_row_abscissa, halfMinThresholdAbove, previous_row[columnCounter], current_row[columnCounter])
                    abscissaToInsertAbove.append(interpolatedAbscissaValue)
            
            #Step 3.b.ii
            #Sorts the abscissa values so they are isnerted in the right order (increasing) and the rest of the row can be interpolated based off of the abscissa value.
            abscissaToInsertAbove=sorted(abscissaToInsertAbove)
        
            
            #Iterate through the interpolated abscissa values that need to be inserted
            for abscissaToInsertValue in abscissaToInsertAbove:
                
                #Step 3.b.iii
                #Inserts the new row in the dataList containing the interpolated values that correspond to the value in the abscissa.
                interpolatedDataRow=analyticalLinearInterpolator(previous_row, current_row , abscissaToInsertValue, previous_row_abscissa, current_row_abscissa )
                #Step 3.b.iiv
                dataList.insert(rowCounter+rowsAdded, interpolatedDataRow )
                
                #Insertes the values in the abscissa list in the order determined by the sorted abscissa to insert
                abscissaList.insert(rowCounter+rowsAdded, abscissaToInsertValue)
                
                rowsAdded+=1
        
        #This section basically repeats the same process as above, but for inserting rows below if the current_row is not the last row.        
        if rowCounter !=len(abscissa)-1:
            
            next_row=dataList[rowCounter+rowsAdded+1]
            next_row_abscissa=abscissaList[rowCounter+rowsAdded+1]
            
            #Step 3.b.i
            #Loops through the values in the nearZeroIndicesRow to determine if any zeros/ are present in the row
            for columnCounter, nearZeroIndexValue in enumerate (nearZeroIndicesRow):
            
                absDifferenceBottom=abs(current_row[columnCounter]-next_row[columnCounter])
                
                #If a zero/nearZero is present AND the difference between the current row and the above row is significant(greater than the IgnorableDeltaYThreshold)
                #a new abscissa will be calculated based on the halfMinThreshold to be inserted.
                if nearZeroIndexValue==0 and absDifferenceBottom>IgnorableDeltaYThreshold:
                    
                    #Makes the sign of the halfMinThreshold the same as the one above.
                    halfMinThresholdBelow=(IgnorableDeltaYThreshold/2)*signOfElement(next_row[columnCounter])
                    #Interpolates the new abscissa value
                    interpolatedAbscissaValue=analyticalLinearInterpolator(current_row_abscissa,next_row_abscissa, halfMinThresholdBelow, current_row[columnCounter], next_row[columnCounter])
                    abscissaToInsertBelow.append(interpolatedAbscissaValue)
            
            #Step 3.b.ii
            #Sorts the abscissa values so they are isnerted in the right order (increasing) and the rest of the row can be interpolated based off of the abscissa value.
            abscissaToInsertBelow=sorted(abscissaToInsertBelow)
        
            #Iterate through the interpolated abscissa values that need to be inserted
            for counter, abscissaToInsertValue in enumerate (abscissaToInsertBelow):
                
                #Step 3.b.iii
                #Inserts the new row in the dataList containing the interpolated values that correspond to the value in the abscissa.
                interpolatedDataRow=analyticalLinearInterpolator(next_row, current_row , abscissaToInsertValue, next_row_abscissa, current_row_abscissa )
                #Step 3.b.iiv
                dataList.insert(rowCounter+1+rowsAdded, interpolatedDataRow )
                
                #Insertes the values in the abscissa list in the order determined by the sorted abscissa to insert
                abscissaList.insert(rowCounter+1+rowsAdded, abscissaToInsertValue)
                
                rowsAdded+=1

    #Makes the lists back into arrays
    abscissa=numpy.array(abscissaList)
    data=numpy.array(dataList)
    
    #Step 3.c
    #Remakes the zeroIndicesArray and nearZeroIndicesArray so they include the rows inserted
    zeroIndicesArray=(data!=0.0)/1.0
    
    #Step 3.d
    #sets all the zeros in the data set equal to 1/10th of the threshold since the marginalChangeRestrictor cannot handle zeros in its row ratio.                
    data[data==0]=IgnorableDeltaYThreshold/10

    return data,abscissa,zeroIndicesArray

#Commented out because no longer in use.
##Returns true when a sign change occurs between two values
#def signChangePresent(value1,value2):
#    signChangePresent=False
#    if value1*value2<0:
#        signChangePresent=True
#        return signChangePresent

#AdujustForSignChange inserts a new row when a sign change is present and the threshold value is crossed. This row contains a zero at the location of the sign change.
#The rest of the values in the row are linearly interpolated.
def adjustForSignChange(data, abscissa, IgnorableDeltaYThreshold):
		
    dataList=list(data)
    abscissaList=list(abscissa)
    
    #Since for loops do not revaluate the len(dataList) everytime it iterates, a number of rows added has to be included that ensures the commands contianed within
    #the loop only operate on the original rows in the dataList.
    rowsAdded=0
    
    #Step 2.a
    #This for loop iterates throught the rows of the original dataList and adds a row below is there is a sign change occuring between any Y-Columns in the X-Row.
    #This row contains a zero in the location of the sign change and the rest of the values are linearly interpolated.
    for rowCounter in range(len(dataList)-1):
        
        #Keeps track of the current_row and abscissa for the row that is being evaluated
        current_row=dataList[rowCounter+rowsAdded]
        current_row_abscissa=abscissaList[rowCounter+rowsAdded]
        
        #Lists of the interpolated abscissa values that need to be inserted above or below a row based on the inserted 0.
        abscissaToInsert=[]
        
        
        next_row=dataList[rowCounter+rowsAdded+1]
        next_row_abscissa=abscissaList[rowCounter+rowsAdded+1]
        
        #Step 2.a.i
        #The below creates vectors of boolean variables. The signChangeVector will return true at a location of a sign chnage across the row. The 
        #absDifferenceSignificant vector will return true at the location of any differences across the row that are greater than the IgnorableDeltaYThreshold.
        signChangeVector=(current_row*next_row<0)
        absDifferenceVector=abs(current_row-next_row)
        absDifferenceSignificant=absDifferenceVector>IgnorableDeltaYThreshold
        
        if any(signChangeVector) and any(absDifferenceSignificant):
            #Loops through the current row to determine if a sign change is present and the difference larger than the threshold
            for columnCounter, current_rowValue in enumerate (current_row):
                
                #Step 2.a.ii
                #If a sign change is present AND the difference between the current row and the above row is significant(greater than the IgnorableDeltaYThreshold)
                #a new abscissa will be calculated based on the 0 to be inserted.
                
                 if signChangeVector[columnCounter] and absDifferenceSignificant[columnCounter]:   
                    #Interpolates the new abscissa value
                    interpolatedAbscissaValue=analyticalLinearInterpolator(current_row_abscissa,next_row_abscissa, 0, current_rowValue, next_row[columnCounter])
                    abscissaToInsert.append(interpolatedAbscissaValue)
            
            #Sorts the abscissa values so they are isnerted in the right order (increasing) and the rest of the row can be interpolated based off of the abscissa value.
            abscissaToInsert=sorted(abscissaToInsert)
            
            #Step 2.a.iii
            #Iterate through the interpolated abscissa values that need to be inserted
            for counter, abscissaToInsertValue in enumerate (abscissaToInsert):
                
                #Inserts the new row in the dataList containing the interpolated values that correspond to the value in the abscissa.
                interpolatedDataRow=analyticalLinearInterpolator(next_row, current_row , abscissaToInsertValue, next_row_abscissa, current_row_abscissa )
                dataList.insert(rowCounter+1+rowsAdded, interpolatedDataRow )
                
                #Insertes the values in the abscissa list in the order determined by the sorted abscissa to insert
                abscissaList.insert(rowCounter+1+rowsAdded, abscissaToInsertValue)
                
                rowsAdded+=1

    #Makes the lists back into arrays
    abscissa=numpy.array(abscissaList)
    data=numpy.array(dataList)
      
    #The analyticalLinearInterpolator will make values extremely close to zero, but not exactly zero, so this changes them to exactly zero.
    data[abs(data)<=(10**-12)]=0
               
    
    return data, abscissa

#The main interpolator inserts abscissa values for both the min and max log ratios. This creates the potential to insert superfluous interpolated rows.
#NOTE: rows refer to the YYY data at each abscissa value.
#removeSuperfluousInterpolatedRows removes the rows corresponding with the abscissa value if the current row and the row after the next row are within 
#the MaxAllowedDeltaYRatio and the IgnorableDeltaYThreshold. 
def superfluousInterpolatedRowCleanup(data, abscissa, insertion_indices, zeroIndicesArray, MaxAllowedDeltaYRatio, IgnorableDeltaYThreshold):
    
    #Make all the zeros in the YYY data equal to 1/10th of the ignorableDeltaYThreshold because a ratio cannot be taken with a zero in
    #any of the rows. This will be set back to zero at the end with the use of the zeroIndicesArray
    data[data==0]=IgnorableDeltaYThreshold/10
    
    #Make the data, abscissa, and zeroIndicesArray so rows can be effeciently deleted
    dataList=list(data)
    abscissaList=list(abscissa)
    zeroIndicesList=list(zeroIndicesArray)
    
    #Initialize rowsChecked variable that keeps track of how many rows were checked, but did not result in deletion
    rowsChecked=-1
    
    #Loops through the insertion indices so original data rows are not deleted
    for insertion_index in insertion_indices[:-2]:
        
        #Initializes the rows that need to be used to decide which rows to delete
        current_row=dataList[insertion_index+rowsChecked]
        next_next_row=dataList[insertion_index+rowsChecked+2]
        
        #Take the ratio of each element in the row in regards to the nex_next_row. This will determine if the two rows fall
        #within the required ratio without need for a row between them.
        row_ratio=next_next_row/current_row
        
        #Creates a boolean array. If an element ratio is between the MaxAllowedDeltaYRatio and the reciprocal of the MaxAllowedDeltaYRatio,
        #True is is put in the location of the element.
        ratioWithinLimits=[(ratio<=MaxAllowedDeltaYRatio and ratio>=MaxAllowedDeltaYRatio**-1.0) for
                           ratio in row_ratio]
        
        #Take the difference between the current_row and the next_next_row
        rowDifference=numpy.abs(current_row-next_next_row)
        
        #Creates a booolean array with True if the difference is significant(larger than the ignorableDeltaYThreshold)
        significantDifference=[(difference>IgnorableDeltaYThreshold) for difference in rowDifference]
        
  
        #If the current_row and the next_next_row are within the ratio and threshold limits and there is no peak present between the rows, delete the next_row from the 
        #data, abscissa, and zeroIndicies
        if all(ratioWithinLimits) and all(significantDifference):
            del dataList[insertion_index+rowsChecked+1]
            del abscissaList[insertion_index+rowsChecked+1]
            del zeroIndicesList[insertion_index+rowsChecked+1]
        else:
            #Since the row was checked, but no row deleted, add to the rowsChecked variable.
            rowsChecked+=1
    
    #Return the data. abscissa, and zeroIndicesArray to lists        
    data=numpy.array(dataList)
    abscissa=numpy.array(abscissaList)
    zeroIndicesArray=numpy.array(zeroIndicesList)
    
    #Remake the originally zero values in the data.
    data=data*zeroIndicesArray
    
    return data,abscissa

# This function is called once an unacceptable gap has been found.
# the limits of that gap are lower_point and upper_point, their
# corresponding abscissa points are the other parameters.
# The function performs the interpolation using Charles' algorithm for
# the data. Then it interpolates that to find the appropriate new abscissa points.
def form_new_abscissa(lower_point, upper_point,
                      abscissa_next_point, abscissa_current_point,
                      log_ratio, factor):

    # Now use Charles' algorithm to find the number of points needed
    # to fill the gap sufficiently.
    ## According to Charles if the largest ratio is exaclty an integer
    ## we should subtract one and it will work perfectly,
    ## else just round down to the nearest integer.
    if log_ratio.is_integer():
        number_needed = int(numpy.abs(log_ratio) - 1.0)
    else:
        number_needed = int(math.floor(numpy.abs(log_ratio)))
            
    # Construct a 1-d list that fills the gap
    # Include the data points from the lower and upper rows
    filled_column_gap = [lower_point]

    for new_point in range(number_needed):
        filled_column_gap.append(filled_column_gap[new_point]*factor)

    filled_column_gap.append(upper_point)
            
    # Now interpolate this data to find the corresponding points
    # for modified abscissa.
    ## NOTE: For some reason numpy interpolate can only function
    ## with increasing first argument, hence the if statements
    if lower_point < upper_point:
        filled_abscissa_gap = numpy.interp(filled_column_gap,
                                        [lower_point,upper_point],
                                        [abscissa_next_point,
                                         abscissa_current_point])
    elif lower_point >  upper_point:
        filled_abscissa_gap = numpy.interp(list(reversed(filled_column_gap)),
                                        [upper_point,lower_point],
                                        [abscissa_current_point,
                                         abscissa_next_point])
        filled_abscissa_gap = list(reversed(filled_abscissa_gap))

    # return only the points inside the gap not the bounding points
    # which we know from the original data
    return list(filled_abscissa_gap[1:-1])
        

##For marginalChangeRestrictor to work properly, the YYY data must be in numpy array form
#The marginalChangeRestrictor prevents two points from differing by more than the factor marginalChangeRestriction by inserting 
#an interpolated abscissa value and the corresponding interpolated YYYdata.
#NOTE: rows refer to the YYY data corresponding to each abscissa value.
def marginalChangeRestrictor(data,abscissa,MaxAllowedDeltaYRatio=2.0, IgnorableDeltaYThreshold=0.0001, cleanupSuperfluousInterpolatedRows=True, extraOutput=False):
    
    # First probably should make sure that the abscissa and the data
    # have the same number of rows (abscissa should index the rows of data)
    if (not data.shape[0] == abscissa.shape[0]):
        raise ValueError("data and abscissa dimensions do not match" +
                         " in the marginalChangeRestrictor function")
    
    # This set of statments runs the data through all of the necessary preprocessing
    # including casting all of the data to floats, removing sign changes with differences
    # greater than the IgnorableDeltaYThreshold, and removing any zeros from the data set.
    #Step 1: Cast data and abscissa to floats
    data=data.astype(float)
    abscissa=abscissa.astype(float)
    
    #Step 2: Remove Sign changes from Data
    signChangeData=adjustForSignChange(data,abscissa,IgnorableDeltaYThreshold)
    data,abscissa=signChangeData[0],signChangeData[1]

    
    #Step 3: Remove the Effect of zero/nearZero values through row insertion
    zeroAdjustment=adjustForNearZerosAndZeros(data,abscissa,IgnorableDeltaYThreshold)
    data,abscissa,zeroIndicesArray=zeroAdjustment[0],zeroAdjustment[1],zeroAdjustment[2]


    #Step 4: Go through the actual data interpolation

    #Step 4.a
    #Initialize list to be added to in the main interpolation
    # Array to eventually be returned with interpolated data
    interpolated_data = data[0,:]
    # list to eventually be the interpolated abscissa
    interpolated_abscissa = [abscissa[0]]
    # list to hold the filled gap chunks of data for
    # one final insert a the end
    data_gaps = []
    # insertion indices will store the row index
    # where the filled data should be inserted to the main data
    # array
    insertion_indices = []
    
    #Step 4.b.
    dataList=list(data)
    abscissaList=list(abscissa)

    #Step 4.c.
    # loop through abscissa, at each step check the data
    # if any gaps are too large determine the data that fills
    # the worst gap sufficiently
    # then interpolate those new data points to the 'abscissa_values_for_gap_combined'.
    # Finally use abscissa_values_for_gap_combined to interpolate to all of the data columns
    # Note: skip iteration of last point so we dont overrun abscissa
    for point_index,point in enumerate(abscissa[0:-1]):

        # get the relevant rows
        current_row = dataList[point_index]
        next_row = dataList[point_index + 1]
        current_zeros_row=zeroIndicesArray[point_index]
        next_zeros_row=zeroIndicesArray[point_index+1]
        
        # abscissa points
        abscissa_current = abscissaList[point_index]
        abscissa_next = abscissaList[point_index+1]

        #Step 4.c.i.
        # get the absolute difference
        # (by absolute I mean both absolute value and
        # arithmetic difference in contrast to the ratios below)
        row_absolute_difference = numpy.abs(next_row - current_row)

        # Determine which gaps exceed the absolute threshold
        # (likely most will/should)

        absolute_gaps = [(significantDifference > IgnorableDeltaYThreshold) for
                             significantDifference in row_absolute_difference]
        
        
        #Step 4.c.i.1
        #Create a row with the locations of any of differences between the 
        #next_row and the current_row that are less than the IgnorableDeltaYThreshold
        negligibleAbsoluteDifferencesArray=[]
        for significantDifference in absolute_gaps:
            if not significantDifference:
                negligibleAbsoluteDifferencesArray.append(0)
            else:
                negligibleAbsoluteDifferencesArray.append(1)
        
        #Step 4.c.ii
        #Take the ratio of the values in the next_row to the corresponding values
        #in the current_row
        row_ratio = next_row/current_row
        
        #This makes the row_ratio of sign changes within the threshold equal to
        #1 so the remainder of the function does not try to interpolate between
        #points having a change in sign.
        row_ratio[row_ratio<0]=1
        
        #Step 4.c.ii.1
        # take the log of the ratio in the desired base
        # Note: I dont think that there is a vectorized numpy log_n function
        ratio_logarithm = numpy.array(
            [math.log(data_point, MaxAllowedDeltaYRatio)
             for data_point in row_ratio]
            )
        
        #Step 4.c.iii
        #Multiplies the ratio_logarrithm by the current_zeros_row and the 
        #nearZeroDifference to prevent the interpolator from interpolating 
        #between points that are in the ignorableDeltaYThreshold
        ratio_logarithm=ratio_logarithm*current_zeros_row*negligibleAbsoluteDifferencesArray

        
        # is there any abs(log(ratio)) that exceed 1.0?
        # that would indicate that the gap in that column is too large
        # for the relative threshold
        relative_gaps = [(numpy.abs(ratio) > 1.0) for ratio in ratio_logarithm]

        # combine the relative and absolute gap
        gaps = [(sig_rel_gap_present and sig_abs_gap_present) for (sig_rel_gap_present,sig_abs_gap_present) in
                    zip(relative_gaps, absolute_gaps)]
        
        #Step 4.c.iii.1
        #returns the zero values in the current and next rows back to zero. Also alters the value in the current row in the actual dataList.
        current_row=current_row*current_zeros_row
        dataList[point_index]=dataList[point_index]*current_zeros_row
        next_row=next_row*next_zeros_row
        
        #Step 4.c.i.2 
        #If there are problem gaps that exceed the absolute threshold
        # and the relative threshold deal with them (i.e. add the points
        # to the modified_abscissa that we need)
        if any(gaps):
            # So now there is a gap (i.e. large gap in at least in one column).
            # If there is a problem gap in both directions (e.g. gap1 == [2,12]
            # and gap2 == [7,1.5]) then we will need to perform the interpolation
            # and backwards interpolation to form the modified_abscissa for both
            # gaps, then just combine the two in an ordered fashion and append
            # like normal
            
            #Step 4.c.iv
            # keep only the column indices (in this pair of rows) which had a gap
            gap_indices = [index for (index,gap)
                                     in enumerate(gaps) if gap]

            gap_log_ratios = [ratio_logarithm[gap] for
                                 gap in gap_indices]

            # Find the maximum and miniumum of the gap log(ratio)
            # dont care if they are unique, if there are two identical extrema
            # they will require identical interpolation anyway
            max_log_ratio = max(gap_log_ratios)
            min_log_ratio = min(gap_log_ratios)

            # Doesn't need to be unique
            max_log_ratio_index = gap_indices[gap_log_ratios.index(max_log_ratio)]
            min_log_ratio_index = gap_indices[gap_log_ratios.index(min_log_ratio)]

            abscissa_values_for_gap_from_rising_data = []
            abscissa_values_for_gap_from_falling_data = []
            
            #Step 4.c.v
            # if there is an increasing gap find the points to add to
            # modified abscissa
            if max_log_ratio > 0:
                abscissa_values_for_gap_from_rising_data = form_new_abscissa(current_row[max_log_ratio_index],
                                                   next_row[max_log_ratio_index],
                                                   abscissa[point_index],
                                                   abscissa[point_index + 1],
                                                   max_log_ratio,
                                                   MaxAllowedDeltaYRatio)
            # if there is a decreasing gap find the points to add to
            # modified abscissa
            if min_log_ratio < 0:
                abscissa_values_for_gap_from_falling_data = form_new_abscissa(current_row[min_log_ratio_index],
                                                   next_row[min_log_ratio_index],
                                                   abscissa[point_index],
                                                   abscissa[point_index + 1],
                                                   min_log_ratio,
                                                   MaxAllowedDeltaYRatio**-1.0)

            #Step 4.c.v.1
            # Combine the two if there is both and remove any duplicates
            abscissa_values_for_gap_combined = abscissa_values_for_gap_from_rising_data + abscissa_values_for_gap_from_falling_data
            abscissa_values_for_gap_combined = list(set(abscissa_values_for_gap_combined))
            
           
            # the abscissa needs to be ordered and increasing or decreasing
            # depending on how whether the boundary abscissa points are
            # increasing or decreasing.
            if abscissa[point_index] < abscissa[point_index + 1]:
                abscissa_values_for_gap_combined = sorted(abscissa_values_for_gap_combined)

            if abscissa[point_index] > abscissa[point_index + 1]:
                abscissa_values_for_gap_combined = list(reversed(sorted(abscissa_values_for_gap_combined)))

            
            #Step 4.c.v.2       
            
            #The following code vectorizes the linear row interpolation based on the abscissa values to be inserted.
            #The slopesVector hold the slopes for each column of the YYY data
            slopesVector=(next_row-current_row)/(abscissa_next-abscissa_current)
            #interceptVector hold the linear intercepts for each column of the YYY data
            interceptVector=current_row-slopesVector*abscissa_current
            
            #The duplicatedAbscissaArray is created to easily perform calculations for the entire list of abscissa values to be inserted at once.
            #There is a numpy array row for each YYY column and a column for each abscissa value to be inserted. The each column in the 
            #duplicatedAbscissaArray will have the same value in each of the rows. This value is the abscissa value to be inserted and the the 
            #columns of the YYY data to be interpolated to
            #Thes
            duplicatedAbscissaArray=numpy.ones((len(slopesVector),len(abscissa_values_for_gap_combined)))*abscissa_values_for_gap_combined
            
            
            #The interpolated rows of the YYY data are calculated with a y=mx+b equation using the slopes and intercepts found above. To calculate
            #correctly, the duplicatedAbscissaArray must be transposed.
            filled_data=numpy.multiply(numpy.transpose(duplicatedAbscissaArray),slopesVector)+interceptVector
            
            #Step 4.c.v.3
            # Add the filled data and abscissa points to the
            # interpolated data and abscissa (we handle the
            # boundary points elsewhere)
            #data_gaps.append(filled_data.reshape(len(abscissa_values_for_gap_combined),len(current_row)))
            for row in filled_data[:,]:
                data_gaps.append(row)
                insertion_indices.append(point_index+1)
            #interpolated_data = numpy.vstack((interpolated_data, filled_data))
            interpolated_abscissa.extend(list(abscissa_values_for_gap_combined))
            
        # END if gap block

        # Regardless of the existence of a gap in this row
        # make sure to add the upper row to the interpolated_data
        # and the upper point to abscissa
        interpolated_data = numpy.vstack((interpolated_data, next_row))
        interpolated_abscissa.extend([abscissa_next])
        
        
    #Else if there are x-rows that need to be inserted, insert the rows and remake the zeroIndicesArray before
    #returning the data
    if len(insertion_indices)!=0:
        #Step 4.d.
        #For the proper data_gaps to be inserted into the proper insertion points in the list,
        #the for loop has to be iterated backwards. This must have to do with the order that the
        #insertion_indices and data_gaps are filled.
        for counter in range(len(insertion_indices)-1,-1,-1):
            dataList.insert(insertion_indices[counter],data_gaps[counter])
            
        #The for loop does not return any originally zero values in the last row of the dataList back to zero.
        #This statement multiplies the last row in the dataList by the last row in the zeroIndiciesArray to return
        #any originally zero values in this row back to zero.
        dataList[-1]=dataList[-1]*zeroIndicesArray[-1]
        
        #return the lists to arrays
        data=numpy.array(dataList)
        abscissa=numpy.array(interpolated_abscissa)
        
        #Step 4.e
        #create a fresh zeroIndicesArray contianing the the zero locations with the inserted rows
        zeroIndicesArray=(data!=0.0)/1.0

        #Step 5
        
        #If the user wants to remove the superfluous rows created from the main interpolator, the cleanupSuperfluousInterpolatedRows will be defined as
        #True. The default value for this is set to True.
        #NOTE: rows refer to the YYY data corresponding to each abscissa value.
        if cleanupSuperfluousInterpolatedRows:
            
            dataWithsuperfluousRowsRemoved=superfluousInterpolatedRowCleanup(data, abscissa, insertion_indices, zeroIndicesArray, MaxAllowedDeltaYRatio, IgnorableDeltaYThreshold)
            data=dataWithsuperfluousRowsRemoved[0]
            abscissa=dataWithsuperfluousRowsRemoved[1]
    

    #Else no rows needed to be interpolated above, and there would be no rows that need to be inserted so the 
    #data does not need to be altered and can be returned.
    else:
        data = numpy.multiply(data,zeroIndicesArray)

    #For unit tester purposes, unit tester will also compare the insertion indices
    if extraOutput:
        return data, abscissa, insertion_indices    
    #Step 6
    return data, abscissa

#interpolateAccompanyingArrays is a preprocessing step that changes the abscissa of the accompanying concentration bounds file to 
#match that of the interpolated data. When the times do not match that of the marginalChangeRestricted_abscissa, it will add a new row full of 
#linearInterpolations across the YYY columns of the bounds and other data contained in the accompanyingArray file.

#The abscissa in the accompanyingArray must match that of the actual data pre-marginalChangeRestrictor.
#The first and last number in the csv file must match that of the marginalChangeRestricted_abscissa. 
#The original intention of this function was to expand the concentrationBoundsFromCSV to have the same abscissa as the data following the
#marginalChangeRestrictor, but the function could be potentially used for additional accompanying arrays that need a matching abscissa.

def interpolateAccompanyingArrays(marginalChangeRestricted_abscissa, accompanyingArray):
    
    #The dataRangeSpecifierInterpolator will progress if points were added in the interpolator. The easiest way to do this is to
    #compare the lengths of the interpolated abscissa to that of the datafromcsv
    if len(marginalChangeRestricted_abscissa)!=len(accompanyingArray):
        accompanyingList=list(accompanyingArray)

        #Loops through the marginalChangeRestricted_abscissa.
        for rowCounter, marginalChangeRestricted_abscissa_value in enumerate(marginalChangeRestricted_abscissa[:-1]):
            
            #Check if the marginalChangeRestricted_abscissa_value matches the abscissa of the same row in the datafromcsvfile
            if marginalChangeRestricted_abscissa_value !=accompanyingList[rowCounter][0]:
                
                #Set the previous and current row to be used to interpolate between them based on the marginalChangeRestricted_abscissa_value that did not
                #match the value in the file
                previous_row=accompanyingList[rowCounter-1]
                current_row=accompanyingList[rowCounter]
                
                #Use the analyticalLinearInterpolator to interpolate 
                interpolatedExpandedRow=analyticalLinearInterpolator(previous_row,current_row,marginalChangeRestricted_abscissa_value,previous_row[0],current_row[0])
                accompanyingList.insert(rowCounter,interpolatedExpandedRow)

        
        expandedArray=numpy.array(accompanyingList)        
                  
    return expandedArray

'''RemoveSignals is used by itrative analysis post processing to subtract the signals for the solved molecules from the experimental data.'''
def RemoveSignals(dataToExtractFrom, totalColumns, dataToExtract, columnsToSubtract):
    #all inputs must be numpy arrays i.e. Func(2D(MxN),1D(N),2D(MxV),1D(V))
    
    #confirm that data lengths are the same
    if (len(dataToExtractFrom[:,0]) != len(dataToExtract[:,0])):
        print("Length of simulated data doesn't equal length of experimental data")
    #copy the row abscissa to be readded later
    rowAbscissa = dataToExtract[:,[0]]
    #Remove row abscissa from each data set
    #dataToExtractFrom = dataToExtractFrom[:,1:]
    dataToExtract = dataToExtract[:,1:]
    
    #confirm that total data sets are the same width
    if (len(dataToExtractFrom[0,:]) != len(totalColumns)):
        print("Number of signals in total experimental data doesn't match header")
    
     #confirm that solved data sets are the same width
    if (len(dataToExtract[0,:]) != len(columnsToSubtract)):
        print("Number of signals in solved experimental data doesn't match header")
        
    #for each solved column
    for solvedColumnNumber in range(len(columnsToSubtract)):
        #Find the matchng total data index
        index = numpy.where(columnsToSubtract[solvedColumnNumber] == totalColumns)[0]
        index= int(index)
        #subtract that solved column from the total data
        dataToExtractFrom[:,index] = dataToExtractFrom[:,index] - dataToExtract[:,solvedColumnNumber]
    #add the row abscissa back into the total data
    subtractedData = numpy.hstack((rowAbscissa,dataToExtractFrom))
    #return the values
    return subtractedData

'''Remove columns at zero or below threshold will delete python array columns (appearing as rows on spreadsheets) where the absolute value of all the values are
either are below a threshold.The absolute value was included so significant negative values will not be deleted. The default of the function will have a default 
threshold of zero. Additionally, a startingRowIndex (appears as columns on spreadsheets) can be used to evaluate a subset of the original array. This becomes 
useful for evaluating a XYYY data set, where the X-rows should not be evaluated.'''
def removeColumnsWithAllvaluesBelowZeroOrThreshold(dataValuesArray, startingRowIndex=0, threshold=0):
    #Turns the array into a list to make deletion easier
    dataValuesList=list(dataValuesArray)
    
    #Since the loop does not re-evaluate the length of the referenceDataList, an offset value must be used to account of the rows that were
    #removed
    columnsDeleted=0

        #Loops through the length of the data list using a rowCounter index
    for columnsCounter in range(len(dataValuesList)):
        #If all the intensity values for a certain mass fragment are below the threshold, the row gets deleted. Since the row is deleted, the
        #row index will have to remain the same to evaluate the next original row, meaning rowsDeleted has to be subtracted from the rowCounter.
        if all(abs(dataValuesList[columnsCounter-columnsDeleted][startingRowIndex:])<=threshold):
            del dataValuesList[columnsCounter-columnsDeleted]
            columnsDeleted+=1
    
    #Makes the list back into an array
    reducedDataValuesArray=numpy.array(dataValuesList)   
     
    return reducedDataValuesArray

'''Takes an XYYY array and another XYYY array (can be just XY), then creates a merged XYYY array.
Note that the abscissa do not have to match. A merged abscissa is made, an any missing values are filled with a third optional argument.
We assume that the abscissa are floats and then use the == comparison to check if they are the same.
With small values  and large values (e.g. 1E-8 and 1E9) there can be problems.
'''
def addXYYYtoXYYY(XYYYarray1, XYYYarray2, defaultValue=0):
    abscissa1 = numpy.array(XYYYarray1[:,0], dtype='float')
    abscissa2 = numpy.array(XYYYarray2[:,0], dtype='float')
    numOrdinatesArray1 = int(len(XYYYarray1[0]) - 1) #subtract 1 for abscissa. #number of ordinated values in array.
    numOrdinatesArray2 = int(len(XYYYarray2[0]) - 1) #subtract 1 for abscissa  #number of ordinated values in array.
    #Make the combined abscissa.
    abscissaCombined = set(abscissa1)|set(abscissa2) #This combines the two sets.
    abscissaCombined = list(abscissaCombined)
    abscissaCombined.sort() #This function returns a None type but sorts the list.
    abscissaCombined = numpy.array(abscissaCombined, dtype='float') #convert to list, then to numpy array.
    abscissaCombinedLength = len(abscissaCombined)
    #for each of the existing XYYY datasets, we will make an extended version (only extended in rows) where the blanks are filled.
    XYYYarray1extended = numpy.full((abscissaCombinedLength ,numOrdinatesArray1+1), defaultValue, dtype="float")
    XYYYarray1extended[:,0] = abscissaCombined*1.0
    for newAbscissaIndex,newAbscissaValue in enumerate(abscissaCombined):
        #Now check if newAbscissaValue occurs in abscissa1. Originaly tried using numpy.where, but it did not work. I did not understand why, but switched to using a loop.
        for oldAbscissaIndex, oldAbscissaValue in enumerate(abscissa1):
            if newAbscissaValue == oldAbscissaValue: #This means the value exists in abscissa1
                XYYYarray1extended[newAbscissaIndex][1:] =  XYYYarray1[oldAbscissaIndex][1:]*1.0#skip first column because it has the abscissa values in it.
            #else: #This else is implied because the array was initialized with the default value.
            #    XYYYarray1extended[newAbscissaIndex,1:] =  defaultValue
    XYYYarray2extended = numpy.full((abscissaCombinedLength ,numOrdinatesArray2+1), defaultValue, dtype="float")
    XYYYarray2extended[:,0] = abscissaCombined*1.0
    for newAbscissaIndex,newAbscissaValue in enumerate(abscissaCombined):
        #Now check if newAbscissaValue occurs in abscissa2. Originaly tried using numpy.where, but it did not work. I did not understand why, but switched to using a loop.
        for oldAbscissaIndex, oldAbscissaValue in enumerate(abscissa2):
            if newAbscissaValue == oldAbscissaValue: #This means the value exists in abscissa1
                XYYYarray2extended[newAbscissaIndex][1:] =  XYYYarray2[oldAbscissaIndex][1:]*1.0#skip first column because it has the abscissa values in it. 
            #else: #This else is implied because the array was initialized with the default value.
            #    XYYYarray2extended[newAbscissaIndex,1:] =  defaultValue
    #now just need to combine the extended arrays with a stacking command.
    combinedArray = numpy.zeros((abscissaCombinedLength,numOrdinatesArray1+numOrdinatesArray2+1), dtype="float")
    combinedArray[:,0] = abscissaCombined
    combinedArray[:,1:numOrdinatesArray1+1] = XYYYarray1extended[:,1:]
    combinedArray[:,numOrdinatesArray1+1:numOrdinatesArray1+numOrdinatesArray2+1] = XYYYarray2extended[:,1:]
    return combinedArray
    
    
'''Takes an XYYY array and a second abscissa (Z), then adds any rows to the first array that were not there.
#TODO: If abscissa2 is not provided, then all integer valuees between the max and min of the first X axis will be filled in.
We assume that the abscissa are floats and then use the == comparison to check if they are the same.
With small values  and large values (e.g. 1E-8 and 1E9) there can be problems.'''
def extendXYYYtoZYYY(XYYYarray1, abscissa2=None, defaultValue=0):
    abscissa1 = numpy.array(XYYYarray1[:,0], dtype='float')
    if type(abscissa2) != type(None):
        abscissa2 = numpy.array(abscissa2, dtype='float')
    #TODO: add an else statement for when abscissa2 is not provided.
    numOrdinatesArray1 = int(len(XYYYarray1[0]) - 1) #subtract 1 for abscissa. #number of ordinated values in array.
    
    #Get a list of new absicssa values to insert.
    valuesToExtendOriginalAbscissaBy = [] #just initializing.
    for abscissa2value in abscissa2:
        if abscissa2value not in abscissa1:
            valuesToExtendOriginalAbscissaBy.append(abscissa2value)

    #Now find the abscissa combined.
    abscissaCombined = set(abscissa1)|set(valuesToExtendOriginalAbscissaBy) #This combines the two sets.
    abscissaCombined = list(abscissaCombined)
    abscissaCombined.sort() #This function returns a None type but sorts the list.
    abscissaCombined = numpy.array(abscissaCombined, dtype='float') #convert to list, then to numpy array.
    abscissaCombinedLength = len(abscissaCombined)
    
    #Now need to extend the new array.
    #for each of the existing XYYY datasets, we will make an extended version (only extended in rows) where the blanks are filled.
    XYYYarray1extended = numpy.full((abscissaCombinedLength ,numOrdinatesArray1+1), defaultValue, dtype="float")
    XYYYarray1extended[:,0] = abscissaCombined*1.0
    for newAbscissaIndex,newAbscissaValue in enumerate(abscissaCombined):
        #Now check if newAbscissaValue occurs in abscissa1. Originaly tried using numpy.where, but it did not work. I did not understand why, but switched to using a loop.
        for oldAbscissaIndex, oldAbscissaValue in enumerate(abscissa1):
            if newAbscissaValue == oldAbscissaValue: #This means the value exists in abscissa1
                XYYYarray1extended[newAbscissaIndex][1:] =  XYYYarray1[oldAbscissaIndex][1:]*1.0#skip first column because it has the abscissa values in it.
            #else: #This else is implied because the array was initialized with the default value.
            #    XYYYarray1extended[newAbscissaIndex,1:] =  defaultValue
    return XYYYarray1extended

    
    