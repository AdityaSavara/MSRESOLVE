"""
Functions for use with MS data. All of these functions assume that the list of Y data lists 
have been turned into numpy arrays and transposed. Then the data series for each mass is a column 
in the new numpy array and the number of columns is the same as the number of elements of 
'choosenMassFragments' in UserInput.py.
"""

import numpy
import copy



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
Then for each inner time (i.e. element of 'timeslist') it 
takes the relevant data and  performs the smoothing and then returns the smoothedData.
'''
def DataSmootherPolynomialSmoothing(timeslist, currentWindowAsArray,  currentTime, polynomialOrder=1):

    # vectorized polyfit/polyval for all Y columns at this timecounter
    smoothedPoints = numpy.polyval(numpy.polyfit(timeslist,currentWindowAsArray,polynomialOrder), currentTime)
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
    for timecounter,_ in enumerate(abscissa):
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
of elements as there are Y in XYY data). headersToConfineTo
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
            extractedColumnTracker[headers.index(relevantHeader)] = 1.0 # this assumes no
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
KeepOnlySelectedYYYYColumns() compares the values of DataAbscissa and AbscissaValuesToKeep.
Values(numbers) that are in DataAbscissa and not in AbscissaValuesToKeep are deleted
from DataAbscissa. Further the columns of YYYYData correspond to the DataAbscissa array,
i.e. shape(YYYData)[1] = len(DataAbscissa). When a value is removed from 
DataAbscissa the corresponding column is removed from YYYYData. 

Parameters:
YYYYData- A 2-d numpy array, shape(YYYYData) = (*, len(DataAbscissa)), i.e. the columns
          of YYYYData correspond to the entries in DataAbscissa
DataAbscissa- List of floats. (or 1-d numpy array)
AbscissaValuesToKeep- List of floats. NOTE: AbscissaValuesToKeep is not necessarily 
                      smaller than DataAbscissa, it may contain values not included
                      in DataAbscissa.
'''
def KeepOnlySelectedYYYYColumns(YYYYData, DataAbscissa, AbscissaValuesToKeep):
    
    # list to track indices that should be deleted
    deletion_indices = []

    # loop through DataAbscissa, record the indices
    # of values not also found in AbscissaValuesToKeep
    # for deletion
    for (valueIndex,value) in enumerate(DataAbscissa):
        if value not in AbscissaValuesToKeep:
            deletion_indices.append(valueIndex)

    # Now remove the unwanted values from DataAbsicssa
    # and the corresponding unwanted columns from YYYYData
    DataAbscissa = numpy.delete(DataAbscissa,deletion_indices)
    YYYYData = numpy.delete(YYYYData,deletion_indices, axis=1)

    return (YYYYData, DataAbscissa)

#TODO: make a function KeepOnlyYYYYRows() that is very similar to this one.
# It may be useful and could replace some of the functionality of ArrayBuilder()
# and UnecessaryMoleculesDeleter()
    

'''
The DataSmoothing function 'smoothes' data over a certain time or datapoint ranges:
it goes through each mass fragment at a certain time and determines a polynomial that modeles the datapoints around
that mass fragment (within point or time radius). Then applies this determined polynomial to recalculate the value of the datapoint.
after all datapoint of a certain time are analyze the function then resets on the next datapoint. 
NOTE: The comments in this function were written in context with mass spectrometry so they are not very general
'''
def DataSmoother(data,abscissa,headers,dataSmootherChoice,dataSmootherTimeRadius,dataSmootherPointRadius,headersToConfineTo,polynomialOrder = 1):
    smoothedData = copy.deepcopy(data)

    ## Option # 1
    #This if statement is for the first two possibilities- if the user does not only want a specific point
    #moved, but rather, all of the abscissa points changed, i.e. all mass series in data smoothed for XYY data
    if headersToConfineTo == []:

        # Get a list of time and data windows for each time in the abscissa
        # these windows are used to perform the smoothing
        (dataWindowsXvaluesInArrays, dataWindowsYYYYvaluesInArrays) = GetDataWindows(data,abscissa, dataSmootherTimeRadius, dataSmootherChoice)


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
        (dataWindowsXvaluesInArrays, dataWindowsYYYYvaluesInArrays) = GetDataWindows(extractedData,abscissa, dataSmootherTimeRadius, dataSmootherChoice)


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
