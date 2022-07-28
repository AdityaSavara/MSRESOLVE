
#Most of the MSRESOLVE dependencies for ExtentOfSLSUniqueSolvable are related to data the MSReference class and data manipulation.
#that class and some functions are included in this file so that the BestMassFragChooser can be used without MSRESOLVE.
#This file was made on Nov 18, 2018, so it may not be the most up to date version of the MSReference class etc.

import pandas
import numpy
import copy
import XYYYDataFunctionsSG as DataFunctions
import bisect


'''
trimDataMassesToMatchChosenMassFragments() and trimDataMassesToMatchReference() are just a wrapper functions for calls to DataFunctions.KeepOnlySelectedYYYYColumns(). 
Both of the functions trim ExperimentData.workingData and ExperimentData.mass_fragment_numbers. 
The first function trims the data according to the mass fragment selections in G.chosenMassFragments.
The second function trims the data to remove any mass fragments for which there is no ReferenceData. 

Parameters:
ExperimentData - of type MSData, the one instantiated in main() named ExperimentData is a good example of one
    that will work here
ReferenceData - of type MSReference, ReferenceData from main() is a good example
chosenMassFragments  - list of integers, like the one created in UserInput  
'''
def trimDataMassesToMatchChosenMassFragments(ExperimentData, chosenMassFragments):
    # If we are only interested in a subset of the MS data
    # and that subset is a subset of the loaded data
    # remove the irrelevant mass data series from ExperimentData.mass_fragment_numbers
    # and the corresponding colums from ExperimentData.workingData
    trimmedExperimentData = copy.deepcopy(ExperimentData)
    #print("MassFragChooser")
    (trimmedExperimentData.workingData, trimmedExperimentData.mass_fragment_numbers) = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedExperimentData.workingData,
                                                                                                            trimmedExperimentData.mass_fragment_numbers,
                                                                                                            chosenMassFragments, header_dtype_casting=float)
    #!!! this line removed from stripped down version of MSRESOLVE functions/class to remove external depenendencies. trimmedExperimentData.ExportCollector("MassFragChooser")
    
    return trimmedExperimentData

def trimDataMassesToMatchReference(ExperimentData, ReferenceData):
    
    trimmedExperimentData = copy.deepcopy(ExperimentData)
    
    # Remove elements of ExperimentData.mass_fragment_numbers for which there is no matching mass in the reference data.
    # Also remove the corresponding mass data column from Experiment.workingData.
    (trimmedExperimentData.workingData, trimmedExperimentData.mass_fragment_numbers) = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedExperimentData.workingData,
                                                                                                        trimmedExperimentData.mass_fragment_numbers,
                                                                                                        ReferenceData.provided_mass_fragments, header_dtype_casting=float)

    return trimmedExperimentData

class MSReference (object):
    def __init__(self, provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies=None, moleculeIonizationType=None, mass_fragment_numbers_monitored=None, referenceFileName=None, form=None, AllMID_ObjectsDict={}):
        self.provided_reference_patterns, self.electronnumbers, self.molecules, self.molecularWeights, self.SourceOfFragmentationPatterns, self.sourceOfIonizationData, self.relativeIonizationEfficiencies, self.moleculeIonizationType, self.mass_fragment_numbers_monitored, self.referenceFileName, self.form = provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form
        #class object variable created to allow class to be used separately from the program. 
        self.ExportAtEachStep = ''
        self.iterationSuffix = ''
        #This loops through the molecules, and removes whitespaces from before and after the molecule's names.
        for moleculeIndex, moleculeName in enumerate(self.molecules):
            self.molecules[moleculeIndex] = moleculeName.strip()     
            
        '''Initializing Export Collector Variables'''
        #start the timer function
        #!!! this line removed from stripped down version of MSRESOLVE to remove external depenendencies. self.previousTime = timeit.default_timer()
        #initalize debugging lists
        #These lists are appended in parallel so a variable of the same index from each list will be related
        self.runTimeAtExport = []
        self.labelToExport = []
        self.dataToExport = []
        self.moleculesToExport = []
        self.exportSuffix = ''
        #self.experimentTimes = []       
        self.provided_mass_fragments = self.provided_reference_patterns[:,0]
        #Get ionization efficiencies and export their values and what method was used to obtain them
        #!!!!!This line is commented out in this stripped down version, in the real MSRESOLVE this line is important. self.populateIonizationEfficiencies(AllMID_ObjectsDict)
        #!!!!self.exportIonizationInfo()
    #TODO exportCollector should be updated to take in a string argument for the data type that it should record (patterns vs various intensities)
    #Additionally, it should take an optional variable to determine the headers that will be used.         
    def ExportCollector(self, callingFunction, use_provided_reference_patterns = False):
        #record current time
        currentTime = timeit.default_timer()
        #add net time to list of run times
        self.runTimeAtExport.append(currentTime - self.previousTime)
        #record current time for next function's use
        self.previousTime = currentTime
        #add the name of the calling function to mark its use
        self.labelToExport.append(callingFunction) 
        
        if self.ExportAtEachStep == 'yes':
            ##record molecules of experiment
            self.moleculesToExport.append(self.molecules.copy())
            #record data of experiment
            if use_provided_reference_patterns:
                self.dataToExport.append(self.provided_reference_patterns.copy())
            elif callingFunction == 'UnnecessaryMoleculesDeleter':
                self.dataToExport.append(self.monitored_reference_intensities.copy())
            elif not use_provided_reference_patterns:
                self.dataToExport.append(self.standardized_reference_patterns.copy())
            
    def ExportFragmentationPatterns(self, verbose=True):
        #Only print if not called from interpolating reference objects
        if verbose:
            print("\n Reference Debugging List:")
        for savePoint in range(len(self.runTimeAtExport)):
            #Only print if not called from interpolating reference objects
            if verbose:
                print(self.labelToExport[savePoint])
                print(self.runTimeAtExport[savePoint])
            if self.ExportAtEachStep == 'yes':
                #inserting the data for a particular savePoint
                filename = 'Exported%s%s.csv'%(savePoint, self.labelToExport[savePoint])
                data = self.dataToExport[savePoint]
                colIndex = ['%s'% y for y in self.moleculesToExport[savePoint]]
                #colIndex = ['%s'% y for y in self.molecules]
                #print(numpy.shape(data),numpy.shape(colIndex))
                ExportXYYYData(filename,data,colIndex, fileSuffix = self.iterationSuffix)

    # This class function removes all rows of zeros from
    # the XYYY sorted provided_reference_patterns, and *also* provided_mass_fragments
    #The logic in the below funtion is badly written, in terms of efficiency. But it seems to work at present.
    #TODO: This is not a good practice, because provided_reference_patterns is getting changed, no longer "Provided".
    #TODO: (continued from previous line) It's more like "zero_trimmed" reference intensities after this.
    def ClearZeroRowsFromProvidedReferenceIntensities(self):
        #initial a counter for the row index, which will be updated during the loop
        currentRowIndexAccountingForDeletions = 0
        #provided_reference_patternsOnly is not used, but is made for future use (see below)
        provided_reference_patternsOnly = self.provided_reference_patterns[:,1:]
        for intensitiesOnlyInRow in provided_reference_patternsOnly:
            #This line checks if there are any non-zeros in the row.
            numberOfNonzeros = numpy.count_nonzero(intensitiesOnlyInRow)
            if numberOfNonzeros == 0 :
                #If there are only zeros. we delete a row and adjust the row index to account for that deletion.
                self.provided_reference_patterns = numpy.delete(self.provided_reference_patterns, currentRowIndexAccountingForDeletions, axis=0 ) #axis = 0 specifies to delete rows (i.e. entire abscissa values at the integer of currentRowIndexAccountingForDeletions).
                self.provided_mass_fragments = numpy.delete(self.provided_mass_fragments, currentRowIndexAccountingForDeletions, axis=0 )
                currentRowIndexAccountingForDeletions = currentRowIndexAccountingForDeletions -1
            #whether we deleted rows or not, we increase the counter of the rows.
            currentRowIndexAccountingForDeletions = currentRowIndexAccountingForDeletions + 1
            
#This class function converts the XYXY data to an XYYY format
    def FromXYXYtoXYYY(self):
        masslists = [] #future lists must be must empty here to append in the for loops
        relativeintensitieslists = [] #future list
        #this loops gathers all the mass fragment numbers for each molecule in one list of arrays, while a second
        #list is made, gathering the relative intensities so that they were indexed the same as their mass fragment
        #numbers in the other list
        #this for loop grabs stuff from the reference array, whose orientation and identity is shown in the flow chart arrays document
        for referenceBy2Index in range(0,len(self.provided_reference_patterns[0,:]),2):#array-indexed for loop, only gets every other value, as half the indexes are mass lists, and the other half are relative intensity
            masslists.append(self.provided_reference_patterns[:,referenceBy2Index])#these are lists of arrays
            relativeintensitieslists.append(self.provided_reference_patterns[:,referenceBy2Index+1])#the relative intensities are after every counter, so there is a +1 (it is array indexed so since the first column is a mass list all the +1's are relative intensities)
        masslist = [] #future list
        #This for loop gets all of the mass fragments from the first index of the list, basically by not adding the 
        #'nan's or empty spaces after the numbers
        for referenceIndex in range(len(self.provided_mass_fragments)): #array-indexed for loop
            if str(masslists[0][referenceIndex]) != 'nan': #we do not want nan's in our array, the genfromtxt function calls empty boxes in excel (might be in .csv as well)'nan'.
                masslist.append(masslists[0][referenceIndex])
        #this second nested for loop gathers all the other mass fragment numbers that have not already been added to
        #the masslist, basically obtaining all the masses in the reference data and then after the loop they are sorted
        #using .sort, then an empty array of zeros is made to accommodate the output array
        for masslistIndex in range(1,len(masslists)):#array-indexed for loop, starts at one because it's checking against all arrays besides itself
            for referenceIndex in range(len(self.provided_mass_fragments)):#array-indexed for loop
                if str(masslists[masslistIndex][referenceIndex]) != 'nan':
                    if sum(masslists[masslistIndex][referenceIndex] == numpy.array(masslist)) == 0:#if the value being looked at is not equal to anything in our masslist already
                        masslist.append(masslists[masslistIndex][referenceIndex])
        masslist.sort()#puts the list in order
        reference_holder = numpy.zeros([len(masslist),len(self.provided_reference_patterns[0,:])/2+1])#makes an array that is full of zeros to hold our future reference array
        reference_holder[:,0:1] = numpy.vstack(numpy.array(masslist))#This puts the mass list in the first column of our new reference array
        #Finally, the for loop below makes a list each revolution, comparing each list of mass fragments (for each molecule)
        #and adding the relative intensities (from the identically indexed array) when they numbers were equal, and otherwise
        #adding a zero in its place. It then adds this list to the array (using numpy.vstack and numpy.array)
        for massListsIndex in range(len(masslists)):#array-indexed for loop
            relativeintensitieslist = [] #empties the list every loop
            for massListIndex in range(len(masslist)):
                placeholder = 0 #after the next for loop finishes, this is reset
                for specificMassListIndex in range(len(masslists[massListsIndex])):#array-indexed for loop, each column of .csv file being checked
                    if masslists[massListsIndex][specificMassListIndex] == masslist[massListIndex]:#This is when the index for the correct mass fragment is found
                        relativeintensitieslist.append(relativeintensitieslists[massListsIndex][specificMassListIndex])#relative intensities lists index the same way
                        placeholder = 1 #so that the next if statement will know that this has happened
                    if specificMassListIndex == len(masslists[massListsIndex])-1 and placeholder == 0:#If it comes to the end of the for loop, and there's no match, then the relative intensity is zero
                        relativeintensitieslist.append(0)
                if massListIndex == len(masslist)-1:#Once the larger for loop is done the 
                    reference_holder[:,(massListsIndex+1):(massListsIndex+2)] = numpy.vstack(numpy.array(relativeintensitieslist)) #the list is made into an array and then stacked (transposed)
        self.provided_reference_patterns = reference_holder

#populateIonizationEfficiencies is an MSReference function that populates a variable, ionizationEfficienciesList, that contains the ionization factors used in CorrectionValuesObtain
#If the ionization factor is known and in the reference data, then that value is used
#If the ionization factor is unknown the the function will look in the MID Dictionary and check if the molecule exists in the ionization data.  If it does the the ionization average of the ionization factors for that particular molecule in the data is used
#If the ionization factor is unknown and the particular molecule does not exist in the MID Data, then the function checks the molecule's ionization type(s).  The function will take all molecules from the MID data that have the same type and will perform a linear fit on the data.  The ionization factor for this molecule is determined based on the linear fit and number of electrons
#If the ionization factor is unknown, the molecule does not exist in the MID data, and the molecule's ionization type is unknown, then the function defaults to the Madix and Ko equation
    def populateIonizationEfficiencies(self, AllMID_ObjectsDict={}):
        self.ionizationEfficienciesList = numpy.zeros(len(self.molecules)) #initialize an array the same length as the number of molecules that will be populated here and used in CorrectionValuesObtain
        self.ionizationEfficienciesSourcesList = copy.copy(self.molecules) #initialize an array to store which method was used to obtain a molecule's ionization factor
        for moleculeIndex in range(len(self.molecules)): #loop through our initialized array
            if isinstance(self.relativeIonizationEfficiencies[moleculeIndex],float): #if the knownIonizationFactor is a float, then that is the value defined by the user
                self.ionizationEfficienciesList[moleculeIndex] = self.relativeIonizationEfficiencies[moleculeIndex]
                self.ionizationEfficienciesSourcesList[moleculeIndex] = 'knownIonizationFactorFromReferenceFile' #the molecule's factor was known
            else: #Ionization factor is not known so look at molecular ionization data from literatiure 
                #Initialize three lists
                MatchingMID_Objects = []
                MatchingMID_RS_Values = []
                MatchingMID_ElectronNumbers = []
                #Initialize a flag to overwrite if a molecule is in both self.molecules and the MID_ObjectDict
                matchingMolecule = False
                if self.relativeIonizationEfficiencies[moleculeIndex] == None or self.relativeIonizationEfficiencies[moleculeIndex] == 'unknown': #the ionization factor is not known so look in the AllMID_ObjectsDict
                    currentMoleculeKeyName = self.molecules[moleculeIndex] + '_IE' #holder variable to store the current molecule name + '_IE' to match the keys of AllMID_ObjectsDict (e.g. Acetaldehyde_IE)
                    if currentMoleculeKeyName in AllMID_ObjectsDict.keys(): #If the current molecule is in the dictionary, use its RS_Value
                        MatchingMID_Objects.append(currentMoleculeKeyName) #append the key
                        MatchingMID_RS_Values.append(numpy.mean(AllMID_ObjectsDict[currentMoleculeKeyName].RS_ValuesList)) #append the average RS_Value
                        MatchingMID_ElectronNumbers.append(AllMID_ObjectsDict[currentMoleculeKeyName].electronNumber) #append the electron number
                        matchingMolecule = True #set the flag to be true
                if matchingMolecule == True: #If the molecule matches a molecule in the MID dictionary, use the average RS_Value
                    self.ionizationEfficienciesList[moleculeIndex] = MatchingMID_RS_Values[0]
                    self.ionizationEfficienciesSourcesList[moleculeIndex] = 'knownIonizationFactorFromProvidedCSV' #A molecule in the reference data is also in the ionization data
                elif matchingMolecule == False: #Otherwise matchingMolecule is False which means its not in the data from literature.  So we will approximate the ionization factor based on a linear fit of the data from literature that share the molecule's type or use the Madix and Ko equation
                    if self.moleculeIonizationType[moleculeIndex] != None and self.moleculeIonizationType[moleculeIndex] != 'unknown': #IF the user did not manually input the ionization factor and none of the molecules in the MID_Dict matched the current molecule
                        #Then get an estimate by performing a linear fit on the data in the MID Dictionary
			#TODO:The program currently only takes in one type but it is a desired feature to allow users to put in multiple types such as type1+type2 which would make a linear fit of the combined data between the two types
			#TODO continued:The user should also be able to put in type1;type2 and the program would find the ionization factor using a linear fit of data from type1 and using a linear fit of data from type2.  The largest of the two ionization factors would be used.
			#TODO continued:Then doing type1+type2;type3 would take the larger value between the linear fit of the combined type1 and type2 data or the value from the linear fit of type3 data
                        for key in AllMID_ObjectsDict: #Loop through the MID Dictionary
                            for MID_MoleculeType in AllMID_ObjectsDict[key].moleculeIonizationType: #Loop through the ionization types to get all the ionization types of a particular molecule (e.g. Ethanol is both an alcohol and a hydrogen non-metal-ide so its RS value(s) will be included if the user has a molecule that is either an alcohol or a hydrogen non-metal-ide)
                                #Use stringCompare to check if a molecule in the MID Dictionary matches a molecule in the reference data since casing and spacing may differ between the two (e.g. reference data may have carbon dioxide while MID Dictionary may have Carbon Dioxide)
                                if parse.stringCompare(self.moleculeIonizationType[moleculeIndex],MID_MoleculeType): #If the knownMoleculeType matches an MID object's molecule type
                                    MatchingMID_Objects.append(key) #Append the key
                                    MatchingMID_RS_Values.append(numpy.mean(AllMID_ObjectsDict[key].RS_ValuesList)) #Append the average of the RS values
                                    MatchingMID_ElectronNumbers.append(AllMID_ObjectsDict[key].electronNumber) #append the electron number
                        if len(MatchingMID_Objects) == 1: #If only one data point add (0,0) to the data
                            MatchingMID_Objects.append(MatchingMID_Objects[0]) #append the molecule type to the objects list
                            MatchingMID_RS_Values.append(0.0)
                            MatchingMID_ElectronNumbers.append(0.0)
                        if len(MatchingMID_Objects) > 1: #When we have more than one value in the data, find a linear fit
                            #Now we use polyfit, poly1d, and polyval to fit the data linearly and find an approximate ionization factor
                            #TODO: We think we should use a power law with y = mx^b (this implies an intercept of 0 and retains the type of curvature we see in the data)
                            #TODO continued: Link to do so: https://scipy-cookbook.readthedocs.io/items/FittingData.html
                            polynomialCoefficients = numpy.polyfit(MatchingMID_ElectronNumbers,MatchingMID_RS_Values,1) #Electron numbers as the independent var, RS_values as the dependent var, and 1 for 1st degree polynomial
                            poly1dObject = numpy.poly1d(polynomialCoefficients) #create the poly1d object
                            self.ionizationEfficienciesList[moleculeIndex] = numpy.polyval(poly1dObject,self.electronnumbers[moleculeIndex]) #use polyval to calculate the ionization factor based on the current molecule's electron number
                            self.ionizationEfficienciesSourcesList[moleculeIndex] = 'evaluatedInterpolationTypeFit' #ionization factor was determined via a linear fit based on a molecule's ionization type
                if len(MatchingMID_Objects) == 0: #Otherwise use the original Madix and Ko equation
                    self.ionizationEfficienciesList[moleculeIndex] = (0.6*self.electronnumbers[moleculeIndex]/14)+0.4        
                    self.ionizationEfficienciesSourcesList[moleculeIndex] = 'MadixAndKo' #ionization efficiency obtained via Madix and Ko equation

#Export the ionization efficiencies used and their respective method used to obtain them (known factor, known molecule, known ionization type, or Madix and Ko)    
    def exportIonizationInfo(self):
        ionizationData = numpy.vstack((self.molecules,self.ionizationEfficienciesList,self.ionizationEfficienciesSourcesList)) #make a 2d array containing molecule names (for the header), the ionization efficiencies, and which method was chosen
        ionizationDataAbsicca = numpy.array([['Molecule'],
                                             ['Ionization Efficiency'],
                                             ['Method to Obtain Ionization Efficiency']]) #create the abscissa headers for the csv file
        ionizationDataToExport = numpy.hstack((ionizationDataAbsicca,ionizationData)) #use hstack to obtain a 2d array with the first column being the abscissa headers
        numpy.savetxt('ExportedIonizationEfficienciesSourcesTypes.csv',ionizationDataToExport,delimiter=',',fmt='%s') #export to a csv file

#readReferenceFile is a helper function that reads the reference file in a certain form and returns the
#variables and data that are used to initialize the class. It can read files both in XYYY and XYXY form.
def readReferenceFile(referenceFileName, form):        
     #This function converts the XYXY data to an XYYY format
    def FromXYXYtoXYYY(provided_reference_patterns):
        print("Warning: FromXYXYtoXYYY for converting data patterns has not been tested in a long time. A unit test should be created and checked prior to use. Then this warning updated (this warning appears in two parts of the code." )
        masslists = [] #future lists must be must empty here to append in the for loops
        relativeintensitieslists = [] #future list
        #this loops gathers all the mass fragment numbers for each molecule in one list of arrays, while a second
        #list is made, gathering the relative intensities so that they were indexed the same as their mass fragment
        #numbers in the other list
        #this for loop grabs stuff from the reference array, whose orientation and identity is shown in the flow chart arrays document
        for referenceBy2Index in range(0,len(provided_reference_patterns[0,:]),2):#array-indexed for loop, only gets every other value, as half the indexes are mass lists, and the other half are relative intensity
            masslists.append(provided_reference_patterns[:,referenceBy2Index])#these are lists of arrays
            relativeintensitieslists.append(provided_reference_patterns[:,referenceBy2Index+1])#the relative intensities are after every counter, so there is a +1 (it is array indexed so since the first column is a mass list all the +1's are relative intensities)
        masslist = [] #future list
        #This for loop gets all of the mass fragments from the first index of the list, basically by not adding the 
        #'nan's or empty spaces after the numbers
        provided_mass_fragments = provided_reference_patterns[:,0] 
        for referenceIndex in range(len(provided_mass_fragments)): #array-indexed for loop
            if str(masslists[0][referenceIndex]) != 'nan': #we do not want nan's in our array, the genfromtxt function calls empty boxes in excel (might be in .csv as well)'nan'.
                masslist.append(masslists[0][referenceIndex])
        #this second nested for loop gathers all the other mass fragment numbers that have not already been added to
        #the masslist, basically obtaining all the masses in the reference data and then after the loop they are sorted
        #using .sort, then an empty array of zeros is made to accommodate the output array
        for masslistIndex in range(1,len(masslists)):#array-indexed for loop, starts at one because it's checking against all arrays besides itself
            for referenceIndex in range(len(provided_mass_fragments)):#array-indexed for loop
                if str(masslists[masslistIndex][referenceIndex]) != 'nan':
                    if sum(masslists[masslistIndex][referenceIndex] == numpy.array(masslist)) == 0:#if the value being looked at is not equal to anything in our masslist already
                        masslist.append(masslists[masslistIndex][referenceIndex])
        masslist.sort()#puts the list in order
        reference_holder = numpy.zeros([len(masslist),len(provided_reference_patterns[0,:])/2+1])#makes an array that is full of zeros to hold our future reference array
        reference_holder[:,0:1] = numpy.vstack(numpy.array(masslist))#This puts the mass list in the first column of our new reference array
        #Finally, the for loop below makes a list each revolution, comparing each list of mass fragments (for each molecule)
        #and adding the relative intensities (from the identically indexed array) when they numbers were equal, and otherwise
        #adding a zero in its place. It then adds this list to the array (using numpy.vstack and numpy.array)
        for massListsIndex in range(len(masslists)):#array-indexed for loop
            relativeintensitieslist = [] #empties the list every loop
            for massListIndex in range(len(masslist)):
                placeholder = 0 #after the next for loop finishes, this is reset
                for specificMassListIndex in range(len(masslists[massListsIndex])):#array-indexed for loop, each column of .csv file being checked
                    if masslists[massListsIndex][specificMassListIndex] == masslist[massListIndex]:#This is when the index for the correct mass fragment is found
                        relativeintensitieslist.append(relativeintensitieslists[massListsIndex][specificMassListIndex])#relative intensities lists index the same way
                        placeholder = 1 #so that the next if statement will know that this has happened
                    if specificMassListIndex == len(masslists[massListsIndex])-1 and placeholder == 0:#If it comes to the end of the for loop, and there's no match, then the relative intensity is zero
                        relativeintensitieslist.append(0)
                if massListIndex == len(masslist)-1:#Once the larger for loop is done the 
                    reference_holder[:,(massListsIndex+1):(massListsIndex+2)] = numpy.vstack(numpy.array(relativeintensitieslist)) #the list is made into an array and then stacked (transposed)
        provided_reference_patterns = reference_holder
        return provided_reference_patterns
    
    #read the csv file into a dataframe
    if '.csv' in referenceFileName:
        dataFrame = pandas.read_csv('%s' %referenceFileName, header = None)
    elif '.tsv' in referenceFileName:
        try: #no easy way to assess utf16 vs utf8, so try both.
            dataFrame = pandas.read_csv('%s' %referenceFileName, header = None, delimiter = '\t', encoding = 'utf8') #need to specify encoding for cases of tab delimited files.
        except: #no easy way to assess utf16 vs utf8, so try both.
            dataFrame = pandas.read_csv('%s' %referenceFileName, header = None, delimiter = '\t', encoding = 'utf16') #need to use utf16 for some cases of tab delimited files.
    
    if form == 'xyyy':
        for rowIndex in range(len(dataFrame)): #Loop through each row and check the abscissa value
            try: #Try to convert the abscissa title to a float
                float(dataFrame.iloc[rowIndex][0]) #if successful, then this rowIndex is the first index of provided reference intensities
                dfreference = dataFrame.iloc[rowIndex:][:] #remove the rows of headers
                reference = dfreference.values #convert to matrix
                provided_reference_patterns = reference.astype(float) #convert the matrix to floats
                provided_reference_patterns = DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(provided_reference_patterns,startingRowIndex=1) #clear row of zeros
                break #exit the for loop
            except: #Otherwise the row consists of other information
                if (dataFrame.iloc[rowIndex][0] == 'SourceOfFragmentationPatterns') or (dataFrame.iloc[rowIndex][0] == 'Source:'): #if the abscissa titles the source (both old and new reference files)
                    dfSourceOfFragmentationPatterns = dataFrame.iloc[rowIndex][1:] #select the row of names
                    SourceOfFragmentationPatterns = dfSourceOfFragmentationPatterns.values #convert to matrix
                    SourceOfFragmentationPatterns = SourceOfFragmentationPatterns.astype(str) #save as class object with type string
                elif dataFrame.iloc[rowIndex][0] == 'sourceOfIonizationData':
                    dfsourceOfIonizationData = dataFrame.iloc[rowIndex][1:] #Select the row of names
                    sourceOfIonizationData = dfsourceOfIonizationData.values #convert to matrix
                    sourceOfIonizationData = sourceOfIonizationData.astype(str) #save as class object with type string
                elif dataFrame.iloc[rowIndex][0] == 'Molecules': #if the abscissa titles the molecule names
                    dfmolecules = dataFrame.iloc[rowIndex][1:] #select the row of names
                    molecules = dfmolecules.values #convert to matrix
                    molecules = molecules.astype(str) #save as class object with type string
                elif dataFrame.iloc[rowIndex][0] == 'Electron Numbers': #if the abscissa titles the electron numbers
                    dfelectronnumbers = dataFrame.iloc[rowIndex][1:] #select the row of names
                    electronnumbers = dfelectronnumbers.values #convert to matrix
                    electronnumbers = electronnumbers.astype(int) #save as class object with type int
                elif dataFrame.iloc[rowIndex][0] == 'Molecular Mass': #if the abscissa titles the molecular weights
                    dfmolecularWeights = dataFrame.iloc[rowIndex][1:] #select row of names
                    molecularWeights = dfmolecularWeights.values #convert to matrix
                    molecularWeights = molecularWeights.astype(numpy.double) #save as class object with type double
                elif dataFrame.iloc[rowIndex][0] == 'moleculeIonizationType':
                    dfmoleculeIonizationType = dataFrame.iloc[rowIndex][1:] #select row of names
                    moleculeIonizationType = dfmoleculeIonizationType.values #convert to matrix
                    moleculeIonizationType = moleculeIonizationType.astype(str) #save as class object with type string
                elif dataFrame.iloc[rowIndex][0] == 'relativeIonizationEfficiencies':
                    dfrelativeIonizationEfficiencies = dataFrame.iloc[rowIndex][1:] #select row of names
                    relativeIonizationEfficiencies = dfrelativeIonizationEfficiencies.values #convert to matrix
                    for index in range(len(relativeIonizationEfficiencies)):
                        try: #try to convert to a float
                            relativeIonizationEfficiencies[index] = float(relativeIonizationEfficiencies[index])
                        except: #if not possible, the value is probably None or 'unknown' so leave as a string
                            pass
#                    relativeIonizationEfficiencies = relativeIonizationEfficiencies.astype(float) #save as class object with type float
        
        
        sourceOfIonizationData = None #To remove MSRESOLVE dependencies.
        relativeIonizationEfficiencies = None  #To remove MSRESOLVE dependencies.
        moleculeIonizationType = None
        '''
        Below is removed to remove external dependencies.
        try: #if using an older reference file, it will not have ionization factors so the elif statement never gets entered meaning knownIonizationFactors does not exist
            relativeIonizationEfficiencies #Try calling this variable, if it exists there will be no error
        except: #if it does not exist, populate it with unknown
            relativeIonizationEfficiencies = ['unknown'] #initialize as a list of len(1)
            relativeIonizationEfficiencies = parse.parallelVectorize(relativeIonizationEfficiencies,len(molecules)) #parallel vectorize to length of molecules
            relativeIonizationEfficiencies = numpy.array(relativeIonizationEfficiencies) #convert to matrix
            
        try: #if using an old reference file, it will not have ionization types so the elif statement never gets entered meaning moleculeIonizationType does not exist
            moleculeIonizationType #Try calling this variable, if it exists there will be no error
        except: #if it does not exist, populate it with unknown
            moleculeIonizationType = ['unknown'] #initialize as a list of len(1)
            moleculeIonizationType = parse.parallelVectorize(moleculeIonizationType,len(molecules)) #parallel vectorize to length of molecules
            moleculeIonizationType = numpy.array(moleculeIonizationType) #convert to matrix
        
        try: #If using an older reference file, it will not have SourceOfIonizationInfo so the elif statement never gets entered meaning the variable does not exist
            sourceOfIonizationData #try calling the variable, if it exists there will be no error
        except: #If it does not exist, populate with empty strings
            sourceOfIonizationData = [''] #initialize as a list of len(1)
            sourceOfIonizationData = parse.parallelVectorize(sourceOfIonizationData,len(molecules)) #parallel vectorize to length of molecules
            sourceOfIonizationData = numpy.array(sourceOfIonizationData) #convert to matrix
        '''  
                
#        ''' generate reference matrix'''
#        #remove top 4 rows
#        dfreference = dataFrame.iloc[4:][:]
#        #convert to matrix
#        reference = dfreference.values
#        #convert the matrix to floats
#        provided_reference_patterns = reference.astype(float)
#        #clear rows of zeros
#        provided_reference_patterns=DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(provided_reference_patterns,startingRowIndex=1)
#    
#        '''generate electron number list'''
#        #select row of electron numbers
#        dfelectronnumbers = dataFrame.iloc[2][1:]
#        #convert to matrix
#        electronnumbers = dfelectronnumbers.values
#        #save as class object with type int
#        electronnumbers = electronnumbers.astype(int32)
#   
#        '''generate list of molecule names'''
#        #select row of names
#        dfmolecules = dataFrame.iloc[1][1:]
#        #convert to matrix
#        molecules = dfmolecules.values
#        #save as class object with type string
#        molecules = molecules.astype(str)
#        
#        '''generate list of molecular weights'''
#        #select row of names
#        dfmolecularWeights = dataFrame.iloc[3][1:]
#        #convert to matrix
#        molecularWeights = dfmolecularWeights.values
#        #save as class object with type float
#        molecularWeights = molecularWeights.astype(float)
#        
#        '''generate list of source information'''
#        #select row of names
#        dfsourceInfo = dataFrame.iloc[0][1:]
#        #convert to matrix
#        sourceInfo = dfsourceInfo.values
#        #save as class object with type string
#        sourceInfo = sourceInfo.astype(str)
        
        '''list of massfragments monitored is not part of reference file'''
        mass_fragment_numbers_monitored = None
        
    elif form == 'xyxy':
        for rowIndex in range(len(dataFrame)): #Loop through each row and check the abscissa value
            try: #Try to convert the abscissa title to a float
                float(dataFrame.iloc[rowIndex][0]) #if successful, then this rowIndex is the first index of provided reference intensities
                dfreference = dataFrame.iloc[rowIndex:][:] #remove the rows of headers
                reference = dfreference.values #convert to matrix
                provided_reference_patterns = reference.astype(float) #convert the matrix to floats
                print("Warning: FromXYXYtoXYYY for converting data patterns has not been tested in a long time. A unit test should be created and checked prior to use. Then this warning updated (this warning appears in two parts of the code." )
                provided_reference_patterns = FromXYXYtoXYYY(provided_reference_patterns) #convert reference from XYXY to XYYY
                provided_reference_patterns = DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(provided_reference_patterns,startingRowIndex=1) #clear row of zeros
                break #exit the for loop
            except: #Otherwise the row consists of other information
                if dataFrame.iloc[rowIndex][0] == 'Source:': #if the abscissa titles the source
                    dfsourceInfo = dataFrame.iloc[rowIndex][1::2] #select the row of names
                    sourceInfo = dfsourceInfo.values #convert to matrix
                    sourceInfo = sourceInfo.astype(str) #save as class object with type string
                elif dataFrame.iloc[rowIndex][0] == 'Molecules': #if the abscissa titles the molecule names
                    dfmolecules = dataFrame.iloc[rowIndex][1::2] #select the row of names
                    molecules = dfmolecules.values #convert to matrix
                    molecules = molecules.astype(str) #save as class object with type string
                elif dataFrame.iloc[rowIndex][0] == 'Electron Numbers': #if the abscissa titles the electron numbers
                    dfelectronnumbers = dataFrame.iloc[rowIndex][1::2] #select the row of names
                    electronnumbers = dfelectronnumbers.values #convert to matrix
                    electronnumbers = electronnumbers.astype(int32) #save as class object with type int
                elif dataFrame.iloc[rowIndex][0] == 'Molecular Mass': #if the abscissa titles the molecular weights
                    dfmolecularWeights = dataFrame.iloc[rowIndex][1::2] #select row of names
                    molecularWeights = dfmolecularWeights.values #convert to matrix
                    molecularWeights = molecularWeights.astype(float) #save as class object with type float
#        '''generate reference matrix'''
#        #remove top 4 rows
#        dfreference = dataFrame.iloc[4:][:]
#        #convert to matrix
#        reference = dfreference.values
#        #convert the matrix to floats 
#        provided_reference_patterns = reference.astype(float)
#        #convert reference from XYXY to XYYY
#        print("Warning: FromXYXYtoXYYY for converting data patterns has not been tested in a long time. A unit test should be created and checked prior to use. Then this warning updated (this warning appears in two parts of the code." )
#        provided_reference_patterns=FromXYXYtoXYYY(provided_reference_patterns)
#        #clear rows of zeros
#        provided_reference_patterns=DataFunctions.removeColumnsWithAllvaluesBelowZeroOrThreshold(provided_reference_patterns,startingRowIndex=1)
#        
#        '''generate electron numbers list'''
#        #create data frame of electron numbers
#        dfelectronnumbers = dataFrame.iloc[2,1::2]
#        #convert to matrix
#        electronnumbers = dfelectronnumbers.values
#        #save as class object with type int
#        electronnumbers = electronnumbers.astype(int32)
#        
#        '''generate list of molecule names'''
#        #select matrix of names
#        dfmolecules = dataFrame.iloc[1,1::2]
#        #convert to matrix
#        molecules = dfmolecules.values
#        #save as class object with type string
#        molecules = molecules.astype(str)
#        
#        '''generate list of molecular weights'''
#        #select row of names
#        dfmolecularWeights = dataFrame.iloc[3][1::2]
#        #convert to matrix
#        molecularWeights = dfmolecularWeights.values
#        #save as class object with type float
#        molecularWeights = molecularWeights.astype(float)
#        
#        '''generate list of source information'''
#        #select row of names
#        dfsourceInfo = dataFrame.iloc[0][1::2]
#        #convert to matrix
#        sourceInfo = dfsourceInfo.values
#        #save as class object with type string
#        sourceInfo = sourceInfo.astype(str)

        '''list of massfragments monitored is not part of reference file'''
        mass_fragment_numbers_monitored = None
        
    return provided_reference_patterns, electronnumbers, molecules, molecularWeights, SourceOfFragmentationPatterns, sourceOfIonizationData, relativeIonizationEfficiencies, moleculeIonizationType, mass_fragment_numbers_monitored, referenceFileName, form

#This function operates in a parallel way to trimDataMasses, but it operates on the reference data and all of it's constituent variables  
def trimDataMoleculesToMatchChosenMolecules(ReferenceData, chosenMolecules):
    
    print("MoleculeChooser")
    #getting a list of all molecules (a header) to compare to during trimming.	    
    allMoleculesList = ReferenceData.molecules
    
    #initializing object that will become the trimmed copy of ReferenceData
    trimmedRefererenceData = copy.deepcopy(ReferenceData)
    
    #trim the reference fragmentation patterns to only the selected molecules 
    #unused trimmed copy molecules is just a place holder to dispose of a function return that is not needed
    trimmedReferenceIntensities, trimmedMoleculesList = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.provided_reference_patterns[:,1:],
                                                                                                                    allMoleculesList, chosenMolecules, header_dtype_casting=str)  
    #add a second dimension to the reference data
    trimmedReferenceMF = numpy.reshape(trimmedRefererenceData.provided_mass_fragments,(-1,1))
    
    #TODO: The below line works with provided_reference_patterns. This is because trimDataMoleculesToMatchChosenMolecules
    #TODO continued: is currently working prior to standardized Reference patterns existing, and also because it is occurring
    #TODO continued: Before we have the monitored mass fragments (which also occurs later data analysis).
    #TODO continued: The best solution is probably to do the standardization earlier and then do this trimming after that.
    #Add the abscissa back into the reference values
    trimmedRefererenceData.provided_reference_patterns = numpy.hstack((trimmedReferenceMF,trimmedReferenceIntensities))
    
    #Shorten the electronnumbers to the correct values, using the full copy of molecules. Do the same for molecularWeights and sourceInfo
    trimmedRefererenceData.electronnumbers, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.electronnumbers, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    trimmedRefererenceData.molecularWeights, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.molecularWeights, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    trimmedRefererenceData.moleculeIonizationType, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.moleculeIonizationType, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    trimmedRefererenceData.relativeIonizationEfficiencies, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.relativeIonizationEfficiencies, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    trimmedRefererenceData.SourceOfFragmentationPatterns, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.SourceOfFragmentationPatterns, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    trimmedRefererenceData.sourceOfIonizationData, trimmedMoleculesList  = DataFunctions.KeepOnlySelectedYYYYColumns(trimmedRefererenceData.sourceOfIonizationData, allMoleculesList, chosenMolecules, Array1D = True, header_dtype_casting=str)
    
    trimmedRefererenceData.molecules = trimmedMoleculesList
    
    #remove any zero rows that may have been created
    trimmedRefererenceData.ClearZeroRowsFromProvidedReferenceIntensities()
    
    #!!!trimmedRefererenceData.ExportCollector("MoleculeChooser", use_provided_reference_patterns=True)    #this line commented out o remove external depenendencies
    return trimmedRefererenceData

    
'''
Standardize is a simple algorithim that reads in one numpy array
(oneArray) and scales the whole array so that the
greatest number is 100 within the array
'''

def StandardizeTo100(a1DArray,n):
    maxNumber=float(numpy.max(a1DArray))
    multiplier=100/maxNumber
    for index2 in range(0,len(a1DArray)):
        a1DArray[index2]=(a1DArray[index2]*multiplier)
    return a1DArray



'''
StandardizeReferencePattern uses StandardizeTo100 to modify reference values so that each column is scaled to
100. NOTE THAT THIS FUNCTION ASSUMES THAT THE FIRST COLUMN in reference contains the 
mass framgnet numbers. 
Parameters: 
standardizedReference -  a name chosen for the numpy array that contains reference values
num_of_molecues-  an integer describing the number of  molecues that contributed to the reference file
'''
def StandardizeReferencePattern(referenceUnstandardized,num_of_molecules):
    # preallocate new array for standardized values
    standardizedReference = copy.deepcopy(referenceUnstandardized)

    # standardize
    for moleculeIndex in range(1,num_of_molecules+1):
        standardizedReference[0:,moleculeIndex]=StandardizeTo100(referenceUnstandardized[0:,moleculeIndex],1)

    return standardizedReference
    
#this function finds the significance of a specified value in an array to the array as a whole
def IndElemSignificanceCalculator(rowDataArray, specifiedColumnIndex, moleculesLikelihood, minThreshold=None, maxIntensityPossible=100):
    length = len(rowDataArray)
    #variable to hold the terms in the summation
    allSummationTerms = []
    #for each value in the array

    if minThreshold == 0: minThreshold=None #the code is written to have minThreshold as None or not None, but 0 counts as None  for this purpose.
    if minThreshold != None: #assume the minThreshold is some kind of number, and use it to calculate and cap the indSummationTerm.
        indSummationTermCap = maxIntensityPossible/minThreshold - 1 #note that the 100 is the default maxIntensityPossible, with the assumption that standardized intensities  are capped at 100.
    
    for moleculecounter in range(length):
        if rowDataArray[moleculecounter] != 0: #if the value is not zero, calculate 
            if minThreshold == None: #just calculate the indSummationTerm as normal.
                #calculates the unweighted ratio of each value, scaled by the likelihood of that molecule 
                indSummationTerm = abs((moleculesLikelihood[moleculecounter]*rowDataArray[moleculecounter])**float(-1)*(moleculesLikelihood[specifiedColumnIndex]*rowDataArray[specifiedColumnIndex]-1))
                allSummationTerms.append(indSummationTerm)       
            if minThreshold != None: #assume the minThreshold is some kind of number, and use it to calculate and cap the indSummationTerm.
                indSummationTerm = abs((moleculesLikelihood[moleculecounter]*rowDataArray[moleculecounter])**float(-1)*(moleculesLikelihood[specifiedColumnIndex]*rowDataArray[specifiedColumnIndex]-1))
                if indSummationTerm > indSummationTermCap: indSummationTerm = indSummationTermCap
                allSummationTerms.append(indSummationTerm)
        if rowDataArray[moleculecounter] == 0: 
            if minThreshold == None: 
                pass #the pass is like allSummationTerms.append(0) #if the intensity/value is zero, and we don't have a minThreshold to use, then the final value will be zero as well.
            if minThreshold != None: #else assume the minThreshold is some kind of number.
                indSummationTerm = indSummationTermCap #use the cap directly if the term for rowDataArray[moleculecounter] == 0, because that means a denominator of 0 which is like infinity. 
                allSummationTerms.append(indSummationTerm)            
    #the following line can be replace with code such as "significance = (sum(allSummationTerms)**SumCoeffient)*(array[specifiedColumnIndex]**ValueCoefficent)"
    # if you would like to add coefficents to increase or decrease the weighting of each term
    significance = sum(allSummationTerms)*rowDataArray[specifiedColumnIndex]*moleculesLikelihood[specifiedColumnIndex]
    return significance

#This function compiles a list of the significances of each row to a particular column 
def ElemSignificanceCalculator(anArray,specifiedColumnIndex, moleculesLikelihood, minThreshold=None, maxIntensityPossible=100):
    #find the number of rows
    row_num = len(anArray)
    #empty list to store values
    sigValuesList = []
    #for each row...
    for rowcounter in range(row_num):
        # the "Significance" of that row to the column is calculated
        sigValue = IndElemSignificanceCalculator(anArray[rowcounter], specifiedColumnIndex, moleculesLikelihood, minThreshold=minThreshold, maxIntensityPossible=maxIntensityPossible)
        # the significance is stored in a list
        sigValuesList.append(sigValue)
        
    return sigValuesList

#This function is made to conduct prelimiary checks. The checks are designed to
#fail only if there is not chance of the the chosen reference data passing the
#SLS method. There are 2 cases where this occurs :1) a molecule does not 
#contain reference intensities for any mass fragments in the combination or 2)
#The reference data for the mass fragments does not contain any zeros.
def passesRowsSumChecks(rowSumsList, massFragCombination, allOverlappingPatterns, numberOfMassFragments=None, numberOfRows=None):
    if numberOfMassFragments==None:
        numberOfMassFragments=len(massFragCombination)##It is slightly faster to pass in the length variable to make it shorter      
    passesRowsSumChecks=True#Initialize return variable as true
    if 0 in rowSumsList: #Check if any row is entirely full of zeros
        passesRowsSumChecks=False
    elif numpy.array(rowSumsList).all() == numberOfMassFragments: #Check if the array passed to it is full. If it is full, it appends to the allOverlapping patterns list
        allOverlappingPatterns.append(massFragCombination)###May be an error with the use of allOverlapping patterns
        #FIXME: #TODO: should be printed to a file.
        #print(allOverlappingPatterns)
        passesRowsSumChecks=False

    return passesRowsSumChecks, allOverlappingPatterns #Return true if the rowsSumsList passes the check

#The function maintains two lists: 1)that contains the objective function values that 
#need to be kept in a sorted order 2) a parallel list where a value needs to be inserted 
#according to the sorting of the first one. It also takes in an integer value, N, that limits 
#the the length of the two lists to N values. Using a binary search method, it will 
#find the location where the value to insert will be inserted(if possible). The
#value will be inserted there and the last value removed from the list (if
#applicable).  The bisect method used requires the list be presorted in ascending order.
#Thus, the algorithm here inserts values in a way to create an ascending ordered list.

#By default storeAndPop is used to keep the best N values based on minimizing
#the objective function values. If instead it is desirable to retain values 
#with the objective function maximized, the optional argument 'optimumType'
#should be set to `Maximum`.

#Alternatively, the objective function values could be multiplied by -1.
#The function supports multidimensional objective functions in nested objects
#such as tuples or lists, e.g., (A,B,C) in which case it will be sorted
#by A, then B, then C.  For multidimensional objective functions, it is
#necessary to have the dimensions as either all "maximum" or all "minimum" type.
#Using a -1 factor can be helpful in this regard (e.g., passing in (A,-B,C) etc.

#If the values in the sample space for parallelList are not unique it is 
#possible that this repeated calls to this function could lead to 
#a parallelList of a particular value repeated many times. If repeated values
#are undesired, then excludeDuplicates can be set to False.
def storeAndPop(objectiveFunctionValuesList, objectiveFunctionValueToInsert, 
                parallelList, valueToInsertInParallelList, maxItemsAllowed,
                optimumType="Minimum", excludeDuplicates=True):
    
    #Find the insertion index where the value will be inserted by using a binary
    #search
    insertionIndex=bisect.bisect(objectiveFunctionValuesList,
                                 objectiveFunctionValueToInsert)

    #Initialize a variable to keep track of if a value was inserted into the
    #list.
    valueStoredInList=False
    
    #If it is a duplicate exit now without checking everything else
    #Note that we only check the value in parallel list to the left of
    #the insertion index. This is because bisect.bisect() will specify 
    #and insertion index such that a duplicate is inserted
    #to the right of the original.
    if (excludeDuplicates and len(parallelList) and 
        (parallelList[insertionIndex-1] == valueToInsertInParallelList)):
        
        #This value is a duplicate. Return the original lists.
        return objectiveFunctionValuesList, parallelList, valueStoredInList    
    
    #If the list isn't yet filled, the value will inherently be in the top N
    #value in the list. This value can just be inserted at the insertionIndex.
    if len(objectiveFunctionValuesList)<maxItemsAllowed:    
        objectiveFunctionValuesList.insert(insertionIndex,
                                           objectiveFunctionValueToInsert)
        parallelList.insert(insertionIndex, valueToInsertInParallelList)
        valueStoredInList=True
    #If the list already contains N elements, a new element could either be 
    #inserted in the list or at the end of the list. For optimumType == Minimum,
    #an item at the end would be worse and thus nothing should be added. 
    #However, for optimumType == Maximum an object at the end would be the best
    #and thus should be added.
    elif (len(objectiveFunctionValuesList) == maxItemsAllowed and 
            (insertionIndex<maxItemsAllowed or optimumType=="Maximum")):
        #insert the value to insert in the location found through the binary
        #search
        objectiveFunctionValuesList.insert(insertionIndex,
                                           objectiveFunctionValueToInsert)
        parallelList.insert(insertionIndex, valueToInsertInParallelList)
        valueStoredInList=True
        
        if optimumType == 'Minimum':
            #delete the last element since something was added to the list
            del objectiveFunctionValuesList[-1]
            del parallelList[-1]
        elif optimumType == 'Maximum':
            #delete the first element since something was added to the list
            del objectiveFunctionValuesList[0]
            del parallelList[0]   
        else:
            raise  ValueError("optimumType must be either 'Maximum' " +
                              "or 'Minimum'")

    return objectiveFunctionValuesList, parallelList, valueStoredInList

'''This function exports all XYYY Data to the User specified document (usually a CSV)'''
#TODO Future Development: This function could benefit from creating folders to store the different
#runs of the program if it is run in the same directory as before
## NOTE: This function leaves a tailing comma on each line of the
## created csv file. When pandas is used to read the csv this
## will result in the creation of a column of 'nan'
## ImportWorkingData() has been modified to remove this column
## If the data header and the data are the same length, the abscissa header is disregarded.
def ExportXYYYData(outputFileName, data, dataHeader, abscissaHeader = 'Mass', fileSuffix = '', dataType = None, rowIndex = [], units = None): 
    formatedDataHeader = dataHeader
    if dataType == 'preProcessed' or dataType == 'simulated' or dataType == 'Experiment':
        formatedDataHeader = ['m%s' % MFNumber for MFNumber in dataHeader]
    if dataType == 'scaled':
        formatedDataHeader = ['%s Concentration Relative to CO' % molecule for molecule in dataHeader]
    if dataType == 'concentration':
        label = ' Concentration(%s)' % units
        formatedDataHeader = [molecule + label for molecule in dataHeader]
    #extraLine is used to create CSV files that conform to MSRESOLVE's import requirements i.e. having a row for comments at the top
    extraLine = False
    if dataType == 'Experiment':
        extraLine = len(data[0,1:])
        
#If future applications of Export XYYY are desired, the new formats can be 
#specified by additional keywords and if statements.

#if iterative analysis is being used and the suffix is wanted
    if not fileSuffix =='':
        #then the filename will have a suffix attached
        outputFileName = outputFileName[:-4] + fileSuffix + outputFileName[-4:]

    #testing if file is open, and rename if it is
    #create list of name options
    nameOptions = [''] + list(range(100))
    for x in nameOptions:
        # create new name
        filename = (outputFileName[:-4] + '%s' + outputFileName[-4:]) %x
        #Test if it can be opened
        try:
            open(filename, 'w')  
            break
        except(IOError):
            pass
    #combine the column headers and data into one array
    try:
        fullArrayToExport = numpy.vstack((formatedDataHeader,data))
    #occasionally, abscissaHeader needs to be inserted without rowIndex being used
    except ValueError: 
        formatedDataHeader = numpy.hstack((abscissaHeader,formatedDataHeader))
        fullArrayToExport = numpy.vstack((formatedDataHeader,data))
        
    #if the row index isn't included in the data, then add it 
    if rowIndex != []:    
        abscissaHeader = numpy.transpose((numpy.array([[abscissaHeader]])))
        rowIndex = numpy.transpose([rowIndex])
        abscissaArrayToExport =  numpy.vstack((abscissaHeader,rowIndex))
        fullArrayToExport = numpy.hstack((abscissaArrayToExport,fullArrayToExport))
    #insert an extra line with a header of the data type. Included to allow exported files to be uploaded during iterative analysis.
    if not extraLine == False:
        lineToInsert = "%s,%s" %(dataType, ',' * (extraLine))
        lineToInsert = numpy.array(lineToInsert.split(','))
        fullArrayToExport = numpy.vstack((lineToInsert, fullArrayToExport))
    #save the file to the correct name
    numpy.savetxt(filename, fullArrayToExport, delimiter = ',', fmt ="%s")
  
