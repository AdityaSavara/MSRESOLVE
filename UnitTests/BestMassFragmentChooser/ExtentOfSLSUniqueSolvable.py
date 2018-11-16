import sys
import os
baseDir = os.getcwd()
sys.path.insert(1, os.path.join(baseDir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(baseDir, os.pardir))
sys.path.insert(1, os.path.join(baseDir, os.pardir, os.pardir))

import numpy as np
import MSRESOLVE
import bisect

def ExtentOfSLSUniqueSolvable(massFrags, species, 
    onesAndZeros, significanceMatrix):
    """
    Arguments:
        massFrags - An n by 1 numpy array that holds the mass fragment 
                    masses. This is an index on the vertical axis of 
                    the onesAndZeros matrix.
    
        species - An m by 1 numpy array that holds the names of the
                  species being considered as strings. This is an index
                  for the horizontal axis of the onesAndZeros matrix.
                    
        onesAndZeros - This is the n by m 2-d numpy array of one and zeros.
                       Each row corresponds to a particular mass fragment
                       and each column to a species. Should be 
                       cleaned by this point to remove irrelevant
                       mass fragment rows.
                       
        significanceMatrix - A n by m 2-d numpy array that mirrors onesAndZeros.
                             The elements of this matrix indicate the 
                             significance of the corresponding intensity. This 
                             is used to select the best unique value/row in 
                             situations where more than one are available.
                       
    Returns a 3-tuple (extentSolvable, solvedSpecies, massFragsUsed):
    
        extentSolvable - An integer that indicates how many species were 
                         solved through unique SLS. It is the primary 
                         metric to judge the suitability of a particular
                         mass fragment combination
                         
        summedSignificance - An scalar that is the sum of the significances 
                             of the intensities used for SLS. It serves
                             as a secondary metric for selecting mass 
                             fragement combinations. 
                             
        
        solvedSpecies - 1-d numpy Array, shadows the species argument to this 
                        function. Initially set to all zeros. When a species is 
                        solvable through unique SLS we change the appropriate 
                        element from 0 to 1 to indicate that fact.
                         
        massFragsUsed - A list that specifies, in order, which 
                        mass fragments were used for the SLS unique 
                        solutions.
    """
        
    extentSolvable = 0 
    summedSignificance = 0
    solvedSpecies = np.zeros(len(species))
    massFragsUsed = []
    
    for molecule in range(len(species) + 1):
            
        rowSum = np.sum(onesAndZeros, axis=1)
    
        #Is there a 1 in the rowSum numpy array?
        #Note: See https://stackoverflow.com/questions/432112/is-there-a-numpy-function-to-return-the-first-index-of-something-in-an-array
        #For info on np.where function
        onesIndices = np.where(rowSum == 1)[0]
        if len(onesIndices) > 0: #If there is at least 1, then there is at least one unique row.
                        
            #get a list of the unique value coordinates within onesAndZeros
            #where those coordinates are stored as tuples (row,col)
            #Below, for (rowOfUnique, colOfUnique), the rowOfUnique is the index of the unique mass fragment *within* the massFrags combination provided, and colOfUnique is the index of them molecule.
            uniqueCoordinates = []
            for rowOfUnique in onesIndices:
                colOfUnique = np.where(onesAndZeros[rowOfUnique,:] == 1)
                uniqueCoordinates.append((rowOfUnique, colOfUnique))
            
            #If there is more than one unique coordinate
            #pick the best one based on the significanceMatrix
            #Below, for (r,c), the r is the index of the unique mass fragment *within* the massFrags combination provided, and c is the index of them molecule.
            if len(uniqueCoordinates) > 1:
                coordSignificances = [significanceMatrix[r,c] for 
                    (r,c) in uniqueCoordinates]
                
                #reorder uniqueCoordinates by significance
                orderedUniqueCoordinates = [
                    x for (molecule,x) in sorted(zip(
                        coordSignificances, uniqueCoordinates),
                        key=lambda pair:pair[0])]
                
                #Take the coordinates with the highest significance
                bestUniqueRowIdx = orderedUniqueCoordinates[-1][0]
                bestUniqueColIdx = orderedUniqueCoordinates[-1][1]
            
            # if there is only a single unique value then it's the best
            else:
                bestUniqueRowIdx = uniqueCoordinates[0][0]
                bestUniqueColIdx = uniqueCoordinates[0][1]
            
            
            
            #indicate that the mass fragment corresponding to this row
            #was used. 
            massFragsUsed.append(massFrags[bestUniqueRowIdx])
            
            #indicate that this species is now solved for
            solvedSpecies[bestUniqueColIdx] = 1
            
            #Set the column that corresponds to the species we have just 
            #solved to zeros. Now we can sum rows again
            onesAndZeros[:,bestUniqueColIdx] = 0
        
            extentSolvable += 1
            summedSignificance += float(
                significanceMatrix[bestUniqueRowIdx, bestUniqueColIdx])
            #print(massFrags)
            #print(extentSolvable, solvedSpecies, summedSignificance, massFragsUsed)
                
        #Else there was no unique one in the rowSum. 
        #We can proceed no further with unique SLS.
        else:
            #print(extentSolvable, solvedSpecies, summedSignificance, massFragsUsed)
            return (extentSolvable, summedSignificance,
                    solvedSpecies, massFragsUsed)   
                      
                      
def generateSignificanceMatrix(intensityMatrix, moleculesLikelihood=None, minThreshold=None, maxIntensityPossible=100): 
    """
    Given a matrix of reference intensities, with a horizontal
    axis of species and a vertical axis of mass fragments (excluding 
    label columns/rows), this function
    determines the significance of each intensity in the matrix
    
    Arguments:
        intensityMatrix - An n by m 2-d numpy array that contains the 
                          reference intensities corresponding to 
                          n mass fragments and m species.
                          
      moleculesLikelihood - this is a one-dimensional array of values between 0 and 1
                         which are basically the probabilities of each of the m species.
                         
        minThreshold -  this is basically the lower cutoff for intensity of a fragment, and is necessary
                         to have in order to make the best possible significance matrix.
                         
     minThreshold -  this is basically the maximum intensity possible of a fragment, and is necessary
                         to have in order to make the best possible significance matrix.
                          
    Returns:
        significanceMatrix - An n by m 2-d numpy array that mirrors 
                             intensityMatrix. Each element in significanceMatrix
                             is the significance of the corresponding intensity 
                             from intensityMatrix.
    """

    numberOfMolecules = len(intensityMatrix[0])
    if moleculesLikelihood == None: #if it is not provided, assume a likelihood of 1 for each molecule.
        np.ones(numberOfMolecules)
    
    #We need to form  2-d array that shadows the intensity data
    #from intensityMatrix.
    #Initially we will set it to zeros and then populate in loop.
    significanceMatrix = np.zeros(
        intensityMatrix.shape)
        
    #loop through columns (i.e. species)
    ##TODO this loop was designed to work with the ElemSignificanceCalculator
    ## and IndElemSignificanceCalculator functions in MSRESOLVE that are
    ## apparently designed for a different application.
    ## If it proves to be a performance bottleneck (which I doubt since 
    ## we aren't looping here), we can rewrite those to cut down on the 
    ## required number of row sums. As it is now we sum a row for every element
    ## in each column, and we do this for every column.
    for (columnIdx, _) in enumerate(significanceMatrix.T):
        #determine, for every element in this column (i.e. for this
        #species) the element's significance
        colSignificances = MSRESOLVE.ElemSignificanceCalculator(
            intensityMatrix, columnIdx, moleculesLikelihood, 
            minThreshold=minThreshold, maxIntensityPossible=maxIntensityPossible)
            
        #Now record these significances
        significanceMatrix[:,columnIdx] = colSignificances
        
    return significanceMatrix
    
if __name__ == "__main__":
    massFrags = np.array([15,25,26])
    species = ['Ethylene (Ethene)', 'Ethanol']
    onesAndZeros = np.array([
        [0,1],
        [1,0],
        [1,1]
    ])
    significanceMatrix = np.array([
        [0,3],
        [7,0],
        [600,4]
    ])
    
    
    
    # print(ExtentOfSLSUniqueSolvable(massFrags, species, 
    #     onesAndZeros, significanceMatrix))
        
    #obj = [(5,1), (4,5) ,(2,3), (2,2)]
    obj = [(2,2), (2,3), (4,5), (5,1)]
    
    #parallel = ['51', '45', '23', '22']
    parallel = ['22', '23', '45', '51']

    new_obj = (5,2)
    new_p = '52'
    
    max_length = 4
    
    #print(bisect.bisect(obj, new_obj))
    
    print(MSRESOLVE.storeAndPop(
        obj, new_obj, parallel, new_p, max_length,
        excludeDuplicates=True, optimumType="Maximum" )
        )
    
    