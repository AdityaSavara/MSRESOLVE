INTRODUCTION TO THE MSRESOLVE PROGRAM: 
A program for converting mass spectrometry signals to concentrations. Baseline corrections, smoothing, solving overlapping patterns, mass spectrometer tuning, and more. 

The Need for MSRESOLVE and what it offers
During “solving” of collected mass spectrometry spectra to extract concentrations, there are several sources of challenges.
A)	For time dependent signal collection, pre-processing may be required (such as baseline corrections, smoothing, etc.)
B)	There may be overlapping signals, making solution difficult.
C)	In many cases a particular molecule’s calibration may not be possible or practical at the  instrument where the collection is being done. This creates two challenges.  Firstly, that one may need to rely upon an externally collected reference pattern that may not exactly match how the molecule fragments in one’s own instrument, and Secondly that one requires a method for converting the signals for that molecule into (approximate) concentrations even in the absence of a calibration.
MSRESOLVE addresses each of the above challenges, as follows:
A)	MSRESOLVE has various pre-processing functions, including baseline correction, and smoothing.
B)	MSRESOLVE has several methods for resolving the signals: Sequential Linear Subtraction (SLS), the Matrix Inverse Method, and also brute force grid (regression).  
a.	For the SLS and the inverse methods, it is possible to extract error bars for the final solved concentrations.
C)	MSRESOLVE is able to convert mass spectrometry signals into a common concentration scale even for uncalibrated signals, provided that reference concentration patterns are available. Furthermore, MSRESOLVE can also correct for differences in mass spectrometers tuning.
Suggested Procedure for Solving Time Series Data of Mass Signals

1)	Create the files for the reference spectra and collected data – don’t delete any signals.
2)	Create an MSRESOLVE run with SLS Unique and see what masse MSRESOLVE chooses, with ExportedSLSUniqueMassesUsedInSolvingMolecules.
3)	If unsatisfied, start narrowing things down with chosen masses.
4)	Also start using some reference file threshold filtering . 
a.	UserChoices['minimalReferenceValue']['referenceValueThreshold'] = [1.0]  #this is what I am suggesting that you use.
5)	If the application warrants doing so, include more sophisticated features of MSRESOLVE, such as mass spectrum tuning correction.


There is a manual in the documentation directory, and also there is an example analysis under ExampleAnalysis.  