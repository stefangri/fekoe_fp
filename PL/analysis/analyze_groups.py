from spectrum import *

# possible peaks (breit_wigner == fano)
# implemented are: breit_wigner, lorentzian, gaussian, voigt
peaks = ['lorentzian']

# select folder you want to analyze and initialize everything
# it doesn't matter if there is one or more files in the folder
spec = spectrum('smallmap')

# Select the interesting region in the spectrum,
# by clicking on the plot
spec.SelectSpectrum()

# find all Muons and remove them
spec.DetectAllMuons()
spec.RemoveAllMuons()

spec.GroupSpectra()

# Function opens a window with the data,
# you can select the regions that do not belong to
# the third degree polynominal background signal
# by clicking in the plot
spec.SelectBaseline()

# fit the baselines
spec.FitAllBaselines()

# Function that opens a Window with the data,
# you can choose initial values for the peaks by clicking on the plot.
# You have to choose peaks for all spectra to get the proper starting
# values. -> Improvement needed
spec.SelectGroupedPeaks(peaks)

# Fit all spectra with initial values provided by SelectBaseline()
# and SelectAllPeaks()
# if you need to see the fit results set report to True,
# otherwise set it to false
spec.FitAllGroupedSpectra(peaks, report=False)
#spec.FitGroups(peaks, groups = [1, 2])

# Save the results of the fit in txt-files
spec.SaveAllFitParams(peaks)
