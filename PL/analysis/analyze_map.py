from mapping import *

# colors for clustering (add more if you use more than 10 clusters)
colorlist = ['red', 'orange', 'yellow',
             'blue', 'lightblue', 'turquoise',
             'green', 'lightgreen', 'white', 'purple']

# xdim:     the number of Spectra in x direction
# ydim:     the number of Spectra in y direction
# stepsize: the interval at which the mapping was collected in µm
map = mapping('smallmap', 2, 4, 10)

# plot mapping in different ways
# use PlotMapping(xmin, xmax) to integrate over fitted or raw data
# xmin:     the lowest wavenumber to be used in the mapping
# xmax:     the highest wavenumber to be used in the mapping
#map.PlotMapping(1550,1650)

# use PlotMapping(maptype='peak_fit_value_file') to map the fitted parameter
#map.PlotMapping(maptype='lorentzian_p1_sigma')

# or use PlotAllMappings for all fit parameters to be mapped
map.PlotAllMappings()

# use PlotMapping(top='file1', bot='file2') to plot a mapping of
# top/bot
map.PlotMapping(top='lorentzian_p1_height.dat',
                bot='breit_wigner_p1_height.dat')

# use PlotMapping(top='file1', bot='file2', distance=True) to plot a mapping of
# top - bot
map.PlotMapping(top='breit_wigner_p1_center.dat',
                bot='lorentzian_p1_center.dat',
                distance=True)

# Use cluster algorithms to identify something in the mapping
map.PlotClusteredPCAMapping(colorlist=colorlist)
