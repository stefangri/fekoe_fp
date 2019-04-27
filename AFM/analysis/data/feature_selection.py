import matplotlib.pyplot as plt
import numpy as np
def select_features(fig, ax):

        xpeak = []  # x and
        ypeak = []  # y arrays for peak coordinates
        global line
        line, = ax.plot(xpeak, ypeak, 'ro', markersize = 10)


        def onclickpeaks(event):

            if event.button == 1: # left mouse click to add data point
                xpeak.append(event.xdata)               # append x data and
                ypeak.append(event.ydata)               # append y data
                line.set_xdata(xpeak)
                line.set_ydata(ypeak)
                plt.draw()                       # and show it


            if event.button == 3: #right mouse click to remove data point
                if xpeak != []:
                    xdata_nearest_index = ( np.abs(xpeak - event.xdata) ).argmin() #nearest neighbour
                    del xpeak[xdata_nearest_index] #delte x
                    del ypeak[xdata_nearest_index] #and y data point
                    line.set_xdata(xpeak) #update x
                    line.set_ydata(ypeak) #and y data in plot
                    plt.draw()   # and show it


        # actual execution of the defined function oneclickpeaks
        cid = fig.canvas.mpl_connect('button_press_event', onclickpeaks)
        figManager = plt.get_current_fig_manager()  # get current figure
        figManager.window.showMaximized()           # show it maximized

        return xpeak, ypeak


fig, ax = plt.subplots()
x1, y1, x2, y2 = np.genfromtxt('quadrate_profile_vertikal', unpack = True)
ax.plot(x1, y1)
select_features(fig, ax)
