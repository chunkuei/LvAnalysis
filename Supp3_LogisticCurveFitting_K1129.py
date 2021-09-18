#############################################################################################################
##### This program is used for negative simgmoid curve fitting. It reads 2-columns text data as (x, y)  #####
#############################################################################################################
import tkinter as tk
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import statistics as stat
import math
import tkinter.filedialog
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
##### Select the 2-column text file for curve fitting ########################################################
application_window = tk.Tk()
my_filetypes = [('all files', '.*'), ('text files', '.txt')]
filename = tkinter.filedialog.askopenfilename(parent = application_window,
                                              initialdir = 'f:/mydatafiles/SpikingTimes/',
                                              title = 'select a txt files with 2-column data',
                                              filetypes = my_filetypes)
application_window.destroy()
###############################################################################################################
##### Assign x-value and y-value from the data in the text file ###############################################
a = np.loadtxt(filename,dtype='float', delimiter='\t')
x, y = a[:,0], np.log10(a[:,1])
x_now = np.arange(np.min(x),np.max(x), 0.01, dtype=float)  # x_now values used for simulated curves
###############################################################################################################
###### Inital values and bounds of parametrical values for curve fitting ######################################
ao, aoLo, aoHi = 100, 0., np.inf
b, bLo, bHi = 200, -np.inf, np.inf
c, cLo, cHi = 0.2649, 0., np.inf
d, dLo, dHi = -0.7165, -np.inf, 0
###############################################################################################################
##### Define logistic function as y = {ao / [1 + b * c^(-x)]} + d    ##########################################
def Logistic(x, ao, b, c, d):
    return ao/(1+b*np.power(c,-x)) + d
###############################################################################################################
##### Cruve fitting to find optimal parametrical values #######################################################
pars, cov = curve_fit(f = Logistic, xdata = x, ydata = y, p0 = [ao,b,c,d],
                      bounds = ([aoLo,bLo,cLo, dLo],[aoHi,bHi,cHi,dHi]))
###############################################################################################################
##### Simulated curves and statistics #########################################################################
y_expc = Logistic(x, *pars)
y_now = Logistic(x_now, *pars)
stdevs = np.sqrt(np.diag(cov))
chisq, ptail = chisquare(y, y_expc, ddof = len(pars))

CovOE = np.cov(y, y_expc) # Covariance of observed y-values and expected y-values
CorrOE = CovOE[0,1]/np.sqrt(CovOE[0,0]*CovOE[1,1]) # coefficient of correlation
Rsquare = np.round((CorrOE**2), 5) # coefficient of determination
##### Report and plot #########################################################################################
print(filename)
print("chi-square = ", chisq)
print("p-tail = ", ptail)
print("a,b,c,d: ", pars)
print("stdevs: ", stdevs)
print("Coefficient of Determiantion: ", Rsquare)

fig = plt.figure(1)
plt.plot(x, y, 'ko', label='data', markersize = 6)
plt.plot(x_now, y_now, 'r-', linewidth = 2)
plt.show()

