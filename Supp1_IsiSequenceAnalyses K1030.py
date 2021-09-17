#### This program is designed for analyses of spiking time series data, taking in one data file ######
#         Acquireing all the metrics describing the rates and the variation of spiking activity      #
#         Producing figures of IFR, dLv, Cwt, Wes, and mPDE                                          #
######################################################################################################

import tkinter as tk
import numpy as np
import pywt
import scipy.signal
from scipy.signal import find_peaks
from scipy.stats import f
from scipy.stats import t
from scipy.stats import gmean
from scipy.stats import hmean
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import statistics as stat
import tkinter.filedialog
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
############ Enter your choices here....############################
wave = 'morl' #'mexh' or 'morl' or 'gaus8'

####################################################################


############ Read a data file for analyses #########################
application_window = tk.Tk()
my_filetypes = [('all files', '.*'), ('text files', '.txt')]
filename = tkinter.filedialog.askopenfilename(parent = application_window,
                                              initialdir = 'd:/mydatafiles/SpikingTimes/OligoTxtFiles/',
                                              title = 'select multiple txt files',
                                              filetypes = my_filetypes)
application_window.destroy()

#####################################################################


############## Statistics for control ###############################
a = np.loadtxt(filename,dtype='float', delimiter=' ')
ISI = a[1:]-a[0:-1]
a = a - a.min()
whr = np.where(a > 1200)
ISI = ISI[0:(whr[0][0]-1)]

LvS, LvT = 0, 0  # Acquireing static Lv
for j in range(len(ISI)-1):
    LvT = ((ISI[j+1]-ISI[j])/(ISI[j+1]+ISI[j]))**2
    LvS = LvT + LvS
Lv = round(LvS *3/(len(ISI) - 1), 5)   

mean, g_mean, h_mean, std, cv = round(np.mean(ISI), 5), round(gmean(ISI), 5), round(hmean(ISI),5), round(np.std(ISI),5), round(np.std(ISI)/np.mean(ISI),5)
median, mode = round(np.median(ISI),5), round(stat.mode(ISI),5)
mIFR_Cont = np.round(np.mean(1/ISI),5)
print("mean = ", mean)
print("std = ", std)
print("Cv = ", cv)
print("Lv = ", Lv)
print("median = ", median)
print("mode = ", mode)
print("mIFR_Cont = ", mIFR_Cont)
print()
mIsiCont = mean

############## Parametrical values to be recorded ######################
MeanCont = str(mean)
CvCont = str(cv)
LvCont = str(Lv)
mIFR_Cont = str(mIFR_Cont)
########################################################################

############### Statistics for treatments ##############################
a = np.loadtxt(filename,dtype='float', delimiter=' ')
ISI = a[1:]-a[0:-1]
nT = len(ISI)
a = a - a.min()
whr = np.where(a > 1800)
ISI = ISI[whr[0][0]:nT]

LvS, LvT = 0, 0 # Acquireing static Lv
for j in range(len(ISI)-1):
    LvT = ((ISI[j+1]-ISI[j])/(ISI[j+1]+ISI[j]))**2
    LvS = LvT + LvS
Lv = round(LvS *3/(len(ISI) - 1), 5)
    
mean, g_mean, h_mean, std, cv = round(np.mean(ISI), 5), round(gmean(ISI), 5), round(hmean(ISI),5),round(np.std(ISI),5), round(np.std(ISI)/np.mean(ISI),5)
median, mode = round(np.median(ISI),5), round(stat.mode(ISI),5)
mIFR_Treat = np.round(np.mean(1/ISI),5)
print("mean = ", mean)
print("std = ", std)
print("Cv = ", cv)
print("Lv = ", Lv)
print("median = ", median)
print("mode = ", mode)
print("mIFR_Treat = ", mIFR_Treat)
print(filename)

############### Parametrical values to be recorded #######################
MeanTret = str(mean)
CvTret = str(cv)
LvTret = str(Lv)
mIFR_Treat = str(mIFR_Treat)
############################################################################

############### Time course of firing and Lv ###############################
TimeSeg = 32  # User control of Time used to calculate Lv time course
Nisi = int(round((TimeSeg / mIsiCont) / 2, 0) * 2) # get an even number of ISI to calculate LV in 'TimeSeg' seconds
if Nisi < 4:
    Nisi = 4
print("Time per episode in s for Lv calcualtion = ", Nisi * mIsiCont)
print("Nisi = ", Nisi)
a = np.loadtxt(filename,dtype='float', delimiter=' ')
ISI = a[1:]-a[0:-1]
nT = len(ISI)
a = a - a.min()
x = a[1:]
y = 1/ISI

fig = plt.figure(1)
ax = fig.add_subplot(111)

plt.plot(x/60, y, 'k-', linewidth = 1.5)
#plt.title('Instantaneous firing')
plt.xlabel('Time course', fontsize = 24)
plt.yscale('log')
plt.ylabel('IFR (Hz)', fontsize = 24)
#plt.axis(linewidth = 5)

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(3)
plt.xticks(fontsize =18)
plt.yticks(fontsize = 18)
ax.tick_params(which='both', width=3)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_minor_locator(MultipleLocator(5))
plt.show()

############## Acquire dynamic Lv from Nisi consecutive ISI, 50% overlap #####
ShiftIndx = int(Nisi/2)
LvX = a[ShiftIndx:-(ShiftIndx+1)] #Time tags of dynamic Lv in 's'
LvY = []
for i in range(len(ISI)-Nisi):
    LvS, LvT = 0, 0
    for j in range(Nisi):
        LvT = ((ISI[i+j+1]-ISI[i+j])/(ISI[i+j+1]+ISI[i+j]))**2
        LvS = LvT + LvS
    Lv = round(LvS *3/(Nisi-1), 5)
    LvY.append(Lv) # dynamic Lv series
    
fig = plt.figure(2)
ax = fig.add_subplot(111)

plt.plot(LvX/60, LvY, 'k-', linewidth = 1.5)
plt.xlabel('Time course', fontsize = 24)
plt.ylabel('LV', fontsize =24)

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(3)
plt.xticks(fontsize =18)
plt.yticks(fontsize = 18)
ax.tick_params(which='both', width=3)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_minor_locator(MultipleLocator(5))
plt.show()

############# Calculate dLv mean and CV ###########################################

whrXcont = np.where(LvX > 1200) # get the array index where the time tags > 1200-s
whrXtret = np.where(LvX >= 1800) # get the array index where the time tags >= 1800-s
LvYcont = LvY[0:(whrXcont[0][0]-1)]
LvYtret = LvY[whrXtret[0][0]:len(LvY)]
dLvNcont = len(LvYcont)
dLvNtret = len(LvYtret)
dLvMeanCont, dLvStdCont, dLvCvCont = round(stat.mean(LvYcont),4), round(stat.stdev(LvYcont),4), round(stat.stdev(LvYcont)/stat.mean(LvYcont),4)
dLvMeanTret, dLvStdTret, dLvCvTret = round(stat.mean(LvYtret),4), round(stat.stdev(LvYtret),4), round(stat.stdev(LvYtret)/stat.mean(LvYtret),4)
t_dLv = abs(dLvMeanCont-dLvMeanTret)/np.sqrt(dLvStdCont**2/dLvNcont + dLvStdTret**2/dLvNtret)
rFactor = (dLvStdCont**2/dLvNcont) / (dLvStdCont**2/dLvNcont + dLvStdTret**2/dLvNtret)
df_t = 1 / ((rFactor**2/(dLvNcont-1) + (1- rFactor**2)/(dLvNtret-1)))
Ptail_dLv = np.round((1 -t.cdf(t_dLv, df_t)),5)
            
print("dLvMeanCont, dLvStdCont, dLvCvCont = ", dLvMeanCont, dLvStdCont, dLvCvCont)
print("dLvMeanTret, dLvStdTret, dLvCvTret = ", dLvMeanTret, dLvStdTret, dLvCvTret)
print("Ptail_dLv = ", Ptail_dLv)

############### CWT section ##########################################################
LvReplot = "y"
Period = 600. # in s for the longest scale calculation
reInterval =0.05 # resampling interval = sampling period
N_coef = 400. # number of coef used to describe frequency range

n = int((max(LvX)-min(LvX))/reInterval)
xre = np.linspace(min(LvX),max(LvX),n, endpoint = False)
yre = np.interp(xre,LvX,LvY)

if LvReplot == "y":
    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    plt.plot(LvX/60, LvY, 'ko')
    plt.plot(xre/60,yre,'bo', markersize = 2)
    #plt.title('LV')
    plt.xlabel('Time course', fontsize = 24)
    plt.ylabel('LV', fontsize =24)
    #plt.axis(linewidth = 5)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(3)
    plt.xticks(fontsize =18)
    plt.yticks(fontsize = 18)
    ax.tick_params(which='both', width=3)
    ax.tick_params(which='major', length=7)
    ax.tick_params(which='minor', length=4, color='k')
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    plt.show()


#### Look into the freq range of 0.004-0.04 Hz (i.e., 250-25 s/ period)
scale = np.arange(0.,N_coef, 1.) #by the following equation, f = 0.1-0.001 Hz
for i in np.arange(int(N_coef)):
    scale[i] = 32.5 + (np.exp(np.log(Period/reInterval)/N_coef *i)-1)

f = pywt.scale2frequency(wave, scale)/reInterval #'mexh' or 'gaus1' or 'morl'
ScaleUp = max([i for i,x in enumerate(f) if x>0.001667])# set low freq band
ScaleLow = max([i for i,x in enumerate(f) if x>0.05])# set high freq band
ScaleNum = round(((ScaleUp -ScaleLow)/20 +0.5),0) * 20
scale = scale[int(ScaleUp-ScaleNum+1):int(ScaleUp+1)]
coef, freqs=pywt.cwt(yre,scale,wave,sampling_period = reInterval, method = 'fft')

############### Calculate Cwt Coef energy ##############################################
CoefSqSum = np.sum(coef**2,axis = 0)
energy = coef**2 / CoefSqSum  # acquire probability density of Coef energy

############### Plot cwt coefficient ###################################################
fig, ax = plt.subplots()
cRatio = 0.75 # Enter a numer <1 to increase contrast but lose resolution
cLo = np.min(coef) * cRatio
cHi = np.max(coef) * cRatio
ax.imshow(coef, cmap = 'viridis', vmax = cHi, vmin = cLo)
ax.set_aspect(50)
plt.xlabel("Time course (min)", fontsize = 24)
plt.ylabel("Frequency (Hz)", fontsize = 24)


locs, labels = plt.xticks()
x_Big = int(round(n*reInterval /300)+1)
locsBig = np.arange(x_Big)
for i in np.arange(x_Big):
    locsBig[i] = i*(300/reInterval)
    if i < len(locs):
        labels[i] = plt.text(locsBig[i],-100000,str(5*i))
    else:
        temp = plt.text(locsBig[i],-100000,str(5*i))
        labels.append(temp)
    
    
ylocs, ylabels = plt.yticks()

cc = 1
for i in ylocs[1:]:
    if i == 0:
        f_now = round(freqs[int(i)], ndigits = 3)
    else:
        f_now = round(freqs[int(i)-1], ndigits = 3)
    ylabels[cc] = plt.text(-50000,ylocs[cc],str(f_now))
    cc = cc+1
   
plt.tick_params(axis='x', left='off', top='off',
                right='off', bottom='off', labelleft='off',
                labeltop='off', labelright='off', labelbottom='off')
plt.xticks(locsBig,labels,fontsize = 18)
plt.yticks(ylocs[1:],ylabels[1:],fontsize =18)
ax.xaxis.set_minor_locator(MultipleLocator(60/reInterval))

ax.tick_params(which='both', width=2)
ax.tick_params(which='major', length=6)
ax.tick_params(which='minor', length=3, color='k')

[row,col] = np.shape(coef)
plt.xlim(left = 0., right = col)

plt.show()

###########################################################################################

############## Plot Coef indicator #############################################################
indicator = "y"  #User input 'y' to plot the cmap indicator
if indicator == "y":
    coef_indicator = np.ndarray((len(freqs),200),dtype = 'float')
    cDiv = (cHi-cLo)/len(freqs)
    for i in np.arange(0,len(freqs)):
        coef_indicator[i,:] = cHi - cDiv * i

    fig, ax = plt.subplots()
    ax.imshow(coef_indicator, cmap = 'viridis')
    ax.set_aspect(50)
    ylocs, ylabels = plt.yticks()
    cc = 0
    for i in ylocs:
        ylabels[cc] = plt.text(500,i,str(round(cHi-cDiv*i, 2)))
        cc = cc+1
    plt.yticks(ylocs[1:],ylabels[1:],fontsize = 14)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=True,# ticks along the top edge are off
    labelleft=False, # ticks along the top edge are off
    labelright=False, width = 2) # labels along the bottom edge are off
    plt.show()
    
############## Plot cwt energy ##############################################################
fig, ax = plt.subplots()
cRatio = 0.75 # Enter a numer <1 to increase contrast but lose resolution
cLo = np.min(energy) * cRatio
cHi = np.max(energy) * cRatio
ax.imshow(energy, cmap = 'gnuplot2', vmax = cHi, vmin = cLo)
ax.set_aspect(50)
plt.xlabel("Time course (min)", fontsize = 24)
plt.ylabel("Frequency (Hz)", fontsize = 24)

locs, labels = plt.xticks()
x_Big = int(round(n*reInterval /300)+1)
locsBig = np.arange(x_Big)
for i in np.arange(x_Big):
    locsBig[i] = i*(300/reInterval)
    if i < len(locs):
        labels[i] = plt.text(locsBig[i],-100000,str(5*i))
    else:
        temp = plt.text(locsBig[i],-100000,str(5*i))
        labels.append(temp)
    
    
ylocs, ylabels = plt.yticks()

cc = 1
for i in ylocs[1:]:
    if i == 0:
        f_now = round(freqs[int(i)], ndigits = 3)
    else:
        f_now = round(freqs[int(i)-1], ndigits = 3)
    ylabels[cc] = plt.text(-50000,ylocs[cc],str(f_now))
    cc = cc+1
   
plt.tick_params(axis='x', left='off', top='off',
                right='off', bottom='off', labelleft='off',
                labeltop='off', labelright='off', labelbottom='off')
plt.xticks(locsBig,labels,fontsize = 18)
plt.yticks(ylocs[1:],ylabels[1:],fontsize = 18)
ax.xaxis.set_minor_locator(MultipleLocator(60/reInterval))

ax.tick_params(which='both', width=2)
ax.tick_params(which='major', length=6)
ax.tick_params(which='minor', length=3, color='k')

[row,col] = np.shape(energy)
plt.xlim(left =0.,right = col)
plt.show()

###############################################################################################

############## Plot energy indicator ##########################################################
indicator = "y"  #User input 'y' to plot the cmap indicator
if indicator == "y":
    energy_indicator = np.ndarray((len(freqs),200),dtype = 'float')
    cDiv = (cHi-cLo)/len(freqs)
    for i in np.arange(0,len(freqs)):
        energy_indicator[i,:] = cHi - cDiv * i

    fig, ax = plt.subplots()
    ax.imshow(energy_indicator, cmap = 'gnuplot2')
    ax.set_aspect(50)
    ylocs, ylabels = plt.yticks()
    cc = 0
    for i in ylocs:
        ylabels[cc] = plt.text(500,i,str(round(cHi-cDiv*i, 4)))
        cc = cc+1
    plt.yticks(ylocs[1:],ylabels[1:],fontsize = 14)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=True,# ticks along the top edge are off
    labelleft=False, # ticks along the top edge are off
    labelright=False, width = 2) # labels along the bottom edge are off
    plt.show()


############## Calcualte means of probabiltiy density of Coef energy under control and during drug application #####################
start_time = Period / 2
start_time = int(start_time / reInterval)
duration = 1200. - Period
duration = int(duration / reInterval)
minus_end = n-(start_time + duration)

briefC = (energy[:,start_time:-minus_end])
briefT = (energy[:,int((1200+Period/2)/reInterval):-300]) # end at 15-s earlier
mean_energyC = np.mean(briefC, axis = 1)
mean_energyT = np.mean(briefT,axis = 1)
mean_energyD = mean_energyT-mean_energyC
std_energy = np.std(energy[:,start_time:-minus_end], axis = 1)

s = stat.stdev(mean_energyD)
thre = 2.57583 * s / (len(mean_energyD)**0.5)

peaks, p_ = find_peaks(mean_energyD, height = thre, distance = 10, width =5)
valleys, v_ = find_peaks(-mean_energyD, height = thre, distance = 10, width =5)

PeaksAF =np.transpose(np.append([mean_energyD[peaks]],[freqs[peaks]],axis=0))
ValleysAF =np.transpose(np.append([mean_energyD[valleys]],[freqs[valleys]],axis=0))
vpAF = np.append(ValleysAF,PeaksAF)

print('peaks at Hz = ', np.round(freqs[peaks], decimals =5))
print('peak coefficients = ', np.round(mean_energyD[peaks], decimals = 5))
print('valleys at Hz = ', np.round(freqs[valleys], decimals =5))
print('valleys coefficients = ', np.round(mean_energyD[valleys], decimals = 5))
print(filename)

fig = plt.figure(22)
ax = fig.add_subplot(111)

plt.plot(freqs,mean_energyC,'ko', markersize = 4)
plt.plot(freqs,mean_energyT,'r-', linewidth = 2)
plt.plot(freqs,mean_energyD, 'b:')
plt.xscale("log") # log or linear
plt.xlabel("Frequency (Hz)", fontsize = 24)
plt.ylabel("energy mean", fontsize = 24)

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(3)
plt.xticks(fontsize =18)
plt.yticks(fontsize = 18)
ax.tick_params(which='both', width=3)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')

plt.show()

CovCT = np.cov(mean_energyC, mean_energyT)
CorrCT = CovCT[0,1]/np.sqrt(CovCT[0,0]*CovCT[1,1]) 
CoD = np.round((1 - CorrCT**2), 5)
print("Coefficient of Deviation caused by the treatment: ", CoD)
############### Parametrical values to be recorded ################################################

CoD = str(CoD)

############ Append record to a file ##########################

MyRecord = [filename, MeanCont, mIFR_Cont, CvCont, LvCont,
             MeanTret, mIFR_Treat, CvTret, LvTret, CoD,
             dLvMeanCont, dLvCvCont, dLvNcont,
             dLvMeanTret, dLvCvTret, dLvNtret, Ptail_dLv]

# Open the file in append & read mode ('a+')
RecordFilename = "d:/MyDataFiles/SpikingTimes/sLv_dLv_Records.txt" # Enter the file path and filename
def MyAppendList(file_name, list_to_append):
    """Append given list as a new line at the end of file"""
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        # Move read cursor to the start of file.
        file_object.seek(0)
        # If file is not empty then append '\n'
        data = file_object.read(100)
        if len(data) > 0:
            file_object.write("\n")
        # Append list at the end of file
        for listItem in list_to_append:
            file_object.write('%s\t' % listItem)



MyAppendList(RecordFilename,MyRecord)


