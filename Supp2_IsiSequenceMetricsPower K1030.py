##### This program is designed for analyses of spiking time series data, taking in multiple data files   #####
#          Acquiring all the metrics describing varaition of ISI sequence, rhythms of dLv, and power         #
##############################################################################################################

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
filenames = tkinter.filedialog.askopenfilenames(parent = application_window,
                                              initialdir = 'd:/mydatafiles/SpikingTimes/OligoTxtFiles/',
                                              title = 'select multiple txt files',
                                              filetypes = my_filetypes)
application_window.destroy()
#####################################################################

filenames = np.sort(filenames)

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

def MyAppendList2(file_name, list_to_append):
    """Append given list as a new line at the end of file"""
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        # Move read cursor to the start of file.
        #file_object.seek(0)
        # If file is not empty then append '\n'
        #data = file_object.read(100)
        #if len(data) == 0:
        #    file_object.write("\n")
        # Append list at the end of file
        for listItemRC in list_to_append:
            for listItem in listItemRC:
                file_object.write('%s\t' % listItem)
            file_object.write("\n")

for filename in filenames:
    ############## Statistics for control ###############################
    a = np.loadtxt(filename,dtype='float', delimiter=' ')
    ISI = a[1:]-a[0:-1]
    a = a - a.min()
    whr = np.where(a > 1200)
    ISI = ISI[0:(whr[0][0]-1)]
    LvS, LvT = 0, 0
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
    LvS, LvT = 0, 0
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
    #############################################################################
    Nisi = int(round((TimeSeg / mIsiCont) / 2, 0) * 2) # get an even number of ISI to calculate LV in 'TimeSeg' seconds
    if Nisi < 4:
        Nisi = 4
    print("Time per episode in s for Lv calcualtion = ", Nisi * mIsiCont)
    print("Nisi = ", Nisi)
    #user inputs an even number as the number of ISI segments to calcualte Lv
    a = np.loadtxt(filename,dtype='float', delimiter=' ')
    ISI = a[1:]-a[0:-1]
    nT = len(ISI)
    a = a - a.min()
    x = a[1:]
    y = 1/ISI
    # count Lv by Nisi consecutive ISI, 50% overlap
    ShiftIndx = int(Nisi/2)
    LvX = a[ShiftIndx:-(ShiftIndx+1)]
    LvY = []
    for i in range(len(ISI)-Nisi):
        LvS, LvT = 0, 0
        for j in range(Nisi):
            LvT = ((ISI[i+j+1]-ISI[i+j])/(ISI[i+j+1]+ISI[i+j]))**2
            LvS = LvT + LvS
        Lv = round(LvS *3/(Nisi-1), 5)
        LvY.append(Lv)
    ############################ Calculate dLv mean and CV ########################
    whrXcont = np.where(LvX > 1200)
    whrXtret = np.where(LvX >= 1800)
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
    ############### CWT section ######################################################
    Period = 600. # in s for the longest scale calculation
    reInterval =0.05 # resampling interval = sampling period
    N_coef = 400. # number of coef used to describe frequency range
    n = int((max(LvX)-min(LvX))/reInterval)
    xre = np.linspace(min(LvX),max(LvX),n, endpoint = False)
    yre = np.interp(xre,LvX,LvY)
    #look into the freq range of 0.004-0.04 Hz (i.e., 250-25 s/ period)
    scale = np.arange(0.,N_coef, 1.) #by the following equation, f = 0.1-0.001 Hz
    for i in np.arange(int(N_coef)):
        scale[i] = 32.5 + (np.exp(np.log(Period/reInterval)/N_coef *i)-1)
    f = pywt.scale2frequency(wave, scale)/reInterval #'mexh' or 'gaus1' or 'morl'
    ScaleUp = max([i for i,x in enumerate(f) if x>0.001667])# set low freq band
    ScaleLow = max([i for i,x in enumerate(f) if x>0.05])# set high freq band
    ScaleNum = round(((ScaleUp -ScaleLow)/20 +0.5),0) * 20
    scale = scale[int(ScaleUp-ScaleNum+1):int(ScaleUp+1)]
    coef, freqs=pywt.cwt(yre,scale,wave,sampling_period = reInterval, method = 'fft')
    ################ Calculate Coef energy #############################################
    CoefSqSum = np.sum(coef**2,axis = 0)
    energy = coef**2 / CoefSqSum  #energy probability density
    ################ Calcualte mean energy under control and during drug applied #######
    start_time = Period / 2
    start_time = int(start_time / reInterval)
    duration = 1200. - Period
    duration = int(duration / reInterval)
    minus_end = n-(start_time + duration)
    briefC = (energy[:,start_time:-minus_end])
    briefT = (energy[:,int((1200+Period/2)/reInterval):-300]) # ends at 15-s earlier
    mean_energyC = np.mean(briefC, axis = 1)*freqs*1000
    mean_energyT = np.mean(briefT,axis = 1)*freqs*1000
    mean_energyD = mean_energyT-mean_energyC
    std_energy = np.std(energy[:,start_time:-minus_end], axis = 1)
    s = stat.stdev(mean_energyD)
    work = str(round(np.sum(mean_energyC),5))
    boost = str(round(np.sum(mean_energyD)/np.sum(mean_energyC),5))
    workCf = round(freqs[np.argmax(mean_energyC)]*1000,5)
    workCa = round(np.max(mean_energyC),5)
    workTf = round(freqs[np.argmax(mean_energyT)]*1000,5)
    workTa = round(np.max(mean_energyT),5)
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
    CovCT = np.cov(mean_energyC, mean_energyT)
    CorrCT = CovCT[0,1]/np.sqrt(CovCT[0,0]*CovCT[1,1]) 
    CoD = np.round((1 - CorrCT**2), 5)
    print("Coefficient of Deviation caused by the treatment: ", CoD)
    print("Work = ", work)
    print("Boost = ", boost)
    print("Work maximum Freq and Amp for control & Treat: ", workCf, workCa, "\n",
      workTf,workTa)
    print()
    ############# Parametrical values to be recorded ########################################
    CoD = str(CoD)
    workCf = str(workCf)
    workCa = str(workCa)
    workTf = str(workTf)
    workTa = str(workTa)
    ############## Append records to a file ##################################################
    MyRecord2 = [filename, MeanCont, mIFR_Cont, CvCont, LvCont, MeanTret, mIFR_Treat,
                CvTret, LvTret, CoD, work, boost, workCf,workCa,workTf,workTa]
    R2Len = len(MyRecord2)
    MyRecord3 = [filename, MeanCont, mIFR_Cont, CvCont, LvCont,
                 MeanTret, mIFR_Treat, CvTret, LvTret, CoD, work, boost,
                 dLvMeanCont, dLvCvCont, dLvNcont, dLvMeanTret, dLvCvTret,
                 dLvNtret, Ptail_dLv, workCf,workCa,workTf,workTa]
    rr,cc = np.shape(ValleysAF)
    MyRecordV = np.ndarray((rr,(R2Len+2)), dtype = 'U200') #No. of recording parameters + 2, for Amp & Freq
    for i in np.arange(0,rr):
        [cc0,cc1]=np.round(ValleysAF[i],6)
        MyRecordV[i] = np.append(MyRecord2,[str(cc0),str(cc1)])
    rr,cc = np.shape(PeaksAF)
    MyRecordP = np.ndarray((rr,(R2Len+2)), dtype = 'U100')
    for i in np.arange(0,rr):
        [cc0,cc1]=np.round(PeaksAF[i],6)
        MyRecordP[i] = np.append(MyRecord2,[str(cc0),str(cc1)])
    MyRecordVP = np.append(MyRecordV,MyRecordP,axis=0)
    RecordFilename2 = "d:/MyDataFiles/SpikingTimes/sLvBoostVP.txt" # Enter filepath and filename here 
    RecordFilename3 = "d:/MyDataFiles/SpikingTimes/sLv_dLv_Boost.txt" # Enter filepath and filename here 
    MyAppendList2(RecordFilename2,MyRecordVP)
    MyAppendList(RecordFilename3,MyRecord3)


