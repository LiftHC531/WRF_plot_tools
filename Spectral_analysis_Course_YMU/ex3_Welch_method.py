#-*- coding: utf-8 -*-
import sys
sys.path.append("..")
from package import wrf_tools #local function
from time import time as cpu_time
import os
import copy
import numpy as np
import scipy.signal
import glob
import Ngl, Nio
"""
Follow the Simulation_FFTspectrum.m
https://youtu.be/Kdz13vP0f-g?list=PLx_IWc-RN82uKOdafF4v4U5R_u4qmYaiu&t=1943
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html
https://github.com/scipy/scipy/issues/8045#issuecomment-337319294
!!detrend = False 
"""
def awgn(sig_volts, snr_db):
    """
    Additive White Gaussian Noise (AWGN)
    https://stackoverflow.com/questions/14058340/adding-noise-to-a-signal-in-python
    """
    sig_watts = sig_volts ** 2 # Voltage to Power
    # Calculate signal power and convert to dB 
    sig_watts = np.mean(sig_watts)
    sig_db = 10. * np.log10(sig_watts)
    noise_db = sig_db - snr_db
    noise_watts = 10.**(noise_db/10.)
    # Generate an sample of white noise
    mean_noise = 0.0
    noise_volts = np.sqrt(noise_watts)
    noise_volts = np.random.normal(mean_noise,noise_volts,len(sig_volts))
    #print(noise_volts.shape)
    return sig_volts + noise_volts 

def nextpow2(p):
    n = 2
    while p > n:
       n *= 2
    return np.log2(n)

def wks_setting(res,image_name):
    wks_type = "png"
    print("\033[44m---Plot: {}.{} ---\033[0m".format(image_name,wks_type))
    rlist = Ngl.Resources()
    rlist.wkWidth  = 3000 #page resolution
    rlist.wkHeight = 3000
    wks = Ngl.open_wks(wks_type,image_name,rlist)
    res.nglMaximize = True
    res.nglDraw  = False; res.nglFrame = False
    res.vpWidthF = 0.74; res.vpHeightF = 0.2#0.54
    res.tmBorderThicknessF = 10.0
    res.tmXBMajorThicknessF = 10.0; res.tmXBMinorThicknessF = res.tmXBMajorThicknessF
    res.tmYLMajorThicknessF = res.tmXBMajorThicknessF; res.tmYLMinorThicknessF = res.tmXBMajorThicknessF
    res.tmXBLabelFont = 26; res.tmYLLabelFont = res.tmXBLabelFont
    res.tmXBLabelFontHeightF = 0.016; res.tmYLLabelFontHeightF = res.tmXBLabelFontHeightF
    res.tiXAxisFontHeightF = 0.020; res.tiYAxisFontHeightF = res.tiXAxisFontHeightF
    res.tiMainFontHeightF = 0.018; res.tiMainFont = 26
    LeftString = ""
    print("\033[44m{} \033[0m".format(LeftString))
    return wks, res, LeftString

#--------------------------------------------------------------------
start_time = cpu_time()
#initialize parameters
samplerate = 1000 # Hz, get 1000 samples per sec. (or kilometer)
N = 1024 # datalength, total samples
zero_padd_log = True
windowlength = 128 # data points, the datalength in the certain window
#--------------------------------------------------------------------
sinefreq1 =  50 # Hz
sinefreq2 = 200 # Hz
SNR = -7 # dB, Signal-to-noise ratio

#generate simulated signals
t = np.arange(1,N+1)/samplerate; #print(t)
sig1 = np.sin(2.*np.pi*sinefreq1*t)
sig2 = np.sin(2.*np.pi*sinefreq2*t)
#Additive White Gaussian Noise (AWGN)
data = awgn(sig1+sig2,SNR)

""" Spectral analysis (FFT) """
#the max. resolved frequency is half sampling rate
if zero_padd_log:
   #Zero-padding could make it faster
   nfft = int(2**nextpow2(N)) # Next power of 2 from length of y 
   faxis = samplerate/2. * np.linspace(0,1,int(nfft/2)+1) # Hz
   data_freq = np.fft.fft(data,nfft) 
else:
   faxis = samplerate/2 * np.linspace(0,1,int(N/2)+1) # Hz
   data_freq = np.fft.fft(data) 
print("the min. resolved period is 2 * {}({}).".format(t[1]-t[0],1/faxis[-1]))
PS = abs(data_freq)**2.
PS = PS/max(PS) # normalize PS
print(PS.shape)

""" Spectral analysis (FFT with Welch method) """
# power spectrum, via scipy welch. 'boxcar' means no window, nperseg=len(y) so that fft computed on the whole signal.
(faxis_W,PS_W) = scipy.signal.welch(data, samplerate, window='hamming',nperseg=windowlength,scaling='spectrum', \
                 axis=-1, average='mean', detrend=False)
PS_W = PS_W/max(PS_W) # normalize PS
print(PS_W.shape)

""" Plot 1 """
image_name = "FFT_ex3_1"
res = Ngl.Resources()
(wks,res,dummy) = wks_setting(res,image_name)
res.tiXAxisString = "~F25~Time (s)"
res.tiYAxisString = "~F25~"
res.xyLineColors = ["blue"]
res.xyLineThicknesses = [20]
res.trXMinF = t[0]    # Limits for X axis
res.trXMaxF = t[len(t)//4]#t[-1]
res.tmXBMode   = "Explicit"
res.tmXBValues = np.arange(0,0.25+0.05,0.05); #print(res.tmXBValues)
res.tmXBLabels = res.tmXBValues; #print(res.tmXBLabels)
res.tmXBMinorValues = np.arange(0.01,t[-1],0.01)
res3 = copy.deepcopy(res)

res.trYMinF = -1.0    # Limits for Y axis
res.trYMaxF =  1.0
res.tmYLMode   = "Explicit"
res.tmYLValues = np.round(np.arange(-1.,1.+0.2,0.2),1); #print(res.tmYLValues)
ylabels = res.tmYLValues; ylabels[(len(ylabels)-1)//2] = 0.0
res.tmYLLabels = ylabels; #print(res.tmYLLabels)
res.tmYLMinorValues = np.arange(-0.9,0.9+0.2,0.2)

leftname = "Data from {} to {} seconds".format(t[0],t[-1])
res.tiMainString = str(sinefreq1)+"-Hz Sine Wave"+ \
                   "~C~~B~"+leftname

res2 = copy.deepcopy(res)
res2.xyLineColors = ["red"]
res2.tiMainString = str(sinefreq2)+"-Hz Sine Wave"+ \
                   "~C~~B~"+leftname

res3.trYMinF = -10.0    # Limits for Y axis
res3.trYMaxF =  10.0
res3.xyLineColors = ["navy"]
res3.tiMainString = "Simulated Signal of Mixed Sine Waves With White Noise, "+ \
                    "~C~SNR = {} dB".format(SNR)
plot = []
plot.append(Ngl.xy(wks,t,sig1,res))
plot.append(Ngl.xy(wks,t,sig2,res2))
plot.append(Ngl.xy(wks,t,data,res3))

panelres = Ngl.Resources()
Ngl.panel(wks,plot[0:3+1],[3,1],panelres)
Ngl.destroy(wks); del res, panelres, plot
""" Plot 2-1 """
image_name = "FFT_ex3_2"
resf = Ngl.Resources()
(wks,resf,dummy) = wks_setting(resf,image_name)
resf.xyLineColors = ["navy"]
resf.xyLineThicknesses = [25]
resf.vpWidthF = 0.74; resf.vpHeightF = 0.54
resf.tiMainFontHeightF = 0.025
resf.tiXAxisString = "~F25~Frequency (Hz)"
resf.tiYAxisString = "~F25~Power Spectrum (dB)"
leftname = "~F25~Period from ~F34~%~F25~ to {} seconds".format(1./faxis[-1])
resf.tiMainString = "Spectral Analysis (FFT)"+ \
                    "~C~~B~"+leftname

resf.trYMinF = -100    # Limits for Y axis
resf.trYMaxF =    0

resf.trXMinF = faxis[0]    # Limits for X axis
resf.trXMaxF = faxis[-1]
""" Plot 2-2 """
resW = copy.deepcopy(resf)
resW.trYMinF = -30    # Limits for Y axis
resW.trYMaxF =   0
resW.trXMinF = faxis_W[0]    # Limits for X axis
resW.trXMaxF = faxis_W[-1]
resW.xyLineColors = ["darkgreen"]
leftname = "~F25~Period from ~F34~%~F25~ to {} seconds".format(1./faxis_W[-1])
resW.tiMainString = "Spectral Analysis (FFT with Welch method)"+ \
                    "~C~~B~"+leftname

plotf = []
PS = PS[range(int(len(PS)/2)+1)]; print("y:{}, x:{}".format(PS.shape,faxis.shape))
PS = 20.*np.log10(PS) # dB
print("Welch method, y:{}, x:{}".format(PS_W.shape,faxis_W.shape))
PS_W = 20.*np.log10(PS_W) # dB
plotf.append(Ngl.xy(wks,faxis,PS,resf))
plotf.append(Ngl.xy(wks,faxis_W,PS_W,resW))
panelres = Ngl.Resources()
Ngl.panel(wks,plotf[0:1+1],[2,1],panelres)
Ngl.destroy(wks); del resf, resW, panelres, plotf

Ngl.end()
end_time = cpu_time(); end_time = (end_time - start_time)/60.0
print(os.path.basename(__file__)+" has done!\nTime elapsed: {:.2f}".format(end_time), "mins.")
