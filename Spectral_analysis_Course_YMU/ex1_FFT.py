#-*- coding: utf-8 -*-
import sys
sys.path.append("..")
from package import wrf_tools #local function
from time import time as cpu_time
import os
import copy
import numpy as np
import glob
import Ngl, Nio
""" 
Follow the FFTspectrum.m
https://youtu.be/Kdz13vP0f-g?list=PLx_IWc-RN82uKOdafF4v4U5R_u4qmYaiu&t=509
"""
def nextpow2(p):
    n = 2
    while p > n:
       n *= 2
    return np.log2(n)

def abs2(x):
    return x.real**2 + x.imag**2

def wks_setting(res):
    image_name = "FFT_ex1"
    wks_type = "png"
    print("\033[44m---Plot: {}.{} ---\033[0m".format(image_name,wks_type))
    rlist = Ngl.Resources()
    rlist.wkWidth  = 3000 #page resolution
    rlist.wkHeight = 3000
    wks = Ngl.open_wks(wks_type,image_name,rlist)
    res.nglMaximize = True
    res.nglDraw  = False; res.nglFrame = False
    res.vpWidthF = 0.74; res.vpHeightF = 0.54
    res.tmBorderThicknessF = 10.0
    res.tmXBMajorThicknessF = 10.0; res.tmXBMinorThicknessF = res.tmXBMajorThicknessF
    res.tmYLMajorThicknessF = res.tmXBMajorThicknessF; res.tmYLMinorThicknessF = res.tmXBMajorThicknessF
    res.tmXBLabelFont = 26; res.tmYLLabelFont = res.tmXBLabelFont
    res.tmXBLabelFontHeightF = 0.016; res.tmYLLabelFontHeightF = res.tmXBLabelFontHeightF
    res.tiXAxisFontHeightF = 0.020; res.tiYAxisFontHeightF = res.tiXAxisFontHeightF
    res.tiMainFont = 25
    LeftString = ""
    print("\033[44m{} \033[0m".format(LeftString))
    return wks, res, LeftString

#--------------------------------------------------------------------
start_time = cpu_time()
#initialize parameters
samplerate = 250 # Hz, get 250 samples per sec. (or kilometer)
N = 512 # datalength, total samples
#print(nextpow2(1023))
zero_padd_log = True
#--------------------------------------------------------------------
sinfreq = 10 # Hz, frequency of the signal (or wavenumber km-1)
#generate a sine singals
t = np.arange(1,N+1)/float(samplerate); #print(t)
sig = np.sin(2.*np.pi*sinfreq*t); #print(sig)

""" Spectral analysis (FFT) """
#the max. resolved frequency is half sampling rate
if zero_padd_log:
   #Zero-padding could make it faster 
   nfft = int(2**nextpow2(N)) # Next power of 2 from length of y
   print(nextpow2(N))
   print(nfft)
   faxis = samplerate/2. * np.linspace(0,1,int(nfft/2)+1) # Hz
   sig_freq = np.fft.fft(sig,nfft)
else:
   faxis = float(samplerate)/2.*np.linspace(0,1,int(N/2)+1) # Hz
   sig_freq = np.fft.fft(sig)
print("the min. resolved period/wavelength is 2 * {}({}).".format(t[1]-t[0],1/faxis[-1]))
   
#PS = abs2(sig_freq)**2.
PS = abs(sig_freq)**2.
PS = PS/PS.max() # normalize PS 


""" Plot """
res = Ngl.Resources()
(wks,res,title) = wks_setting(res)
res.xyLineColors = ["blue"]
res.xyLineThicknesses = [25]
res2 = copy.deepcopy(res)
res.trXMinF = t[0]    # Limits for X axis
res.trXMaxF = t[-1] 
res.tmXBMode   = "Explicit"
res.tmXBValues = np.array([0.5,1.0,1.5,2.0]); #print(res.tmXBValues)
res.tmXBLabels = res.tmXBValues; #print(res.tmXBLabels)
res.tmXBMinorValues = np.arange(0.1,t[-1],0.1)

res.trYMinF = -1.0    # Limits for Y axis
res.trYMaxF =  1.0
res.tmYLMode   = "Explicit"
res.tmYLValues = np.round(np.arange(-1.,1.+0.2,0.2),1); #print(res.tmYLValues)
ylabels = res.tmYLValues; ylabels[(len(ylabels)-1)//2] = 0.0 
res.tmYLLabels = ylabels; #print(res.tmYLLabels)
res.tmYLMinorValues = np.arange(-0.9,0.9+0.2,0.2)

res.tiXAxisString = "~F25~Time (s)"
res.tiYAxisString = "~F25~"
leftname = "Data from {} to {} seconds".format(t[0],t[-1])
res.tiMainString = str(sinfreq)+"-Hz Sine Wave"+ \
                   "~C~~B~"+leftname


res2.trXMinF = faxis[0]    # Limits for X axis
res2.trXMaxF = faxis[-1] 


res2.tiXAxisString = "~F25~Frequency (Hz)"
res2.tiYAxisString = "~F25~Power Spectrum"
leftname = "Period from ~F34~%~F25~ to {} seconds".format(1./faxis[-1])
res2.tiMainString = "Spectral Analysis (FFT)"+ \
                    "~C~~B~"+leftname

plot = []
plot.append(Ngl.xy(wks,t,sig,res))
plot.append(Ngl.xy(wks,faxis,PS,res2))

#wrf_tools.ngl_Strings(wks, plot[0], left=leftname, center="", right="")

panelres = Ngl.Resources()

Ngl.panel(wks,plot[0:1+1],[2,1],panelres)  

#Ngl.draw(plot[0])
#Ngl.frame(wks)
Ngl.destroy(wks); del res, plot

Ngl.end()

end_time = cpu_time(); end_time = (end_time - start_time)/60.0
print(os.path.basename(__file__)+" has done!\nTime elapsed: {:.2f}".format(end_time), "mins.")

