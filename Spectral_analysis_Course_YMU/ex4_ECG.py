#-*- coding: utf-8 -*-
import sys
sys.path.append("..")
from package import wrf_tools,colorbar_tables #local function
from time import time as cpu_time
import os
import copy
import numpy as np
import scipy.signal
import glob
import Ngl, Nio
""" 
Follow the ./ECG/ECG_spectrum.m
https://github.com/scipy/scipy/issues/8045#issuecomment-337319294
"""

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
    res.vpWidthF = 0.74; res.vpHeightF = 0.3#0.54
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
zero_padd_log = True 
Path = "./ECG/"
#initialize parameters
samplerate = 500 # Hz, get 500 sample per 1 sec.
File = "ECG_500SampleRate.TXT" 
with open(Path+File, "r") as f:
     data = f.readlines()
     N = len(data)
     for i,datum in enumerate(data):
         datum = datum.replace("\n","")
         data[i] = float(datum) 
     print(data[0:10+1])
     print(N)
data = np.array(data)
#--------------------------------------------------------------------
""" raw data """
t = np.arange(1,N+1)/samplerate

"""
detect R peak and store R-R interval
perform double-threshold method to detect R peak
"""
# threshold 1: magnitude
#np.argwhere() --> find() in Matlab
indmags = np.argwhere(data >= 1.5)[:,0] # find the index with value larger than 1.5
nind = len(indmags)
print(indmags.shape)
print(indmags[1-1:nind-1+1]) # total
# threshold 2: define window
diffind = indmags[2-1:(nind-1)+1] - indmags[1-1:(nind-1-1)+1]
indgaps = np.argwhere(diffind > 1)[:,0]
print(indgaps.shape)
print(indgaps)
indmax = []   # the location of index with maximal value in each cycle(gaps)
for i in range(len(indgaps[:])+1): 
    if i == 0:
       # the first part in loop
       period = indmags[0:indgaps[0]]; #print(period) 
    elif i == len(indgaps[:])+1-1:
       # the last part in loop
       period = indmags[indgaps[i-1]+1:(nind-1)+1] 
    else:
       period = indmags[indgaps[i-1]+1:indgaps[i]]
    ind = np.argmax(data[period])
    if i < len(indgaps[:])+1-1:
       indmax.append(period[ind])

# heart rate variability (HRV)
#print(np.diff(indmax))
RRinterval = np.diff(indmax,axis=0)/samplerate
HRV = 60./RRinterval
# rewrite data set by HRV
data2 = RRinterval; print(data2.shape)
#data2 = HRV

""" Spectral Analysis (FFT) """
samplerate = 1./np.mean(RRinterval) # Hz
print("\033[95mget one sample of R-R interval per {:.3f} sec.\033[0m".format(samplerate))
N = len(RRinterval)
if zero_padd_log:
   nfft = int(2**nextpow2(N))
   print(nfft)
   faxis = samplerate/2. * np.linspace(0,1,int(nfft/2)+1)
   data2_freq = np.fft.fft(data2,nfft)
else:
   faxis = samplerate/2. * np.linspace(0,1,int(N/2)+1)
   data2_freq = np.fft.fft(data2)

PS = abs(data2_freq)**2.
PS = PS/max(PS)

""" Spectral analysis (FFT with Welch method) """
segmentNo = 20 #if larger, then smoother 
windowlength = int(round(N/segmentNo,0)); #print(windowlength)
# power spectrum, via scipy welch. 'boxcar' means no window, nperseg=len(y) so that fft computed on the whole signal.
(faxis_W,PS_W) = scipy.signal.welch(data2, samplerate, window='hamming',nperseg=windowlength,scaling='spectrum', \
              average='mean',detrend=False)
PS_W = PS_W/max(PS_W) # normalize PS
print(PS_W.shape)



""" Plot 1 x:t, y:data """
image_name = "ECG_ex4_1"
res = Ngl.Resources()
(wks,res,title) = wks_setting(res,image_name)
res.xyLineColors = ["blue"]
res.xyLineThicknesses = [15]
#res2 = copy.deepcopy(res)
res.trXMinF =  0 #t[0]    # Limits for X axis
res.trXMaxF = 10 #20 
res.tmXBMode   = "Explicit"
#res.tmXBValues = np.array([0.5,1.0,1.5,2.0]); #print(res.tmXBValues)
res.tmXBValues = np.arange(res.trXMinF,res.trXMaxF+2,2); #print(res.tmXBValues)
res.tmXBLabels = res.tmXBValues; #print(res.tmXBLabels)
res.tmXBMinorValues = np.arange(res.trXMinF+1,res.trXMaxF+1+2,2)

res.trYMinF = 0.8    # Limits for Y axis
res.trYMaxF = 2.2
res.tmYLMode   = "Explicit"
res.tmYLValues = np.arange(res.trYMinF,res.trYMaxF+0.2,0.2); #print(res.tmYLValues)
res.tmYLLabels = res.tmYLValues; #print(res.tmYLLabels)
res.tmYLMinorValues = np.arange(res.trYMinF+0.1,res.trYMaxF+0.1+0.2,0.2)

res.tiXAxisString = "~F25~Time (s)"
res.tiYAxisString = "~F25~"
leftname = "Data from {:.3f} to {:.3f} seconds".format(t[0],t[-1])
res.tiMainString = "Electrocardiogram (ECG)"+ \
                   "~C~~B~"+leftname

plot = []
plot.append(Ngl.xy(wks,t,data,res))

Ngl.draw(plot[0])
Ngl.frame(wks)
Ngl.destroy(wks); del res, plot
""" Plot 2 x:beats, y:R-R interval """
image_name = "ECG_ex4_2"
res = Ngl.Resources()
(wks,res,title) = wks_setting(res,image_name)
res.xyLineColors = ["blue"]
res.xyLineThicknesses = [15]
res.trYMinF = 0.5    # Limits for Y axis
res.trYMaxF = 1.3
res.tmYLMode   = "Explicit"
res.tmYLValues = np.arange(res.trYMinF,res.trYMaxF+0.1,0.1); #print(res.tmYLValues)
res.tmYLLabels = res.tmYLValues; #print(res.tmYLLabels)
res.tiXAxisString = "~F25~Beats(*{:.2f} sec.)".format(samplerate)
res.tiYAxisString = "~F25~R-R interval (s)"
leftname = "Beats from {} to {}".format(1,len(data2))
res.tiMainString = "Heart rate variability (HRV)"+ \
                   "~C~~B~"+leftname
plot = []
plot.append(Ngl.y(wks,data2,res))

Ngl.draw(plot[0])
Ngl.frame(wks)
Ngl.destroy(wks); del res, plot
""" Plot 3-1 x:frequency, y:Power spectrum """
image_name = "ECG_ex4_3"
resf = Ngl.Resources()
(wks,resf,title) = wks_setting(resf,image_name)
resf.xyLineColors = ["blue"]
resf.xyLineThicknesses = [15]
resf.vpWidthF = 0.74; resf.vpHeightF = 0.54
resf.tiMainFontHeightF = 0.025
resf.tiXAxisString = "~F25~Frequency (Hz)"
resf.tiYAxisString = "~F25~Power Spectrum (dB)"
leftname = "Period from ~F34~%~F25~ to {:.3f} seconds".format(1./faxis[-1])
resf.tiMainString = "Spectral Analysis (FFT)"+ \
                    "~C~~B~"+leftname
""" Plot 3-2 x:frequency, y:Power spectrum """
resW = copy.deepcopy(resf)
resW.xyLineColors = ["darkgreen"]
leftname = "Period from ~F34~%~F25~ to {:.3f} seconds".format(1./faxis_W[-1])
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
