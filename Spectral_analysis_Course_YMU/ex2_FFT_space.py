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
""" Follow the ex1_FFT.py """

def wks_setting(res):
    image_name = "FFT_ex2"
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
samplerate = 1./3. # km-1, get 1 sample per 3 km
N = 150 # the size of x/y grid, total samples
#--------------------------------------------------------------------
sinfreq = 1./120. # km-1, wavenumber, i.e. wavelength = 120 km 
#generate a sine singals
x = np.arange(N)/float(samplerate); #print(x)
sig = np.sin(2.*np.pi*sinfreq*x); #print(sig)

""" Spectral analysis (FFT) """
#the max. resolved wavenumber is half sampling rate
faxis = float(samplerate)/2.*np.linspace(0,1,int(N/2)+1) # km-1

faxis_same = np.fft.fftfreq(x.size, x[1] - x[0]); #print(len(faxis_same)) 
faxis_same = abs(faxis_same[range(int(len(faxis_same)/2)+1)])

#print(1/faxis) 
#print(1/faxis_same)

sig_freq = np.fft.fft(sig)
PS = abs(sig_freq)**2.
PS = PS/PS.max() # normalize PS 

""" Plot """
res = Ngl.Resources()
(wks,res,title) = wks_setting(res)
res.xyLineColors = ["blue"]
res.xyLineThicknesses = [25]
res2 = copy.deepcopy(res)
res.trXMinF = x[0]    # Limits for X axis
res.trXMaxF = x[-1] 
res.tmXBMode   = "Explicit"
res.tmXBValues = np.arange(0,400+1,50); #print(res.tmXBValues)
res.tmXBLabels = res.tmXBValues; #print(res.tmXBLabels)
res.tmXBMinorValues = np.arange(0,x[-1]+1,10)

res.trYMinF = -1.0    # Limits for Y axis
res.trYMaxF =  1.0
res.tmYLMode   = "Explicit"
res.tmYLValues = np.round(np.arange(-1.,1.+0.2,0.2),1); #print(res.tmYLValues)
ylabels = res.tmYLValues; ylabels[(len(ylabels)-1)//2] = 0.0 
res.tmYLLabels = ylabels; #print(res.tmYLLabels)
res.tmYLMinorValues = np.arange(-0.9,0.9+0.2,0.2)

res.tiXAxisString = "~F25~x (km)"
res.tiYAxisString = "~F25~"
leftname = "Data from {} to {} km".format(x[0],x[-1])
res.tiMainString = str(1./sinfreq)+"-km wavelength Sine Wave"+ \
                   "~C~~B~"+leftname


res2.trXMinF = faxis[0]    # Limits for X axis
res2.trXMaxF = faxis[-1] 


res2.tiXAxisString = "~F25~Wavenumber (km~S~-1~N~)"
res2.tiYAxisString = "~F25~Power Spectrum"
leftname = "Wavelength from ~F34~%~F25~ to {} km".format(1./faxis[-1])
res2.tiMainString = "Spectral Analysis (FFT)"+ \
                    "~C~~B~"+leftname

plot = []
plot.append(Ngl.xy(wks,x,sig,res))
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

