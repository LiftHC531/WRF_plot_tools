# -*- coding: utf-8 -*-
import numpy as np

def cnMode(res,mode="ManualLevels",Min=0,Max=10,delta=1,\
           levels=[]):
    """ PyNgl setting 
        example:
        cnLevels = np.array([1,2,6,10,15,20,30,40,50,70,90,110,130,150,200,300], 'f') 
        res = colorbar_tables.cnMode(res,mode="ExplicitLevels",levels=cnLevels)
                                       or
        res = colorbar_tables.cnMode(res,mode="ManualLevels",Min=-5,Max=5,delta=1)
    """
    res.cnLevelSelectionMode = mode 
    if mode == "ManualLevels":
       res.cnMinLevelValF  = Min
       res.cnMaxLevelValF  = Max
       res.cnLevelSpacingF = delta
    elif mode == "ExplicitLevels":
       res.cnLevels = levels
    return res

def forplotly(colors):
    colors = np.round(colors*255).astype(int)
    plotly_type = []
    i = 0.0
    colorbin = 1.0/float(len(colors)-2-1)
    for j,color in enumerate(colors):
        if j > 1:
           #print("--{}".format(str(color)[1:12]))
           #print(i) # 0 ~ 1
           plotly_type.append([round(i,6),'rgb({}, {}, {})'\
                  .format(color[0],color[1],color[2])])
           i += colorbin
    print("--{}".format(plotly_type))
    colors = plotly_type
    return colors


def CWB_precip(plotly=False):
    """ 18 colors
       The color bar of rain for Central Weather Bureau, Taiwan.
       example: 
       Ngl.define_colormap(wks, colorbar_tables.CWB_precip())
       res.cnLevelSelectionMode = "ExplicitLevels" 
       res.cnLevels = numpy.array([0.1,1,2,6,10,15,20,30,40,50,70,90,110,130,150,200,300], 'f')
       #res.cnLevels = numpy.array([0.1,10,20,60,100,150,200,300,400,500,600,700,800,900,1000,1200,1500], 'f')
    """
    #The first two elements are background color and frame color
    colors = np.array([
     [1.000,1.000,1.000], [0.000,0.000,0.000], \
     [1.000,1.000,1.000], [0.757,0.757,0.757], [0.608,1.000,1.000], \
     [0.000,0.812,1.000], [0.039,0.596,1.000], [0.039,0.396,1.000], \
     [0.188,0.600,0.039], [0.196,1.000,0.000], \
     [0.973,1.000,0.000], [1.000,0.796,0.000], [1.000,0.603,0.000], \
     [0.980,0.012,0.000], [0.800,0.000,0.012], [0.627,0.000,0.000], \
     [0.596,0.000,0.604], [0.765,0.016,0.800], \
     [0.973,0.020,0.953], [0.996,0.796,1.000], \
    ], 'f')
    if plotly:
       colors = forplotly(colors) 
    return colors

def CWB_radar(plotly=False):
    """ 67 colors
       The color bar of radar reflectivity for Central Weather Bureau, Taiwan.
       Ngl.define_colormap(wks, colorbar_tables.CWB_radar())
       res.cnLevelSelectionMode = "ExplicitLevels" 
       res.cnLevels = numpy.arange(65+1) # 0, 1, ,2 ......, 65
    """
    #The first two elements are background color and frame color
    colors = np.array([
     [1.000,1.000,1.000], [0.000,0.000,0.000], \
     [1.000,1.000,1.000], [0.000,1.000,1.000], [0.000,0.925,1.000], \
     [0.000,0.854,1.000], [0.000,0.784,1.000], [0.000,0.714,1.000], \
     [0.000,0.639,1.000], [0.000,0.569,1.000], [0.000,0.498,1.000], \
     [0.000,0.427,1.000], [0.000,0.357,1.000], [0.000,0.282,1.000], \
     [0.000,0.212,1.000], [0.000,0.141,1.000], [0.000,0.071,1.000], \
     [0.000,0.000,1.000], [0.000,1.000,0.000], [0.000,0.957,0.000], \
     [0.000,0.914,0.000], [0.000,0.871,0.000], [0.000,0.827,0.000], \
     [0.000,0.784,0.000], [0.000,0.745,0.000], [0.000,0.706,0.000], \
     [0.000,0.667,0.000], [0.000,0.627,0.000], [0.000,0.588,0.000], \
     [0.200,0.671,0.000], [0.400,0.753,0.000], [0.600,0.835,0.000], \
     [0.800,0.918,0.000], [1.000,1.000,0.000], [1.000,0.957,0.000], \
     [1.000,0.914,0.000], [1.000,0.871,0.000], [1.000,0.827,0.000], \
     [1.000,0.784,0.000], [1.000,0.769,0.000], [1.000,0.722,0.000], \
     [1.000,0.659,0.000], [1.000,0.596,0.000], [1.000,0.533,0.000], \
     [1.000,0.471,0.000], [1.000,0.282,0.000], [1.000,0.188,0.000], \
     [1.000,0.094,0.000], [1.000,0.000,0.000], [0.957,0.000,0.000], \
     [0.914,0.000,0.000], [0.871,0.000,0.000], [0.827,0.000,0.000], \
     [0.784,0.000,0.000], [0.745,0.000,0.000], [0.706,0.000,0.000], \
     [0.667,0.000,0.000], [0.627,0.000,0.000], [0.588,0.000,0.000], \
     [0.671,0.000,0.200], [0.753,0.000,0.400], [0.835,0.000,0.600], \
     [0.918,0.000,0.800], [1.000,0.000,1.000], [0.918,0.000,1.000], \
     [0.835,0.000,1.000], [0.753,0.000,1.000], [0.671,0.000,1.000], \
     [0.588,0.000,1.000], \
    ], 'f')
    if plotly:
       colors = forplotly(colors) 
    return colors

def idlcorr(plotly=False):
    """ 16 colors
       The color bar of error correlation cite by Chung et al. 2013.
       'Examination of Situation-Dependent Background Error Covariances at
        the Convective Scale in the Context of the Ensemble Kalman Filter'
       example: Ngl.define_colormap(wks, colorbar_tables.idlcorr())
    """
    #The first two elements are background color and frame color
    colors = np.array([
     [1.000,1.000,1.000], [0.000,0.000,0.000], \
     [1.000,1.000,1.000], [0.329,0.560,0.784], \
     [0.286,0.447,0.710], [0.235,0.337,0.643], [0.220,0.247,0.502], \
     [0.224,0.710,0.373], [0.204,0.580,0.337], [0.227,0.424,0.267], \
     [0.953,0.922,0.250], [0.992,0.737,0.251], [0.953,0.522,0.239], \
     [0.894,0.247,0.263], [0.690,0.376,0.643], [0.612,0.329,0.624], \
     [0.596,0.243,0.533], [0.145,0.149,0.157], \
    ], 'f')
    if plotly:
       colors = forplotly(colors) 
    return colors

def NMCCMA_precip(plotly=False):
    """ 13 colors
       The color bar of rain for National Meteorological Center of
       China Meteorological Administration (NMCCMA, 中國中央氣象台).
       http://www.nmc.cn/
       example:
       Ngl.define_colormap(wks, colorbar_tables.NMCCMA_precip())
       res.cnLevelSelectionMode = "ExplicitLevels" 
       res.cnLevels = numpy.array([0,2.5,5,10,25,50,100,250,400,600,1000,1500], 'f')
    """
    #The first two elements are background color and frame color
    colors = np.array([
     [1.000,1.000,1.000], [0.000,0.000,0.000], \
     [1.000,1.000,1.000], [0.718,0.957,0.671], [0.447,0.839,0.435], \
     [0.220,0.733,0.231], [0.180,0.537,0.169], \
     [0.384,0.722,0.996], [0.000,0.000,0.996], [0.980,0.000,0.980], \
     [0.506,0.000,0.251], [1.000,0.667,0.004], [1.000,0.400,0.000], \
     [0.902,0.000,0.000], [0.600,0.004,0.000], \
    ], 'f')
    if plotly:
       colors = forplotly(colors) 
    return colors
def CWB_Vr(plotly=False):
    """ 29 colors
       The color bar of radar radial wind for Central Weather Bureau, Taiwan.
       Ngl.define_colormap(wks, colorbar_tables.CWB_Vr())
    """
    #The first two elements are background color and frame color
    colors = np.array([
     [1.000,1.000,1.000], [0.000,0.000,0.000], \
     [0.627,1.000,0.627], [0.333,0.853,0.333], [0.039,0.706,0.039], \
     [0.020,0.353,0.353], [0.000,0.000,0.667], [0.000,0.059,0.745], \
     [0.000,0.117,0.823], [0.000,0.196,0.902], [0.000,0.274,0.980], \
     [0.000,0.412,0.990], [0.000,0.549,1.000], [0.000,0.775,1.000], \
     [0.000,1.000,1.000], [0.490,0.990,0.990], [0.980,0.980,0.980], \
     [0.990,0.990,0.726], [1.000,1.000,0.471], [1.000,0.873,0.236], \
     [1.000,0.745,0.000], [1.000,0.510,0.000], [1.000,0.275,0.000], \
     [0.912,0.197,0.000], [0.824,0.118,0.000], [0.746,0.059,0.000], \
     [0.667,0.000,0.000], [0.824,0.079,0.490], [0.980,0.157,0.980], \
     [0.990,0.393,0.990], [1.000,0.628,1.000], \
    ], 'f')
    if plotly:
       colors = forplotly(colors)
    return colors

    
if __name__ == '__main__':
   cnMode()
   CWB_precip()
   CWB_radar()
   idlcorr()
   NMCCMA_precip()
   CWB_Vr()
