import numpy as np
from wrf import getvar,ALL_TIMES,get_pyngl,latlon_coords,to_np, \
                interplevel
import Ngl, Nio
import glob
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor

#Li-Hsin Chen 2020

def info(ff):
    var = getvar(ff, "QCLOUD", timeidx=0)
    nz = np.int32(len(var[:,0,0]))
    ny = np.int32(len(var[0,:,0]))
    nx = np.int32(len(var[0,0,:])) 
    (lat,lon) = latlon_coords(var) # Get the latitude and longitude points
    lat = to_np(lat); lon = to_np(lon)
    #Set some map options based on information in WRF output file
    res = get_pyngl(var); del var 
    res.tfDoNDCOverlay = True # required for native projection
    t = getvar(ff, "times", timeidx=ALL_TIMES) #-1 isn't work
    nt = np.int32(len(t))
    print("Domain:{}".format(nt)+"|{}".format(nz)+\
    "|{}".format(ny)+"|{} ".format(nx)+"."*6+"Get Dimensions!! \
    \nGet WRF map-related resources and coordinates(lat. & long.)")
    return nt, nz, ny, nx, t, lat, lon, res

def search_wrf_file(Path="./"):
    FileID = glob.glob(Path+"*wrf*:*") #just like "ls *wrf*:*"
    for i, Filename in enumerate(FileID): print("\033[92m({}){}\033[0m".format(i,Filename))
    ID = input("\033[92mInput File ID[0-{}]: \033[0m".format(len(FileID)-1))
    if ID.replace(".","",1).isdigit() and float(ID) >= 0.0: #integer or float True
       File = FileID[int(round(float(ID),0))]
    else:
       File = FileID[0] #Default: 0
       print("\033[91m\tCan't find ur input ID (Default: 0)\033[0m")
    print("\033[103m\033[30mThe choosen file: "+File+"\033[0m")
    return File

def add_files(File,Path="./"):
    add_file_log = input("\033[92mAdd another experiment ? (Y/N) \033[0m")
    while add_file_log == "Y":
          File.append(search_wrf_file(Path))
          add_file_log = input("\033[92mAdd another experiment ? (Y/N) \033[0m")
    return File

def check_plt_time(nt,times):
    check_time = "N"
    while check_time != "Y":
          tid = []
          tid.append(input("\033[92mSelect start Time [0-{}]: \033[0m".format(nt-1)))
          if tid[0].replace(".","",1).isdigit() == False: tid[0] = 0 #Default
          tid.append(input("\033[92m  Select end Time [{}-{}]: \033[0m".format(tid[0],nt-1)))
          if tid[1].replace(".","",1).isdigit() == False: tid[1] = tid[0] #Default
          tid = [np.int(i) for i in tid]
          if len(times) <= 30:
             for i, time in enumerate(times):
                 time_print = "\t({:02d}){}".format(i,str(time.values)[0:19])
                 if i >= tid[0] and i <= tid[1]: print(time_print)
                 else: print("\033[90m"+time_print+"\033[0m")
          else:
             for i in range(30):
                 time_print = ""
                 for j,time in enumerate(times):
                     if j % 30 == i:
                        if j >= tid[0] and j <= tid[1]:
                           time_print = time_print+"\t({:02d}){}".format(j,str(time.values)[0:19])
                        else: 
                           time_print = time_print+"\033[90m\t({:02d}){}\033[0m"\
                                        .format(j,str(time.values)[0:19])
                 print(time_print)

          check_time = input("\033[92mConfirm for the selected time ? (Y/N) \033[0m").strip()
          check_time = np.where(check_time == "y","Y",check_time)
    return tid

#Shapefile reference: http://www.diva-gis.org/gdata
#                     https://data.gov.tw/dataset/7442
#                     https://data.gov.tw/dataset/13795
#Local Path: /home/WRF/shapefile/(Shimen/)
def add_shapefile_polylines(ff,wks,plot,color="black",thick=10.0):
    """ Attach shapefile polylines to map """
    f_shap = Nio.open_file(ff, "r")
    lon = f_shap.variables["x"][:] #np.ravel()
    lat = f_shap.variables["y"][:]
    lnres = Ngl.Resources()
    lnres.gsLineColor = color
    lnres.gsLineThicknessF = thick
    lnres.gsSegments = f_shap.variables["segments"][:,0]
    return Ngl.add_polyline(wks, plot, lon, lat, lnres)

def get_shp_traces(ff,color="black",thick=10.0):
    """ for 'plotly' attach shapefile polylines to map (gsn_code.ncl) """
    # init. plotting list
    data = dict(
        x=[],
        y=[],
        mode='lines',
        line=dict(color=color),
        name=' '
    )
    f_shap = Nio.open_file(ff, "r")
    #print(f_shap)
    #---Read global attributes
    geom_segIndex = f_shap.geom_segIndex
    geom_numSegs  = f_shap.geom_numSegs
    segs_xyzIndex = f_shap.segs_xyzIndex
    segs_numPnts  = f_shap.segs_numPnts

    lon = f_shap.variables["x"][:] #np.ravel()
    lat = f_shap.variables["y"][:]
    segments = f_shap.variables["segments"][:,:]
    geometry = f_shap.variables["geometry"][:,:]
    numFeatures = len(geometry[:,0])
    for i in range(numFeatures):
        startSegment = np.int(geometry[i, geom_segIndex])
        numSegments  = np.int(geometry[i, geom_numSegs])
        #print("{}, {}, {}".format(i,startSegment, numSegments+1))
        for s in range(startSegment, startSegment + numSegments):
            startPT = np.int(segments[s, segs_xyzIndex])
            endPT   = startPT + np.int(segments[s, segs_numPnts]) - 1
            lat_cc = lat[startPT:endPT]
            lon_cc = lon[startPT:endPT]
            data['x'] = data['x'] + lon_cc.tolist() + [np.nan]
            data['y'] = data['y'] + lat_cc.tolist() + [np.nan]
    return [data]

def plt_marker(wks,plot,lat,lon,log,jc=50,ic=50,\
                   idx=12,sz=20.0,tk=10.0,cr="black"):
    """ Add some polymarkers showing the original locations of the X,Y points."""
    def add_polymarker(wks,plot,xd,yd,index,size,thick,color):
        poly_res = Ngl.Resources()
        poly_res.gsMarkerIndex = index #4   # choose circle as polymarker
        poly_res.gsMarkerSizeF = size       # select size to avoid streaking
        poly_res.gsMarkerThicknessF = thick #30.0
        #poly_res.gsMarkerOpacityF = 0.3
        poly_res.gsMarkerColor = color #purple4 # choose color
        return Ngl.add_polymarker(wks, plot, xd, yd, poly_res)
    #log: if plot by doing loop, only ask first time
    check_plt = ""
    if log == 0:
       #Default: N
       check_plt = input("\033[92mPlot Marker at ({},{}) ? (Y/N) \033[0m".format(jc,ic))
       check_plt = np.where(check_plt == "y","Y",check_plt)
    if check_plt == "Y" or log == 1: 
       xd = lon[jc,ic]; yd = lat[jc,ic]
       plot_marker = add_polymarker(wks,plot,xd,yd,idx,sz,tk,cr)
       log = 1
    else:
       plot_marker = None 
       log = 2
    return log,plot_marker
 
def target_area(res,lat,lon,jc=50,ic=50,jrg=20,irg=20):
    plt_target_area = False #Default: N
    check_plt_area = input("\033[92mPlot target area ? (Y/N) \033[0m")
    check_plt_area = np.where(check_plt_area == "y","Y",check_plt_area)
    if check_plt_area == "Y": 
       pt = [jc,ic] 
       jj = np.array([pt[0]-jrg, pt[0]+jrg]).astype(int) 
       ii = np.array([pt[1]-irg, pt[1]+irg]).astype(int)
       plt_target_area = True
    else:
       jj = np.array([0, len(lat[:,0])-1]).astype(int)
       ii = np.array([0, len(lat[0,:])-1]).astype(int)
    # Zoom in on map area of interest
    res.mpLimitMode = "Corners" 
    res.mpLeftCornerLatF  = lat[jj[0],ii[0]] #bottom y left
    res.mpLeftCornerLonF  = lon[jj[0],ii[0]] #bottom x left
    res.mpRightCornerLatF = lat[jj[1],ii[1]] #Top y left
    res.mpRightCornerLonF = lon[jj[1],ii[1]] #Top x left
    print("start point: ({:3d}){:.3f}, ({:3d}){:.3f}\
          \n  end point: ({:3d}){:.3f}, ({:3d}){:.3f}"\
     .format(jj[0],res.mpLeftCornerLatF,ii[0],res.mpLeftCornerLonF, \
             jj[1],res.mpRightCornerLatF,ii[1],res.mpRightCornerLonF))
        
    return res, ii, jj, plt_target_area

#Test ProcessPoolExecutor & ThreadPoolExecutor
def interp_z(ff,time,varname,thread=True,cores=7):
    """ thread = True may be faster"""
    # global variables 
    z = getvar(ff, "z", timeidx=time) 
    var = getvar(ff, str(varname), timeidx=time)
    """ interpolation to vertical coordinate """
    def delta_z(kbot=100,nk=150):
        kk = kbot # m bottom
        # top 15.0 km
        lev = np.zeros(nk,dtype=np.float32)
        for k in range(nk):
            lev[k] = kk; #print(lev[k])
            kk += 100.0
        return lev, nk

    (lev,nk) = delta_z() # dz for interpolation
    ny = np.int32(len(var[0,:,0]))
    nx = np.int32(len(var[0,0,:]))
    var_m = np.zeros(shape=(nk,ny,nx), dtype=np.float32)
    print("Interpolation dimensions (k,j,i): \033[93m{}\033[0m".format(var_m.shape))
    if thread:
       with ThreadPoolExecutor(max_workers=cores) as executor:
            for k in range(nk):
                var_m[k,:,:] = executor.submit(interplevel,var,z,lev[k]).result()
    else:
       for k in range(nk):
           var_m[k,:,:] = interplevel(var, z, lev[k])
   
    return var_m

#-- Routine ngl_Strings: draw left, right and/or center string
def ngl_Strings(wks, plot, left='', center='', right=''):
    """
    *Reference: https://github.com/NCAR/pyngl/issues/11
       ngl_Strings(wks, plot, left='', center='', right='')

       Add annotations
	- left, right or center string above plot
			
       Correspond to NCL's 
	  gsnLeftString, gsnCenterString, gsnRightString'
    """
    assert str(getattr(wks,'__class__')  == "<class 'int'>"), 'ERROR - 1st parameter is not a Ngl wks'
    assert str(getattr(plot,'__class__') == "<class 'ngl.PlotIds'>"), 'ERROR - 2nd parameter is not a Ngl plot'
	
    vpx = Ngl.get_float(plot,"vpXF")         #-- retrieve value of res.vpXF from plot
    vpy = Ngl.get_float(plot,"vpYF")         #-- retrieve value of res.vpYF from plot
    vpw = Ngl.get_float(plot,"vpWidthF")     #-- retrieve value of res.vpWidthF from plot
	
    txres = Ngl.Resources()
    txres.txFontHeightF = 0.018              #-- font size for left, center and right string
    txres.txFont = 26
	
    y = vpy + 0.025                          #-- y-position
	
    if(left != ''):
       txres.txJust = "CenterLeft"           #-- text justification
       x = vpx                               #-- x-position
       Ngl.text_ndc(wks, left, x, y, txres)  #-- add text to wks
	   
    if(center != ''):
      txres.txJust = "CenterCenter"          #-- text justification	   x = vpx + vpw/2
      x = vpx + vpw/2
      Ngl.text_ndc(wks, center, x, y, txres) #-- add text to wks

    if(right != ''):
       txres.txJust = "CenterRight"          #-- text justification
       x = vpx+vpw                           #-- x-position
       Ngl.text_ndc(wks, right, x, y, txres) #-- add text to wks

def Date_string(yy,mm,dd,hh,mi,TW_LST=False,plotly=False):
    def det_mon_yr_add(y,m,d,h):
        #string to integer
        yy = int(y); mm = int(m)
        dd = int(d); hh = int(h)
        months_day = [31,28,31,30,31,30,31,31,30,31,30,31]
        # leap year or ordinary year
        if yy % 4 == 0 and yy % 100 != 0:
           months_day[2-1] = 29
        elif yy % 400 == 0:
           months_day[2-1] = 29
        #update hour & day & month
        if hh >= 24:
           hh  = hh - 24
           dd = dd + 1
           if dd > months_day[mm-1]:
              dd = dd - months_day[mm-1]
              mm = mm + 1
              if mm > 12:
                 mm = mm - 12
                 yy = yy + 1
        return str(yy), format(mm,'02d'), \
               format(dd,'02d'), format(hh,'02d')
    #--------------------------------------
    time_cord = "UTC"
    if TW_LST:
       time_cord = "LST"
       hh = str(int(hh) + 8) # UTC to LST, GMT+8
       (yy,mm,dd,hh) = det_mon_yr_add(yy,mm,dd,hh)
       #####yy = str(int(yy) - 1911) # the year of the Republic Era
       #print("{} LST {} {}, {}".format(hh,dd,mm,yy))

    #string to integer
    year = int(yy); mon  = int(mm)
    day  = int(dd); hr   = int(hh)
    months = ["Jan.","Feb.","Mar.","Apr.","May","June", \
              "July","Aug.","Sept.","Oct.","Nov.","Dec."]
    for i, month in enumerate(months):
        mm = np.where(mon == i + 1, month, mm)

    days_str = ["st","nd","rd"] # superscript
    sup = np.array([["~S~","~N~"],["<sup>","</sup>"]])
    s = 0
    if plotly: s = 1
    if day <= 3:
       for i, day_str in enumerate(days_str):
           if day == i+1: dd = dd+sup[s,0]+day_str+sup[s,1]
    else:
       dd = dd+sup[s,0]+"th"+sup[s,1]

    date = "{}{} {} {} {}, {}".format(hh,mi,time_cord,dd,mm,yy)

    print(date)
    return date

if __name__ == '__main__':
   info()
   search_wrf_file()
   add_files()
   check_plt_time()
   add_shapefile_polylines()
   get_shp_traces() # plotly
   plt_marker()
   interp_z()
   ngl_Strings()
   Date_string()
