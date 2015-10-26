# Import Nessary Modules 
###################################################################
import sys, shutil
EssCod_PATH='/home/yuzhang/XPCS_Anlysis/notebook_XPCS/EssCode/'
EssCod_PATH in sys.path or sys.path.append(EssCod_PATH)


from Get_Eiger_Pixel_Mask import load_eiger_mask
from setupQ import *
from numpy import array,load,arange
import os
T=True
F=False
########################
'''This is the parameter input file for the visibility analysis
    the parameters include:
    data_path ( DATA_DIR, RES_DIR)
    FILENAME( dataPref,outPrefix )
    frame_format (TIFF, EDF,...)
    frame_dimension (dimx,dimy)
    frame_center (cenx,ceny)
    frame to be processed (begframe,noframes,)
    mask
    DK (for the dark image)
    qrange (qstart,qend,qwidth,noqs,
    
    dpix (pixel size),lambda_ (wavelength),Ldet (sample to dector distance)
    for transform pixel to real q values
    
    dead_time,exposuretime
    
    exposureList,imgtobin,noframe_fast for XSVS
    nobuf, nolev for XPCS analysis
'''
 


inDir= '/XF11ID/data/2015/10/16/'
outDir='/home/yuzhang/XPCS_Anlysis/Results/'

#sid = 'series4M_11_master.h5' 
#FILENAME = 'coralPor_10Hz_'  

sid = 'series4M_19_master.h5' 
FILENAME = 'coralPor_400Hz_'  

sid = 'series4M_13_master.h5' 
FILENAME = 'coralPor_100Hz_'  

sid = 'series4M_12_master.h5' 
FILENAME = 'coralPor12_100Hz_'  

sid = 'series4M_11_master.h5' 
FILENAME = 'coralPor_10Hz_' 



sid = 'series4M_20_master.h5' 
FILENAME = 'coralPor20_250Hz_' 



sid = 'series4M_21_master.h5' 
FILENAME = 'coralPor21_100Hz_' 

sid = 'series4M_19_master.h5' 
FILENAME = 'coralPor_400Hz_'  


sid = 'series4M_14_master.h5' 
FILENAME = 'coralPor_400Hz_'  


sid = 'series4M_11_master.h5' 
FILENAME = 'coralPor11_10Hz_' 


#mask= load( outDir+ 'mask3.npy') #use hotspot >90
#mask= load( outDir+ 'mask4.npy') #for 1400 Hz
mask1= load( outDir+ 'mask2.npy') #use hotspot >90
maskps = load_eiger_mask( sid, inDir)
mask = mask1 | maskps 
 


DK= None
FK= None

#outDir='/home/yugang/CHX_Analysis/Res/2220_50ms/'
#if not os.path.exists(outDir): os.makedirs(outDir)


dimx= 2070  #1030  #1030;
dimy= 2167  #1065  #1065

ceny =  1634.66  #608.0
cenx =  838.6    #820.0

qstart = 24  #20  #18 #48#35 #23  #4          #27 for one ring #24

qstart= 110
#qstart= 110
qend = 520 #60  #36  #56
qwidth = 2  #2
noqs = 10
qlist_= None


begframe = 0 #  
noframes = 1000

#noframes = 100
tmax=noframes

tmax=noframes * 1.0

dpix=0.075  #for eiger 1M and 4M the pixel size is 75X75 um
 
lambda_ = get_Lambda(E = 9, u='A')
qperpixel =qpix(dpix, lambda_=lambda_, Ldet =5000.) #8000 eV
 


exposuretime= 100 * 10**(-3) #unit in second, 100ms, 50 ms
acquisition_period = 100 * 10**(-3)  #10 second
#deadtime= 0   # 60e-6;
#timeperframe = exposuretime+deadtime
timeperframe = acquisition_period  
#################
#for XSVS########
#################

noframe_fast =noframes
#noframe_fast = 1000
flatfield=  FK



#################
#for XPCS########
#################
nobuf= 8  #//must be even!
nolev=  12

#nobuf= 4  #//must be even!
#nolev=  3

xbar=cenx;ybar=ceny

 
#################
#for XSVS########
#################
 
sids = { 's1': [sid, 1.], 
        
          }



print ('Sid:  %s is in processing.'%sid)


    
