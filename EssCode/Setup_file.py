from setupQ import *
from numpy import array,load,arange
import os
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
 

inDir= '/home/yugang/Desktop/Data/CoralPor/'
outDir='/home/yugang/CHX_Analysis/Res/'
sid = 'coralPor_2220_master.h5'
sid = 'coralPor_2232_master.h5'  #for 1 ms
sid = 'coralPor_2224_master.h5'  #for 5 ms
sid = 'coralPor_2220_master.h5'  #for 50 ms
  

sid = 'coralPor_2220_master.h5'   
FILENAME = 'coralPor_2220_'  



mask= load( outDir+ 'mask2.npy')
DK= None
FK= None

outDir='/home/yugang/CHX_Analysis/Res/2220_50ms/'
if not os.path.exists(outDir): os.makedirs(outDir)


dimx= 1030  #1030;
dimy= 1065  #1065

cenx =  608.0
ceny = 820.0

qstart = 24  #20  #18 #48#35 #23  #4          #27 for one ring #24
qstart=33
qend = 120 #60  #36  #56
qwidth = 2  #2
noqs = 10
qlist_= None


begframe = 0 #1
noframes = 3000  
noframes = 1000 
#noframes = 100
tmax=noframes

tmax=noframes * 1.0

dpix=0.055
 
lambda_ = get_Lambda(E = 12.3, u='A')
qperpixel =qpix(dpix, lambda_=lambda_, Ldet =1400.) #8000 eV
 


exposuretime= 50 * 10**(-3) #50 ms

deadtime=60e-6;
timeperframe = exposuretime+deadtime

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
nolev=  9

#nobuf= 4  #//must be even!
#nolev=  3

xbar=cenx;ybar=ceny

 
#################
#for XSVS########
#################

#sid = 'coralPor_2224_master.h5' 
sids = { 's1': [sid, 1.],            
          }



print 'Sid:  %s is in processing.'%sid


    
