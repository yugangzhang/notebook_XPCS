import sys, shutil
EssCod_PATH='/home/yuzhang/XPCS_Anlysis/notebook_XPCS/EssCode/'
EssCod_PATH in sys.path or sys.path.append(EssCod_PATH)
 

from Get_Data import readframe_series
from Process_Data import average_img, cpdump,cpopen 
from Plot_Data import plot,show_img,show_mask,show_img_with_ROI,pd_data_plot
from XPCS_Functions import make_qlist, calqlist,get_trace_and_maxc
from XPCS_Functions import get_pixellist_intensity,azimuthal_integration
from XPCS_Functions import calqlist_regions, get_waterfall
from XPCS import xpcs

from XSVS import xsvs
from Plot_XSVS import xsvs_plot_histogram


import numpy as np
import matplotlib.pyplot as plt


from Setup_file import *

ave= np.load(outDir + 'sid_%s-ave.npy'%sid)


dimy,dimx = img.shape
qmask = mask

#noqs = int( float(qend - qstart)/qwidth +1)
#qend = qstart + noqs*qwidth
#print noqs,qend
qlist,qradi = make_qlist(qstart,qend, qwidth, noqs)
pixellist,qind,nopr,nopixels = calqlist(qlist,qradi, dimx,dimy, cenx,ceny,qmask)


    
