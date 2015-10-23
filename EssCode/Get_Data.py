# -*- coding: utf-8 -*-
######################################################################################
#####################Read  Images ##########################################
#########By Dr. Yugang Zhang, yuzhang@bnl.gov, Octa 18, 2015##############################################
#######################################################################################
######################################################################################
############################################################################################################
#######################################################################################
######################################################################################
################################################################################

'''Functions:
1) readframe_series: load_data


'''
#export HDF5_PLUGIN_PATH=/home/yugang/anaconda/envs/ophyd/lib
######################################
##############readframe_series using pims_reader

def readframe_series(  FileName ='lio3_I_201_master.h5',
                       DataDir ='/XF11ID/data/2015/7/27/'):
    '''Read Eiger frames
       DataDir is the file path, e.g., '/XF11ID/data/2015/7/27/'
       FileName is mster data name, e.g., 'lio3_I_201_master.h5' 
    '''
    #from chxtools.pims_readers import EigerImages as Images
    from eiger_io.pims_reader import EigerImages as Images

    return Images(  DataDir + FileName )








