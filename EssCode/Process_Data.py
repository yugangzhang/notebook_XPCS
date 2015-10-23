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
   1) average_img: do image average
   2) cpdump( data, filename, outDir=None)
   3) cpopen(  filename=None, inDir=None,  )  
   4) trans_data_to_pd(data, label = None, dtype = 'list')   
   5) select_region(img, vert, keep_shape=True, qmask=None,  )
   6) find_max_inten( imgs, thres=10 )


   '''

##average_images

def average_img( imgs, Ns=None,Ne = None ):
    ''' Do imgs average,
        Optiions:
        imgs: the image seriers
        Ns: the start image
        Ne: the last image
        e.g.,
        ave = average_img(imgs)'''
    import numpy as np 
    ave = np.zeros_like(imgs[0],dtype =float)
    if Ns is None:Ns=0
    if Ne is None:Ne=len(imgs)
    for i in range(Ns,Ne):
        ave += imgs[i]
    ave /= (Ne-Ns)
    return ave

def cpdump( data, filename, outDir=None):
    import cPickle
    if outDir!=None:filename=outDir + filename
    fp=file(filename,'wb')   
    cPickle.dump(data,fp)
    fp.close()
    

def cpopen(  filename=None, inDir=None,  ): 
    import cPickle,os    
    if inDir!=None:filename=inDir + filename
    if os.path.isfile(filename):
        fp=file(filename,'rb')   
        data = cPickle.load(fp)
        fp.close()
        return data
    else:
        return None




############################################################################
##data format conversion

def trans_data_to_pd(data, label = None, dtype = 'list'):
    '''
    convert a data to a panda.dataframe
    '''
    
    from numpy import arange, array
    import pandas as pd
    import sys
    if dtype == 'list':
        data = array(data).T
    elif dtype == 'array':
        data = array(data)
    else:
        print "Wrong data type! Now only support 'list' and 'array' tpye"
    (N, M,) = data.shape
    index = arange(N)
    if label is None:
        label = [ 'data%s' % i for i in range(M) ]
    df = pd.DataFrame(data, index=index, columns=label)
    return df 



 




############################################################################
##get the interested region

def select_region(img, vert, keep_shape=True, qmask=None,  ):
    '''Get a pixellist by a rectangular region
        defined by
        verts e.g. xs,xe,ys,ye = vert #x_start, x_end, y_start,y_end
        (dimy, dimx,) = img.shape
       Giving cut postion, start, end, width
     if keep_shape:the output img don't change shape and
     make the non-interested region zeroes
     else: only return the interested region

    '''
    import numpy as np

    xs,xe,ys,ye = vert
    if keep_shape:       
        img_= np.zeros_like( img )
        #img_= np.zeros( [dimy,dimx])
    
        try:
            img_[ys:ye, xs:xe] = True
        except:
            img_[ys:ye, xs:xe,:] = True
        pixellist_ = np.where( img_.ravel() )[0]
        #pixellist_ =  img_.ravel()
        if qmask is not None:
            b=np.where(  qmask.flatten()==False )[0]
            pixellist_ = np.intersect1d(pixellist_,b)
        #imgx = img[pixellist_]
        #imgx = imgx.reshape( xe-xs, ye-ys)
        imgx = img_.ravel()    
        imgx[pixellist_] = img.ravel()[pixellist_]   
        imgx = imgx.reshape(  img.shape )
        
    else:         
        try:
            imgx =img[ys:ye, xs:xe]
        except:
    
            imgx =img[ys:ye, xs:xe,:]
        
    return  imgx


############################################################################
##get the max_intensity positoin

def find_max_inten( imgs, thres=10 ):
    '''find the pos (x,y) of the max intensity for each frame
        imgs is a pims.frames
        thres is the threshold for the max_intensity'''
    import numpy as np
    N = len(imgs)
    pos =[]  #[hor_cut, vert_cut]
    for i in range(N):
        img = imgs[i]
        m = img.max()
        if m<=thres:pos_=[0,0]
        else:
            x,y = np.where( img==img.max())
            pos_ = [x[0],y[0]]
        pos.append( pos_ )
    pos = trans_data_to_pd( pos, ['x','y'],'array' )
    return pos    
    







