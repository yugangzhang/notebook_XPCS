# -*- coding: utf-8 -*-
######################################################################################
#####################XPCS Functions ##########################################
#########By Dr. Yugang Zhang, yuzhang@bnl.gov, Octa 18, 2015##############################################
#######################################################################################
######################################################################################
############################################################################################################
#######################################################################################
######################################################################################
################################################################################
 

'''Functions:
   1) make_qlist( qstart,qend,qwidth,noqs,  )
   2) calqlist( qlist,qradi, dimx,dimy, cenx,ceny,
              qmask=None ,  shape='circle' ):to make qind, pixellist
   3) calqlist_region, select the interested region by giving vertices
   4) calqlist_regions, the interested regions 
   5) azimuthal_integration( image, qstart, qend, qwidth,
                          cenx,ceny, qmask=None,)
   6) get_waterfall(  imgs, pixellist,qind, fs=None, nf=None, 
                   q=None ) : calculate waterfall  
   7) get_trace_and_maxc: calculate trace of image series
   8) get_qRings(  qind, pixellist)
   9) get intensity_distribution(image, pixellist, qind)
   10) get_pixellist_intensity( image, center, q=10, qwidth=2, qmask=None,)
   11) find_center2(image, est_center, qs=10, qe=40,
                    qw=1,  var=5, qmask=None,)
   12) find_center(image, est_center, inner_radius=10, qwidth=1,  var=5, qmask=None,)
 

   '''


def make_qlist( qstart,qend,qwidth,noqs,  ):
    ''' DOCUMENT make_qlist( )
    give qstart,qend,qwidth,noqs
    return a qlist by giving the noqs, qstart,qend,qwidth.
    KEYWORD:  None    ''' 
    import numpy as np 
    qradi = np.linspace(qstart,qend,noqs)
    qlist=np.zeros(2*noqs) 
    qlist[::2]= np.int_(qradi-qwidth/2)  #render  even value
    qlist[1::2]= np.int_(qradi+(1+qwidth)/2) #render odd value
    #if qlist_!=None:qlist=qlist_
    return qlist,qradi

def calqlist( qlist,qradi, dimx,dimy, cenx,ceny,
              qmask=None ,  shape='circle' ):
    ''' DOCUMENT calqlist( qmask=,shape=, )
    calculate the equvilent pixel with a shape,
    return pixellist,nopixels,qind,nopr        
    KEYWORD:  qmask, a mask file;
        shape='circle', give a circle shaped qlist
        shape='column', give a column shaped qlist
        shape='row', give a row shaped qlist             
    '''       
    import numpy as np 
    y, x = np.indices( [dimy,dimx] )
    #print x.shape,y.shape
    #x and y have the same shape as img, which shape is [dimy,dimx]
    
    if shape=='circle':
        y_= y- ceny +1;x_=x-cenx+1        
        r= np.int_( np.hypot(x_, y_)    + 0.5  )#why add 0.5?
        #r= int_( hypot(x_, y_)     )
    elif shape=='column': 
        r= x
    elif shape=='row':
        r=y
    else:pass
    r= r.flatten() 
    #print qlist
    noqrs = len(qlist)    
    qind = np.digitize(r, qlist)        
    #if qmask==None:
    if qmask is  None:
        w_= np.where( (qind)%2 )# qind should be odd;print 'Yes'
        w=w_[0]
    else:
        a=np.where( (qind)%2 )[0]            
        b=np.where(  qmask.flatten()==False )[0]
        #print a,b
        #print a.shape,b.shape
        w= np.intersect1d(a,b)
        #w=a | b
        #w=b
    nopixels=len(w)
    qind=qind[w]/2
    #pixellist= (   y*imgwidth +x ).flatten() [w]
    pixellist= (   y*dimx +x ).flatten() [w]
    nopr,bins=np.histogram( qind, bins= range( len(qradi) +1 ))
    return pixellist,qind,nopr,nopixels


def calqlist_region(  vert,  dimx, dimy,qmask=None , ):
    '''Get a pixellist by a rectangular region
        Input:
            vert: a list, define the region,
                    provided as xs,xe,ys,ye, namely,
                    x_start, x_end, y_start,y_end
            dimx, dimy:
                    the shape of the image,
                    PAY ATTENTION: (dimy, dimx,) = img.shape
            qmask:
                    a optional mask
        Output:
           a pixel list '''
 
    import numpy as np

    ys,ye,xs,xe = vert
    #img_= np.zeros_like( img )
    img_= np.zeros( [dimy,dimx])

    img_[ys:ye, xs:xe] = True
    pixellist_ = np.where( img_.ravel() )[0]
    
    if qmask is not None:

        b=np.where(  qmask.flatten()==False )[0]
        pixellist_ = np.intersect1d(pixellist_,b)
  
    return pixellist_

def calqlist_regions(  verts,  dimx, dimy, qmask=None ,):
    '''Get a pixellist by a rectangular region
        Input:
            verts: a list, [vert1, vert2,...], define the regions,
                    each vert1 is also a list, provided as xs,xe,ys,ye, namely,
                    x_start, x_end, y_start,y_end
            dimx, dimy:
                    the shape of the image,
                    PAY ATTENTION: (dimy, dimx,) = img.shape
            qmask:
                    a optional mask
        Output:
           a pixel list and a index number '''
    import numpy as np
    pixellist = []
    qind=[]
    for i, vert in enumerate(verts):
        p = calqlist_region(  vert,  dimx, dimy, qmask )
        #print p.shape
        pixellist.append(  p  )
        qind.append( [i]*len(p))
    #pixellist = np.array(pixellist)
    return np.concatenate(pixellist),np.concatenate(qind) 


def azimuthal_integration( image, qstart, qend, qwidth,
                          cenx,ceny, qmask=None,):

    '''Do azimuthal_integration,
       image: 2d-data
       give: qstart, qend, qwithd, cenx, ceny, and qmask'''
    import numpy as np 
    dimy,dimx = image.shape
 
    noqs = int( float(qend - qstart)/qwidth +1)
    qend = qstart + noqs*qwidth
    #print noqs,qend
    qlist,qradi = make_qlist(qstart,qend, qwidth, noqs)
    pixellist,qind,nopr,nopixels = calqlist(qlist,qradi, dimx,dimy, cenx,ceny,qmask)
    #qinten = zeros( [noqs,2] )
    #qinten[:,0]=qradi
    inten=[]
    for n in xrange(noqs):
        value = np.ravel(image)[ pixellist[qind == n] ]
        #inten.append( value.mean() )
        L=len(value)
        if len(value)==0:L=1.0
        inten.append( value.sum()/L )
        
    #qinten[:,1]= inten    
    return  np.array(qradi),np.array(inten)


def get_waterfall(  imgs, pixellist,qind, fs=None, nf=None, 
                   q=None ):
    '''get waterfall,
        Giving: imgs: a pims.frames_sequences
                fs, nf: frame_start, number of frames
                pixellist, qind for the selected pixel
                if q is not None, only calculated the q-value waterfall
        Output: a dict{  q: np.array }'''

    import numpy as np

    uniq_q, nopr = np.unique(qind,return_counts=True )
    noqs = len(uniq_q)
    
    
    waterfall_dict={}
    if fs is None:fs=0
    if nf is None:nf = len(imgs)
    if nf>len(imgs):nf = len(imgs)
    nofram = nf
    #print nofra
    if q is None:            
        for qn in range(noqs):
            waterfall_dict['q%i'%qn]=np.zeros( [nofram, nopr[qn]])
    else:
        waterfall_dict['q%i'%q]=np.zeros( [nofram, nopr[q]])            
    for i in range(0, nofram ):            
        n= fs + i
        #print n
        #img_ = readframe_series( sid )
        img = imgs[n].ravel()                      
        img=img[pixellist]
        if q is not None:
            waterfall_dict['q%i'%q][i] = img[qind==q]
        else:
            for qn in range(noqs):
                waterfall_dict['q%i'%qn][i] = img[qind==qn]
    return waterfall_dict

    

def get_trace_and_maxc(  imgs, pixellist,qind, fs=None, nf=None, 
                 ):

    '''get a mean intensity of a qlist as function of frames
        Giving: imgs: a pims.frames_sequences
                fs, nf: frame_start, number of frames
                pixellist, qind for the selected pixel
                #if q is not None, only calculated the q-value trace
        Output: a pandas.DataFrame 

    '''
    import pandas as pd
    import numpy as np

    if fs is None:fs=0
    if nf is None:nf = len(imgs)
    if nf>len(imgs):nf = len(imgs)
    frame_start, nofram= fs,nf
    max_cts = 0
    
    trace_dict={} 
    
    for n in range(0,nofram  ): 
        img = imgs[n+ frame_start].ravel()        
        itd = intensity_distribution( img, pixellist,qind)        
        for k in itd.keys():
            c = itd[k].max()
            mean = itd[k].mean()
            if not trace_dict.has_key(k):trace_dict[k]=[]
            trace_dict[k].append( mean )
            if max_cts < c:
                max_cts = c
    #print (trace_dict.keys())
    keys = sorted( itd.keys() )            
    for k in keys:
        
        if k==min(keys):
            data = np.array(trace_dict[k]).reshape( nofram,1)
        else:
            data= np.hstack( [data, np.array(trace_dict[k]).reshape( nofram,1)   ]) 
    index=xrange( nofram )
    qn = len( trace_dict.keys())
    columns = [ 'q%i'%n for n in range(qn) ]
    trace = pd.DataFrame(data=data, index=index,columns=columns)
    return trace, int( max_cts  )+1



def get_qRings(  qind, pixellist):
    ''' give qind, pixellist, return the qRings
       qRings[0] = a pixel list
       ...
    '''
    import numpy as np
    noqs = len( np.unique( qind ) )
    qRings = list(np.zeros(noqs))
    for j in xrange(noqs):  #noqs is number of qrings
        qRings[j] = pixellist[np.where(qind==j)]
    return qRings

    



 
def intensity_distribution(image, pixellist, qind):
    '''Get a intensity distribution as a dict
       Giving: image, a two-d array
       pixellist, qind
       Output:
       dict{ 0-q0: values, 1-q1: values }
    '''
    import numpy as np
    label_num = np.unique(qind)
    intensity_distribution = {}
    for n in label_num:
        value = np.ravel(image)[ pixellist[qind == n] ]        
        intensity_distribution[n] = value
    return intensity_distribution



def get_pixellist_intensity( image, center, q=10, qwidth=2, qmask=None,):

    '''Get an intensity array of a q of qwidth,
       Giving: image, a two-d array
               center, of the image
               q: the interested q
               qwidth: q-width
               qmaks: the mask of the image
        
       Output:
       pixel_list, intensity
    '''
    import numpy as np
    qlist,qradi = make_qlist(qstart=q,qend=q+500,
                             qwidth= qwidth, noqs=1)     
    dimy,dimx = image.shape
    x,y=center
    pixellist,qind,nopr,nopixels = calqlist(qlist,qradi, dimx,dimy, x,y, qmask)
    intensity_dist = intensity_distribution(image, pixellist, qind)
    #L = len( intensity_dist.values()[0] )
    indices = np.arange( len(pixellist) )
    print ('average of intensity is:  %s'%(intensity_dist.values()[0]).mean())
    return indices, intensity_dist.values()[0]






def find_center2(image, est_center, qs=10, qe=40,
                    qw=1,  var=5, qmask=None,):

    import matplotlib.pyplot as plt
    from img_process import plot
    import numpy as np
    
    """
    This will find the most suitable center for the speckle pattern using an
    estimated center.
    Parameters
    ----------
     image : array
        image data dimensions are: (rr, cc)
    center : tuple
        estimated center
    inner_radius : float, optional
        inner radius of the inner-most ring
    width : float, optional
        ring thickness
    mask : array, optional
        boolen array shape image shape
    var : float, optional
        varying the selected center value to find the suitable center
        
    Returns
    -------
    center : tuple
       center for the speckle pattern
    """
 
    #print qradi
    m_value = float('inf')
     
    dimy,dimx = image.shape

    Imax = 0
    YM = 0
    for x in np.arange(est_center[0]-var, est_center[0]+var,.5):
        for y in np.arange(est_center[1]-var, est_center[1]+var,.5):
            cenx,ceny = x,y
            q,iq=azimuthal_interagate( image, qs, qe, qw,x,y,qmask)
            fx=q
            fy = iq*q**2
            ym = fy.min()
            yM = fy.max() 
             
            #print '@(x: %s,y: %s) the width is:  %s'%(x,y,wd)
            if yM > YM:
                center = [x,y]
                QIQ = fy  
    #plt.plot( Int, marker='o', c='r',ls='--' )
    #plt.show()
    plot( q, QIQ)
    return center


def find_center(image, est_center, inner_radius=10,
                    qwidth=1,  var=5, qmask=None,):

    import matplotlib.pyplot as plt
    import numpy as np
    from numpy.linalg import lstsq
    
    """
    This will find the most suitable center for the speckle pattern using an
    estimated center.
    Parameters
    ----------
     image : array
        image data dimensions are: (rr, cc)
    center : tuple
        estimated center
    inner_radius : float, optional
        inner radius of the inner-most ring
    width : float, optional
        ring thickness
    mask : array, optional
        boolen array shape image shape
    var : float, optional
        varying the selected center value to find the suitable center
        
    Returns
    -------
    center : tuple
       center for the speckle pattern
    """
    qlist,qradi = make_qlist(qstart=inner_radius,qend=inner_radius+500,
                             qwidth= qwidth, noqs=1)
    #print qradi
    m_value = float('inf')
     
    dimy,dimx = image.shape

    Imax = 0
    for x in np.arange(est_center[0]-var, est_center[0]+var,.5):
        for y in np.arange(est_center[1]-var, est_center[1]+var,.5):
            #print x,y
            pixellist,qind,nopr,nopixels = calqlist(qlist,qradi, dimx,dimy, x,y, qmask)
            intensity_dist = intensity_distribution(image, pixellist, qind)
            L = len( intensity_dist.values()[0] )
            indices = xrange( L )
            a = np.vstack([indices, np.ones(len(indices))]).T
            m, c = lstsq(a, intensity_dist.values()[0])[0]
            Im = (intensity_dist.values()[0]).mean()
            #print '@(x: %s,y: %s) the average intensity is:  %s'%(x,y,Im)
            if Im > Imax:
                Imax=Im
                X=x;
                Y=y;
            
            #print m, x, y,c #, a
            if abs(m) < m_value:
                m_value = abs(m) #abs(m)
                center = (x, y)
                Int = intensity_dist.values()[0]
    #print X,Y, Imax
    plt.plot( Int, marker='o', c='r',ls='--' )
    plt.show()
   
    return center






