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
############################################################################

'''
Functions:
1) plot(x,y, ax=None,yerr=None, title='img_',xlabel=None,ylabel=None,
         logx=False,logy=False,vlines=None,hlines=None,
         vmin=0,vmax=1, hmin=0,hmax=1, 
         save=False,show=True,outDir=None, xlim=None,ylim=None)
2) plot_with_fit(x,y, y_fit,yerr=None, title='img_',xlabel=None,ylabel=None,logx=False,logy=False,
         save=False,show=True,outDir=None, xlim=None,ylim=None)

3) def show_fit(data, para, fx = None, title = None, xlim=None,
             ax = None, save = False,
             outDir = None, fit_func = 'Gauss')
             
4) pd_data_plot( pdata, x=None, y=None, title='data',
            xlabel=None,ylabel=None,
         logx=False,logy=False,
         save=False,show=True,outDir=None, xlim=None,ylim=None)


5) show_img( img, ax=None, save = False, vmin = None,
             vmax = None, cmap = 'spectral', fontsize = 24,
             axis_on = True, title_on = True, xlabel = None,
             ylabel = None, aspect = None, title = 'img_',
             show = True, set_cb = False, close= False, 
             outDir = None, logs = False, sizex = 9, sizey = 9, ylim = None,
             xlim = None, xticks = True, yticks = True, extent = None)


6) show_mask( img, mask, only_mask=False,
               vmin=None,vmax=None, xlim=None,ylim=None,title='img',                      
            logs=True,set_cb=True,save = False,outDir=None  )


7) show_img_with_ROI(img,pixellist, qind,vmin=None,vmax=None,
                     xlim=None,ylim=None,title='img',                      
            logs=True,set_cb=True,save = False,outDir=None  )



'''







#####for data plot
###########
def plot(x,y, ax=None,yerr=None, title='img_',xlabel=None,ylabel=None,
         logx=False,logy=False,vlines=None,hlines=None,
         vmin=0,vmax=1, hmin=0,hmax=1, 
         save=False,show=True,outDir=None, xlim=None,ylim=None):
         
    '''
    plot 1D data
    '''
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    if ax is None:fig,ax = plt.subplots()    
    ax.plot( x, y, 'bo', ls='--', linewidth=3, markersize=8 );
    if xlim   is not  None:plt.xlim(xlim)
    if ylim   is not  None:plt.ylim(ylim)
    if title!=None:plt.title(  title ,fontsize=28)
    if xlabel!=None:plt.xlabel(xlabel,fontsize=24)
    if ylabel!=None:plt.ylabel(ylabel,fontsize=24)
    if logy:ax.set_yscale('log')
    if logx:ax.set_xscale('log')

    if vlines is not None:
        for vl in vlines:
            ax.vlines( vl,hmin,hmax )

    if hlines is not None:
        for hl in hlines:
            ax.hlines( hl,vmin,vmax )
    
    if save:
        if outDir!=None:fp=outDir + title + '_.png'
        else:fp= title + '_.png'
        plt.savefig( fp )        
    if show:plt.show()    




############################################################################
###show a data with a fitted curve

def plot_with_fit(x,y, y_fit,yerr=None, title='img_',xlabel=None,ylabel=None,logx=False,logy=False,
         save=False,show=True,outDir=None, xlim=None,ylim=None):

    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    fig,ax = plt.subplots()    
    if yerr is None:ax.plot( x, y, 'bo',  markersize=8 );
    else:ax.errorbar(x,y,yerr,marker='o',c='b',  markersize=8 );
    ax.plot( x,y_fit, c='red', ls='--', linewidth=3,)
    
    if xlim   is not  None:plt.xlim(xlim)
    if ylim   is not  None:plt.ylim(ylim)
    if title!=None:plt.title(  title ,fontsize=36)
    if xlabel!=None:plt.xlabel(xlabel,fontsize=24)
    if ylabel!=None:plt.ylabel(ylabel,fontsize=24)
    if logy:ax.set_yscale('log')
    if logx:ax.set_xscale('log')
    
    if save:
        if outDir!=None:fp=outDir + title + '_.png'
        else:fp= title + '_.png'
        plt.savefig( fp )        
    if show:plt.show() 


############################################################################
###show a data with fitting obtained by a fit_func model

def show_fit(data, para, fx = None, title = None, xlim=None,
             ax = None, save = False,
             outDir = None, fit_func = 'Gauss'):
    '''
    show a data with fitting obtained by a fit_func model
    data: a two column array, e.g. [x,y]
    para: the obtained fitting parameters by the fit_func
    ax: plot
    fit_func: supported func includes Gauss, ploynominal
    

    '''
    
    import matplotlib.pyplot as plt
    import numpy as np
    (x, y,) = data
    x = np.array(x)
    y = np.array(y)
    if fx is None:
        fx = x
    else:
        fx = np.array(fx)
    x0 = np.linspace(fx.min(), fx.max(), 2000)
    if fit_func == 'polyn':
        func = np.polyval
        fit = func(para, x0)
    elif fit_func == 'Gauss':
        func = gauss
        fit = func(x0, *para)
    else:
        fit = func(x0, *para)
    if ax is None:
        (fig, ax,) = plt.subplots()
    if title is not None:
        ax.set_title(title)
    if xlim is not None:ax.set_xlim( xlim )    
    ax.plot(x, y, 'bo')
    ax.plot(x0, fit, 'r', ls='-')
    
    if save:
        if outDir != None:
            fp = outDir + title + '_.png'
        else:
            fp = title + '_.png'
        plt.savefig(fp)





############################################################################
###plot a panda.framedata
def pd_data_plot( pdata, x=None, y=None, title='data',
            xlabel=None,ylabel=None,
         logx=False,logy=False,
         save=False,show=True,outDir=None, xlim=None,ylim=None):
    '''plot a panda data'''
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    if x is None:x= pdata.index
    if y is None:y= pdata.keys()
    N = len( pdata )   
    fig,ax=plt.subplots()
    pdata.plot(x=x,y=y,marker='o',ls='--',logx=logx,logy=logy, ax=ax
               );

    if xlim   is not  None:plt.xlim(xlim)
    if ylim   is not  None:plt.ylim(ylim)    
    if title!=None:plt.title(  title ,fontsize=24)
    if xlabel!=None:plt.xlabel(xlabel,fontsize=24)
    if ylabel!=None:plt.ylabel(ylabel,fontsize=24)
      
    if save:
        if outDir!=None:fp=outDir + title + '_.png'
        else:fp= title + '_.png'
        plt.savefig( fp )        
    if show:plt.show() 



##########################################################
###for images

 
def determine_aspect(shape, extent):
    dx = (extent[1] - extent[0]) / float(shape[1])
    dy = (extent[3] - extent[2]) / float(shape[0])
    return dx / dy

#'YlGnBu'
def show_img( img, ax=None, save = False, vmin = None,
             vmax = None, cmap = 'spectral', fontsize = 24,
             axis_on = True, title_on = True, xlabel = None,
             ylabel = None, aspect = None, title = 'img_',
             show = True, set_cb = False, close= False, 
             outDir = None, logs = False, sizex = 9, sizey = 9, ylim = None,
             xlim = None, xticks = True, yticks = True, extent = None):
    """show a two-D image"""
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np
    from matplotlib.colors import LogNorm
    
    if ax is None:
        fig, ax = plt.subplots()
        ax = fig.gca()    
    fig = ax.figure    
    if vmin == None:
        vmin = img.min()
    if vmax == None:
        vmax = img.max()
    shape = img.shape
    (dy, dx,) = shape
    if xlim is not None:
        ax.set_xlim(xlim)
    if extent is not None:
        (x1, x2, y1, y2,) = extent
        if ylim is None:
            ylim = [y2, y1]
        aspect_ = determine_aspect(shape, extent)
    else:
        aspect_ = None
    if ylim is None:
        ylim = [0, dy]
    ax.set_ylim([ylim[0], ylim[1]])
    if not logs:
        cax = ax.imshow(img, cmap=cmap, vmin=vmin,
                        vmax=vmax, aspect=aspect_, extent=extent)
    if logs:
         
        if vmin is None:
            vmin = img[nonzero(img)].min()
        if vmin ==0.0:vmin+=1e-1   
        cax = ax.imshow(img, cmap=cmap,
                    norm=LogNorm(vmin=vmin, vmax=vmax),
                    aspect=aspect_, extent=extent)

    if aspect is not None:
        im = ax.get_images()
        (x1, x2, y1, y2,) = im[0].get_extent()
        if ylim is not None:
            (y1, y2,) = ylim
        if xlim is not None:
            x1.x2 = xlim
        ax.set_aspect(abs((x2 - x1) / (y2 - y1)) / aspect)
    if title_on:
        plt.title(title, fontsize=fontsize)
    if not axis_on:
        plt.axis('off')
    if xlabel is not None:
        plt.xlabel(xlabel, fontsize=fontsize)
    else:
        plt.xlabel('')
    if ylabel is not None:
        plt.ylabel(ylabel, fontsize=fontsize)
    else:
        plt.ylabel('')
    if xticks is None:
        plt.xticks([])
    elif xticks is True:
        plt.setp(ax.get_xticklabels(), fontsize=fontsize * 0.8, visible=True)
    else:
        plt.xticks(xticks, fontsize=fontsize * 0.8)
    if yticks is None:
        plt.yticks([])
    elif yticks is True:
        plt.setp(ax.get_yticklabels(), fontsize=fontsize * 0.8, visible=True)
    else:
        plt.yticks(yticks, fontsize=fontsize * 0.8)
    if set_cb:
        cbar = fig.colorbar(cax, ticks=[vmin, vmax])
        #print 'here'
    #print 'jj'
    #cbar = fig.colorbar(cax, ticks=[vmin, vmax])
    if save:
        if outDir != None:
            fp = outDir + title + '_.png'
        else:
            fp = title + '_.png'
        plt.savefig(fp)
        
    #if close:
        #print 'close'
        #plt.close()
        #plt.close('all')
    if show:
        plt.show()


def show_mask( img, mask, only_mask=False,
               vmin=None,vmax=None, xlim=None,ylim=None,title='img',                      
            logs=True,set_cb=True,save = False,outDir=None  ):
    '''Show image with mask'''
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np
    from matplotlib.colors import LogNorm


    dimy,dimx= img.shape
    #imgx= np.zeros( [dimy,dimx])
    img_ = img.copy()
    if only_mask:img_[:]=0
    img_[mask]=  0
    show_img( img_, logs=logs,vmin=vmin,vmax=vmax,
              set_cb=set_cb,xlim=xlim,ylim=ylim,
              save = save,outDir = outDir,title=title )
    
 
    

def show_img_with_ROI(img,pixellist, qind,mask=None,vmin=None,vmax=None,
                     xlim=None,ylim=None,title='img',                      
                       logs=True,set_cb=True,save = False,outDir=None  ):
    """give the img, pixellist, qind to show the edf file
        with the interested area defined by pixellist and qind """
    import numpy as np
    
    dimy,dimx= img.shape
    #imgx= np.zeros( [dimy,dimx])
    imgx = img.copy()
    if mask is not None:imgx[mask]=  0
    img_ = imgx.ravel()        
    #img_ = np.array( img_, dtype = np.int32)
    #img.dtype='int32'
    img_[pixellist]= qind
    img_=img_.reshape( [ dimy,dimx] )
    show_img( img_, logs=logs,vmin=vmin,vmax=vmax,
              set_cb=set_cb,xlim=xlim,ylim=ylim,
              save = save,outDir = outDir,title=title )




#def show_ROI(img,pixellist, qind,mask=None,vmin=None,vmax=None,
#                     xlim=None,ylim=None,title='img',                      
#                       logs=True,set_cb=True,save = False,outDir=None  ):
#    """give the img, pixellist, qind to show the edf file
#        with the interested area defined by pixellist and qind """
#    import numpy as np
#    
##    dimy,dimx= img.shape
#   #imgx= np.zeros( [dimy,dimx])
#    imgx = img.copy()
##    
#    if mask is not None:imgx[mask]=  0
#    img_ = imgx.ravel()        
#    #img_ = np.array( img_, dtype = np.int32)
#    #img.dtype='int32'
#    img_[pixellist]= qind
#    img_=img_.reshape( [ dimy,dimx] )
#    show_img( img_, logs=logs,vmin=vmin,vmax=vmax,
#              set_cb=set_cb,xlim=xlim,ylim=ylim,
#              save = save,outDir = outDir,title=title )





