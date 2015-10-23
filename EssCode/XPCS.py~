# -*- coding: utf-8 -*-
######################################################################################
##Revised Based on Yorick_Multitau_Code Obtained from Dr. Andrei Fluerasu###############
#####################Do one-time correlation ####################################
####################   Coded by Dr. Yugang Zhang #################################
####################   631-885-4714   ############################################
####################   Email: yuzhang@bnl.gov     ####################################
########################################################################################
#################### or: yuzhangnew@icloud.com , Octa 18, 2015##     ##################################
#######################################################################################
######################################################################################
################################################################################
##############################################################################

'''
Functions
1) __init
    from XPCS_Functions import readframe_series,make_qlist,calqlist
    self.data = readframe_series
    self.qlist,self.qradi = make_qlist
    self.pixellist,self.qind,self.nopr,self.nopixels = calqlist

2) delays(self,time=1,
            nolevs=None,nobufs=None, tmaxs=None)
            
3)process(self,lev,bufno, n=None)

4) insertimg(self, n, norm=None, print_=False, brute=False)
5) autocor( self, print_=False,  brute=False,)
6) g2_to_pds(self, dly, g2, tscale = None)
7) show(self,g2p, qnum, title=None)
8) showg2(self,data='g2',show=True,show_fit=False,
               fit_para=None, save_=True,filename='g2' )

               
One example use:

xp = xpcs()
dly = xp.delays( time= 1)
g2=xp.autocor(  )
xp.showg2(g2)


'''


from numpy import pi,sin,arctan,sqrt,mgrid,where,shape,exp,linspace,std,arange
from numpy import power,log,log10,array,zeros,ones,reshape,mean,histogram,round,int_
from numpy import indices,hypot,digitize,ma,histogramdd,apply_over_axes,sum
from numpy import around,intersect1d, ravel, unique,hstack,vstack,zeros_like
from numpy import save, load, dot,isnan
from numpy.linalg import lstsq
from numpy import polyfit,poly1d;
import sys
import pandas as pd
import matplotlib.pyplot as plt 
import time

from Setup_file import  *

class xpcs( object):
    def __init__(self,):
        """ DOCUMENT __init__(   )
        the initilization of the XPCS class
          """
        from Get_Data import readframe_series
        from XPCS_Functions import make_qlist,calqlist
        #from Setup_file import qstart,qend,qwidth,noqs,sid, inDir 
	#from Setup_file import dimx,dimy,cenx,ceny, mask

        self.version='version_1'
        self.create_time='Octa_18_2015'
        self.author='Yugang_Zhang_chx11id_nsls2_BNL'
        self.email='yuzhang@bnl.gov'
        self.email2='yuzhangnew@icloud.com'
        
        self.data = readframe_series( sid, inDir   ) #read data with sid, and bg_id
        self.qlist,self.qradi = make_qlist(qstart,qend,qwidth,noqs)
        self.pixellist,self.qind,self.nopr,self.nopixels = calqlist(
            self.qlist,self.qradi, dimx,dimy,
                        cenx,ceny,qmask=mask)
        if min(self.qind)!=0:print ('The qmin is too small!')
        
    def delays(self,time=1,
               nolevs=None,nobufs=None, tmaxs=None): 
        ''' DOCUMENT delays(time=)
        Using the lev,buf concept, to generate array of time delays
        return array of delays.
        KEYWORD:  time: scale delays by time ( should be time between frames)
                  nolevs: lev (a integer number)
                  nobufs: buf (a integer number)
                  tmax:   the max time in the calculation, usually, the noframes
                  
        '''
         
        if nolevs is None:nolevs=nolev #defined by the set-up file
        if nobufs is None:nobufs=nobuf #defined by the set-up file
        if tmaxs is None:tmaxs=tmax    #defined by the set-up file
        if nobufs%2!=0:print ("nobuf must be even!!!"    )
        dly=zeros( (nolevs+1)*nobufs/2 +1  )        
        dict_dly ={}
        for i in range( 1,nolevs+1):
            if i==1:imin= 1
            else:imin=nobufs/2+1
            ptr=(i-1)*nobufs/2+ arange(imin,nobufs+1)
            dly[ptr]= arange( imin, nobufs+1) *2**(i-1)            
            dict_dly[i] = dly[ptr-1]            
        dly*=time 
        #dly = dly[:-1]
        dly_ = dly[: where( dly < tmaxs)[0][-1] +1 ]
        self.dly=dly 
        self.dly_=dly_
        self.dict_dly = dict_dly
        #print nolevs,nobufs
        return dly
 


    ###########################################################################
    ########for one_time correlation function for xyt frames
    ##################################################################
    def process(self,lev,bufno, n=None):
        
        num[lev]+=1  
        if lev==0:imin=0
        else:imin=nobuf/2        
        for i in range(imin, min(num[lev],nobuf) ):
            ptr=lev*nobuf/2+i    
            delayno=(bufno-i)%nobuf 
            IP=buf[lev,delayno]
            IF=buf[lev,bufno]
            G[ptr]+= (  (histogram(qind, bins=noqs, weights= IF*IP))[0]/nopr-G[ptr] )/ (num[lev]-i)
            IAP[ptr]+= (  (histogram(qind, bins=noqs, weights= IP))[0]/nopr-IAP[ptr] )/ (num[lev]-i)
            IAF[ptr]+= (  (histogram(qind, bins=noqs, weights= IF))[0]/nopr-IAF[ptr] )/ (num[lev]-i)

    def insertimg(self, n, norm=None, print_=False, brute=False):
        
        cur[0]=1+cur[0]%nobuf        
        #img =  readframe_series(n )
        img =  self.data[n]
        if print_:print ('The insert image  is %s' %(n))
        buf[0, cur[0]-1 ]= array( img.flatten()[pixellist] , dtype=float)       
        img=[] #//save space    
        self.process(lev=0, bufno=cur[0]-1, n=n )    
        processing=1
        if not brute:
            lev=1
            while processing: 
                if cts[lev]:
                    prev=  1+ (cur[lev-1]-1-1+nobuf)%nobuf
                    cur[lev]=  1+ cur[lev]%nobuf
                    buf[lev,cur[lev]-1] = ( buf[lev-1,prev-1] + buf[lev-1,cur[lev-1]-1] ) /2.
                    cts[lev]=0                 
                    self.process(lev= lev, bufno= cur[lev]-1 , n=n)        
                    lev+=1
                    if lev<nolev:processing = 1
                    else:processing = 0                                
                else:
                    cts[lev]=1      #// set flag to process next time
                    processing=0    #// can stop until more images are accumulated

    def autocor( self, print_=False,  brute=False,):
        global buf,G,IAP,IAF,num,cts,cur,g2,gmax,sigg2
        global Ndel,Npix,pixellist,qind,nopr,nopixels
        
        (pixellist,qind,nopr,nopixels) = (self.pixellist,
                                    self.qind,self.nopr,self.nopixels)
        start_time = time.time()
        #initialize all arrays
        buf=zeros([nolev,nobuf,nopixels])  #// matrix of buffers
        cts=zeros(nolev)
        cur=ones(nolev) * nobuf        
        G=zeros( [(nolev+1)*nobuf/2,noqs])
        IAP=zeros( [(nolev+1)*nobuf/2,noqs])
        IAF=zeros( [(nolev+1)*nobuf/2,noqs])
        num= array(zeros(  nolev ),dtype='int')        
        ttx=0
        print ('Doing g2 caculation of %s frames for file--%s'%(noframes,sid))
        for n in range(1,noframes +1 ):
            #print n
            self.insertimg( begframe+n-1, print_=print_,brute=brute)            
            if  n %(noframes/10) ==0:
                sys.stdout.write("#")
                sys.stdout.flush()
        elapsed_time = time.time() - start_time
        print ('Total time: %.2f min' %(elapsed_time/60.))
        #print G.shape    
        if (len(where(IAP==0)[0])!=0) and ( 0 not in self.nopr):
            gmax = where(IAP==0)[0][0]        
        else:
            gmax=IAP.shape[0]
        #g2=G/(IAP*IAF)
        #print G
        g2=(G[:gmax]/(IAP[:gmax]*IAF[:gmax]))       
        self.g2=g2   
        return g2  #,elapsed_time/60.


    def g2_to_pds(self, dly, g2, tscale = None):
        '''convert g2 to a pandas frame'''        
        if len(g2.shape)==1:g2=g2.reshape( [len(g2),1] )
        tn, qn = g2.shape
        tindex=xrange( tn )
        qcolumns = ['t'] + [ 'g2_q%s'%q for q in range(1,1+qn)  ]
        if tscale is None:tscale = 1.0
        #print tn,qn,tindex
        g2t = hstack( [dly[:tn].reshape(tn,1) * tscale, g2 ])
        #print g2t.shape
        g2p = pd.DataFrame(data=g2t, index=tindex,columns=qcolumns)         
        return g2p
    

    def show(self,g2p, qnum, title=None):
        t = g2p.t  
        N = len( g2p )
        ylim = [g2p['g2_q%s'%qnum].min(), g2p[1:N]['g2_q%s'%qnum].max()  ]
        
        g2p.plot(x=t,y='g2_q%s'%qnum,marker='o',ls='--',logx=True,ylim=ylim);
        plt.xlabel('time delay, s',fontsize=12)
        if title is None:title='Q_%s'%(qnum)
        plt.title(title)
        #plt.savefig( RES_DIR + title +'.png' )       
        plt.show()


    def showg2(self,data='g2',show=True,show_fit=False,
               fit_para=None,
               save_=True,filename='g2' ):
        qradi = self.qradi

        sx= int( round (sqrt(noqs)) )
        if noqs%sx==0:sy=noqs/sx
        else:sy=noqs/sx+1
        
        if data=='g2':data=self.g2 
        m,n=data.shape        
        fig = plt.figure()
        fig.set_size_inches(20,10.5)
        title=filename
        plt.title(title,fontsize=28, y =1.08)
        
        plt.axis('off')
        #fig.tight_layout()
        for sn in range(0, n):
            #ax = fig.add_subplot('%i'%sx,'%i'%sy,'%i'%(sn+1) )
            ax = fig.add_subplot(sx,sy,sn+1 )

##            plt.annotate('q= %.5f A-1'%(qradi[sn]*qperpixel),
##              xy=(0, 1), xytext=(12, -12), va='top',
##             xycoords='axes fraction', textcoords='offset points')
             
            txt = r'$q= %.5f$ $\AA^{-1}$'%(qradi[sn]*qperpixel)            
##            ax.text(.3,1.02, txt,
##            fontsize=20,  horizontalalignment='left',
##                transform=ax.transAxes)            
            #plt.title('q= %.3f'%qradi[sn],fontsize=16)
            
            plt.title(txt,fontsize=20, y =1.02)
            x=self.dly[:m]
            #x=self.dly[1:m+1]
            
            y=data[:,sn]
            if not show_fit:
                plt.plot(x, y, 'o',linewidth=3, ls='--', color='b',markersize=8)            
            #ax.set_xscale('log')
            if show_fit:
                plt.plot(x, y, 'o',linewidth=1, ls='', color='b',markersize=8) 
                yf = g2_tau(x, para[sn]) #x is time delay, sn = the qs
                plt.plot(x, yf, linewidth=3, ls='-', color='r')  
                
            plt.xscale('log')
            #print max(y) 
            plt.ylim([    min(y) , max(y[1:]) ])
            plt.xlim([   min(x)+0*1e-6, max(x)])
            #plt.ylim([  1, 2])
            plt.xlabel(r'$\tau,$ $s$',fontsize=20)
            #ax.xaxis.set_label_coords(0.5, -0.025)
            
            plt.ylabel(r'$g_2$',fontsize=20, )
            #plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
            #plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
            yfmt = plt.ScalarFormatter()
            yfmt.set_useOffset(0)
            ax.yaxis.set_major_formatter(yfmt)            
            #fig.tight_layout()
            #fig.subplots_adjust(hspace=2)
        fig.tight_layout()    
        if save_:
            plt.savefig( outDir + filename +'.png' )            
            #cpdump(data,FOUT+'_g2',RES_DIR)
            #cpdump(sigg2,FOUT+'_sigg2',RES_DIR)
            
        if show:plt.show()













