# -*- coding: utf-8 -*-
######################################################################################
#######################################################################################
#####################Do Visiblity Analysis ####################################
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
Functions:
1) __init
    from XPCS_Functions import make_qlist,calqlist,get_qRings
    self.qlist,self.qradi = make_qlist
    self.pixellist,self.qind,self.nopr,self.nopixels = calqlist
    self.qRings =  get_qRings
2) process(self,lev,bufno, buf, num, qRings, BinEdges,
                Histrogram2_3D_Matrix_samplei,Histrogram_3D_Matrix_samplei)
4) _get_trace_and_maxc(self, sid,fs=None,nf=None ):
    from XPCS_Functions import readframe_series,get_trace_and_maxc
5) XSVS_Single( self, sid,  expose_time = 1.0, maxcts = 50,
            do_bining = True,timebin = 2,timebin_level=None)
    from XPCS_Functions import readframe_series  

6) fit_K_PK_witherr(self,func, initParams, xdata,ydata,yerr,roi,initKValue=None)
7) fit_K_PK(self,func, initParams, xdata,ydata,roi,initKValue=None)
8) FitVisibility(self, bins, his, dev=None, trace_mean=None,
                      fix_K = False, fit_gm=True)

9) model_to_pds( self, model, qv=None)
10) 

one example:

xs = xsvs()
trace, maxcts = xs._get_trace_and_maxc(sid,fs,nf )
bins,his, dev = xs.XSVS_Single(sid,
            expose_time = 1.0, maxcts = maxcts,
            do_bining = True,timebin = 2,)
trace_mean = trace.mean()    
model=xs.FitVisibility(bins, his, dev, trace_mean,
                    fix_K = False,fit_gm=True)

    
'''



from math import isnan
import numpy as np
from Setup_file import *
from XSVS_Functions import *

class xsvs( object):
    def __init__(self):
        """ DOCUMENT __init__(   )
        the initilization of the XSVS class
          """
        from XPCS_Functions import make_qlist,calqlist, get_qRings
        
        self.version='version_1'
        self.create_time='Octa_18_2015'
        self.author='Yugang_Zhang@chx11id_nsls2_BNL'
        self.email='yuzhang@bnl.gov'
        self.email2='yuzhangnew@icloud.com'
        #self.data = readframe_series( sid, inDir   ) #read data with sid, and bg_id

        self.qlist,self.qradi = make_qlist(qstart,qend,qwidth,noqs)
        self.pixellist,self.qind,self.nopr,self.nopixels = calqlist(
            self.qlist,self.qradi, dimx,dimy,
                        cenx,ceny,qmask=mask)
        self.qRings =  get_qRings(  self.qind, self.pixellist)
        

    def process(self,lev,bufno, buf, num, qRings, BinEdges,
                Histrogram2_3D_Matrix_samplei,Histrogram_3D_Matrix_samplei):
        
        num[lev]+=1
        #print num
        #print num[lev]-1
        for j in range(noqs):
            #change noqs to qbins
            #qRings_j = pixellist[where(qind==j)]
            data = buf[lev,bufno][ qRings[j] ]  #/float(nopr[j])                              
            qHist,qBins = np.histogram(data,bins=BinEdges,normed=True)#normalized by pixel number
            
            #qHist,qBins = np.histogram(data,bins=BinEdges,normed=True)#normalized by pixel number
            
            #qHist,qBins = histogram(data,bins=BinEdges,normed=False
            #print BinEdges
            #print data
            #print qHist
            #print qBins
            #qHist *= 2**(lev+1)
            if isnan( qHist[0] ):
                #print ('Here is the nan!!!')
                qHist = np.zeros_like( qHist)
                    #here assume histBins same for all the qs
            Histrogram2_3D_Matrix_samplei[lev,j] += (
                np.power(qHist,2) - Histrogram2_3D_Matrix_samplei[lev,j] )/(num[lev])                
            Histrogram_3D_Matrix_samplei[lev,j] += (
                qHist - Histrogram_3D_Matrix_samplei[lev,j])/(num[lev])
            #trace[si,lev,j,num[lev]-1] = mean(data)
            #print lev,j,num[lev]


 
 
    def _get_trace_and_maxc(self, sid,fs=None,nf=None ):
        
        from Get_Data import readframe_series
        from XPCS_Functions import get_trace_and_maxc
        imgs = readframe_series( sid, inDir )
        print ('Doing trace analysis for sid_%s with %s frames'%(sid,noframes ))
        trace, max_cts = get_trace_and_maxc(imgs,self.pixellist,self.qind,fs,nf)        
        return trace,max_cts
        
       
    def XSVS_Single( self, sid,  expose_time = 1.0, maxcts = 50,
            do_bining = True,timebin = 2,timebin_level=None):

        '''sid:   give the sample name/folder 
           qRings:  a pixel list for different qs
           maxcts: the max counts for the lowest exposure time
           bins:
              if True, bin the images; bin number is timebines,
                       timebin_level gives the max level for bin
              if Flase: just do histgram for the images

          Return: expose time: a list with length as timebin_level 
                  Bins: for each expose time, each Qs, shaped as [time_lev, Qs]
                  Histrogram_3D_Matrix: a histrogram, shaped as [time_lev, Qs]
                  HistStddev_all: a deviation of histrogram, shaped as [time_lev, Qs]
           
        '''
        from Get_Data import readframe_series        
        img = readframe_series( sid, inDir )
        N = len(img)
        
	
        if noframes>N:noframes_=N
        else:noframes_=noframes
        
        qRings = self.qRings        
        if do_bining:        
            if timebin_level is None:
                timebin_level = int( np.log( noframes  )/np.log(2)) #  +1            
            time_list = [ 2**i for i in range(timebin_level)]
        else:
            timebin_level = 1
            time_list = [  expose_time ]
            
        Bins =  np.zeros([timebin_level,noqs],dtype=object)        #GlobalPDFArray_all  is the Histrogram_3D_Matrix
        Histrogram_3D_Matrix = np.zeros([timebin_level,noqs],dtype=object)

        Histrogram2_3D_Matrix = np.zeros([timebin_level,noqs],dtype=object)
        #this function is to calculate the deviation of histrogram        

        BinEdges = np.arange( maxcts )
        #print BinEdges

	#print img      
        for i,n in enumerate(time_list):
            for j in range(noqs): 
                Bins[ i,j ]=  np.arange( maxcts * time_list[i] )
        ###########################################################
     
        print ('Running sample---%s'%(sid))
        print ('Doing histrogram for %s frames'%noframes_         )
        buf=np.zeros([ timebin_level, timebin], dtype=object)  #// matrix of buffers
        cts=np.zeros(timebin_level)
        cur=np.ones(timebin_level) * timebin
        num= np.array(np.zeros(  timebin_level ),dtype='int')

        ###############################################
        #the loop for process of each frame
        #print noframes_,noframes
        for n in range(1,noframes_ +1 ):
            #print n+ begframe-1,
            cur[0]=1+cur[0]%timebin                       
            ifg =  img[ n + begframe-1 ].ravel()           
            ifg[ (np.where(ifg<0))[0] ]=0
            buf[0, cur[0]-1 ]= np.array( ifg, dtype=float)    
            #img=[] #//save space                
            #self.process(lev=0, bufno=cur[0]-1 )                
            self.process(0,cur[0]-1, buf, num, qRings, Bins[ 0,0 ],
                        Histrogram2_3D_Matrix,Histrogram_3D_Matrix)
            
            processing=1        
            lev=1                
            #print cts
            if do_bining==False:processing = 0
            while processing: 
                if cts[lev]:
                    prev=  1+ (cur[lev-1]-1-1+nobuf)%timebin
                    cur[lev]=  1+ cur[lev]%timebin
                    buf[lev,cur[lev]-1] = (
                        buf[lev-1,prev-1] + buf[lev-1,cur[lev-1]-1] ) #/2.
                    cts[lev]=0                 
                    #self.process(lev= lev, bufno= cur[lev]-1 , n=n)
                    self.process(lev,cur[lev]-1, buf, num, qRings, Bins[ lev,0 ],
                                Histrogram2_3D_Matrix,Histrogram_3D_Matrix)                        
                    lev+=1
                    if lev<timebin_level:processing = 1
                    else:processing = 0                                
                else:
                    
                    cts[lev]=1      #// set flag to process next time
                    processing=0    #// can stop until more images are accumulated

            if isnan( Histrogram_3D_Matrix[0,0][0] ):
                print (n)
 
        HistStddev_all = np.power( ( Histrogram2_3D_Matrix  -
                            np.power(Histrogram_3D_Matrix ,2)), .5)  
        #expost = time_list
        self.time_list = np.array( time_list ) *expose_time
        self.timebin_level=timebin_level
        self.noframes = noframes
        return  Bins, Histrogram_3D_Matrix, HistStddev_all


    def fit_K_PK_witherr(self,func, initParams, xdata,ydata,yerr,roi,initKValue=None):
        #logErr = yerr/(ydata*log(10))
        #print   yerr       #logErr[roi] #ydata[roi]        
        from scipy.optimize import leastsq
        
        if initKValue==None:#for not fix K
            plsq = leastsq(func,initParams,\
                    args=( ydata[roi],xdata[roi],yerr[roi] ),full_output=1)
                #args=(log10(ydata[roi]),xdata[roi],logErr[roi],),full_output=1)
        else:#for fix K
            plsq = leastsq(func,initParams,\
                    args=(ydata[roi],xdata[roi],yerr[roi],initKValue),full_output=1)
                #args=(log10(ydata[roi]),xdata[roi],logErr[roi],initKValue),\full_output=1)
            
        s_sq = (plsq[2]['fvec']**2).sum()/\
                    (len(ydata[roi])-len(initParams))
        if hasattr(plsq[1], 'diagonal'):paramStdDev = sqrt(plsq[1].diagonal()*s_sq) 
        else:
            paramStdDev = 0,0
            #print ('here!!!')
        return plsq[0],paramStdDev
    
    def fit_K_PK(self,func, initParams, xdata,ydata,roi,initKValue=None):
        #logErr = yerr/(ydata*log(10))
        #print   yerr       #logErr[roi] #ydata[roi]
        from scipy.optimize import leastsq
        
        if initKValue==None:#for not fix K
            plsq = leastsq(func,initParams,\
                    args=( ydata[roi],xdata[roi] ),full_output=1)
                #args=(log10(ydata[roi]),xdata[roi],logErr[roi],),full_output=1)
        else:#for fix K
            plsq = leastsq(func,initParams,\
                    args=(ydata[roi],xdata[roi],initKValue),full_output=1)
                #args=(log10(ydata[roi]),xdata[roi],logErr[roi],initKValue),\full_output=1)
            
        s_sq = (plsq[2]['fvec']**2).sum()/\
                    (len(ydata[roi])-len(initParams))
        if hasattr(plsq[1], 'diagonal'):paramStdDev = sqrt(plsq[1].diagonal()*s_sq) 
        else:
            paramStdDev = 0,0
            #print ('here!!!')
        return plsq[0],paramStdDev
    
 
    def FitVisibility(self, bins, his, dev=None, trace_mean=None,
                      fix_K = False, fit_gm=True):
        #global qradi,qperpixel
        '''Function calculating visibility for each of the exposures
        listed in the input file.
        fix_K: give a fix K parameter for the fit
        
        Global Varibles:
        pixellist, qind, noqs,qradi, dimx, dimy,
        
        RES_DIR , outPrefix, dataPref,
        
        exposureList: list for the exposure time
        createFastStatic_frame_start =1
        createFastStatic_frame_end =200
        '''
        from Process_Data import cpdump
        
        timebin_level,noqs = bins.shape
        beta = np.zeros((noqs,2))
        
        MArray = np.zeros((noqs,2))
        KArray = np.zeros((noqs,2))
        RArray = np.zeros((noqs,2))
        
        MArray_gm = np.zeros((noqs,2))
        KArray_gm = np.zeros((noqs,2))

        
        dataToSave = np.zeros([ timebin_level,noqs ], dtype=object)# prepare the data container
        #for each point (t,q) save a fitted parameter, K,M
        tra_mean = trace_mean #trace.mean()
        for i in range(timebin_level): #for each timelevel
            for j in range(noqs): #for each q
                #trace.shape = (samples, tlevel, noqs, frames)
                
                #initKValue = trace[:,i,j,:].mean() # inital guess for K               
                #initKValue = trace[:,i,j,:].mean() * 2**(i+1)
                if len(trace_mean.shape)==1:
                    initKValue = tra_mean['q%i'%j] * 2**(i+1)
                else:
                    initKValue =tra_mean['q%i'%j][i]                 
                # Fit a negative binomial distribution pmf to
                #the histograms for each q
                bins[i,j] = np.array( bins[i,j] )
                xdata = bins[i,j][:-1]
                ydata = his[i,j]
                yerr =  dev[i,j]
                logErr = yerr/(ydata*np.log(10))
                
                #if min(ydata)==0:
                    #ydata += 1e-15;
                    #print '#===Here add 1e-15 to the histrogram!===#'
                
                #logErr = yerr/(ydata*log(10))
                
                #x = linspace(xdata[0],xdata[-1],num=50)                
                M0,K0 = 5.0, initKValue                
                # Set a roi for the fitting range
                roi = np.where(ydata>1e-5)
                #print roi[0].shape
                if len(roi[0]) > len(ydata)-2:
                    roi = (np.array(roi[0][:-2]),)                    
                elif len(roi[0]) < 2:
                    roi = np.where(ydata>=0)
                #print roi[0].shape
                #print ydata,roi
                #roi = arange( len(ydata))
                ###########################
                ########## Do fit #########
                ###########################
                if not fix_K:
                    func=residuals;initParams=[initKValue,5.0];
                    #func=residuals_err;initParams=[initKValue,5.0];
                    iKV=None
                    if fit_gm:
                        #func_gamma = residuals_gamma_err;
                        func_gamma = residuals_gamma;
                    #print 'Not fix K'
                    
                else:
                    func=residuals2;initParams=[5.0];iKV=initKValue
                    #func=residuals2_err;initParams=[5.0];iKV=initKValue
                    if fit_gm:
                        #func_gamma = residuals2_gamma_err;
                        func_gamma = residuals2_gamma;
                
                plsq0,paramStdDev=self.fit_K_PK(func, initParams,
                                    xdata,ydata,roi,iKV )           
                                    #xdata,ydata,yerr,roi,iKV )            
                                    #xdata,log10(ydata),logErr,roi,iKV )
                                    #xdata,ydata,yerr,roi,iKV )
                #print i,j, plsq0
                if fit_gm:
                    plsq0_gm,paramStdDev_gm=self.fit_K_PK(func_gamma, initParams,
                                    xdata,ydata,roi,iKV )     
                                    #xdata,ydata,yerr,roi,iKV )
                
                if not fix_K:
                    K,M = plsq0
                    KArray[j,1],MArray[j,1] = paramStdDev #the fit deviation

                    if fit_gm:
                        K_gm,M_gm = plsq0_gm
                        KArray_gm[j,1],MArray_gm[j,1] = paramStdDev_gm #the fit deviation
                    
                else:
                    M = plsq0;K = initKValue
                    MArray[j,1] = paramStdDev
                    KArray[j,1] = np.std(trace[j,:])

                    if fit_gm:
                        M_gm = plsq0_gm;K_gm = initKValue
                        MArray_gm[j,1] = paramStdDev_gm
                        KArray_gm[j,1] = np.std(trace[j,:])
                    
                MArray[j,0],KArray[j,0] = M,K  #the fit res for M,K
                if fit_gm:
                    MArray_gm[j,0],KArray_gm[j,0] = M_gm,K_gm #the fit res for M,K
                
                beta[j,0] = 1./MArray[j,0]
                beta[j,1] = MArray[j,1]/MArray[j,0]**2

                RArray[j,0] = 2.*ydata[2]*(1.-ydata[1])/ydata[1]**2 - 1
                
                dataToSave[i,j]={}
                dataToSave[i,j]['M'] = {'value': MArray[j,0], 'stddev': MArray[j,1]}
                dataToSave[i,j]['K'] = {'value': KArray[j,0], 'stddev': KArray[j,1]}
                if fit_gm:
                    dataToSave[i,j]['M_gm'] = {'value': MArray_gm[j,0], 'stddev': MArray_gm[j,1]}
                    dataToSave[i,j]['K_gm'] = {'value': KArray_gm[j,0], 'stddev': KArray_gm[j,1]}
                
                dataToSave[i,j]['beta'] = {'value': beta[j,0], 'stddev': beta[j,1]}
                dataToSave[i,j]['R-1'] = {'value': RArray[j,0], 'stddev': RArray[j,1]}
                
##        cpdump(dataToSave,  sid
##            +'results_nq%s_qs%s.json'%(noqs,qstart), outDir=outDir)
##        
##        print 'the results file:\n  %s are saved in folder:\n  %s'%(sid
##                            +'results_nq%s_qs%s.json'%(noqs,qstart),outDir)
        return dataToSave



    def model_to_pds( self, model, qv=None):
        '''translate the model into a pandas frame'''
        import pandas as pd
        
        tn, qn = model.shape
        tindex=range( tn )
        
        tx=[2**i * exposuretime for i in range(tn)]
        if qv is None:qx=range(qn)
        else:qx=qv
        #qx=[i*qperpixel for i in qradi]
        
        qindex=range( qn )
        
        m_tcolumns = [ 'M%i'%n for n in range(tn) ]  #for M
        me_tcolumns = [ 'Me%i'%n for n in range(tn) ]

        b_tcolumns = [ 'B%i'%n for n in range(tn) ]  #for beta
        be_tcolumns = [ 'Be%i'%n for n in range(tn) ]

        k_tcolumns = [ 'K%i'%n for n in range(tn) ]  #for K
        ke_tcolumns = [ 'Ke%i'%n for n in range(tn) ]
        
        
        m_qcolumns = [ 'M%i'%n for n in range(qn) ]  #for M
        me_qcolumns = [ 'Me%i'%n for n in range(qn) ]

        b_qcolumns = [ 'B%i'%n for n in range(qn) ]  #for beta
        be_qcolumns = [ 'Be%i'%n for n in range(qn) ]

        k_qcolumns = [ 'K%i'%n for n in range(qn) ]  #for K
        ke_qcolumns = [ 'Ke%i'%n for n in range(qn) ]

        Mt = {}
        Met= {}
        Kt ={}
        Ket={}
        bt={}
        bet={}
        
        for q in range(qn):
            Mt[q]=[];Met[q]=[];Kt[q]=[];Ket[q]=[];bt[q]=[];bet[q]=[]
            for t in range(tn):
                Mt[q].append( model[t,q]['M']['value'] )
                Met[q].append( model[t,q]['M']['stddev'] )
                Kt[q].append( model[t,q]['K']['value'] )
                Ket[q].append( model[t,q]['K']['stddev'] )                                
                bt[q].append( model[t,q]['beta']['value'] )                                
                bet[q].append( model[t,q]['beta']['stddev'] )

        data = np.array( tx).reshape([tn,1])  #zeros( [ tn,1])
        for dicts in [Mt,Met,Kt,Ket,bt,bet]:            
            for q in range(qn):
                #print array(dicts[q]).shape
                data = np.hstack( [data, np.array( dicts[q]).reshape([tn,1]) ] )         
        columns = ['x'] + m_qcolumns + me_qcolumns +k_qcolumns +ke_qcolumns +b_qcolumns +be_qcolumns 
        df_t = pd.DataFrame(data=data, index=tindex,columns=columns)

        Mq = {}
        Meq= {}
        Kq ={}
        Keq={}
        bq={}
        beq={}
        
        for t in range(tn):
            Mq[t]=[];Meq[t]=[];Kq[t]=[];Keq[t]=[];bq[t]=[];beq[t]=[]
            for q in range(qn):
                Mq[t].append( model[t,q]['M']['value'] )
                Meq[t].append( model[t,q]['M']['stddev'] )
                Kq[t].append( model[t,q]['K']['value'] )
                Keq[t].append( model[t,q]['K']['stddev'] )                                
                bq[t].append( model[t,q]['beta']['value'] )                                
                beq[t].append( model[t,q]['beta']['stddev'] )
                               
        data = np.array( qx ).reshape([qn,1])  
        for dicts in [Mq,Meq,Kq,Keq,bq,beq]:            
            for t in range(tn):
                data = np.hstack( [data, np.array( dicts[t]).reshape([qn,1])] )         
        columns = ['x'] + m_tcolumns + me_tcolumns +k_tcolumns +ke_tcolumns +b_tcolumns +be_tcolumns 
        df_q = pd.DataFrame(data=data, index=qindex,columns=columns)
        
        return df_q,df_t

############################################





















