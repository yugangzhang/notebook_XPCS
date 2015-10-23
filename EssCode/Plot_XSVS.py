# -*- coding: utf-8 -*-
######################################################################################
#######################################################################################
#####################Do XSVS_Plot######### ####################################
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
1) xsvs_plot_histogram( bins,his, dev, trace_mean, 
     model=None, title=None, show=True,show_poisson = False,
     show_bn_K = False,show_bn=True,show_gamma=False, qlist=None,
     print_K=False, scalex=False, show_dev=True, xlim=None,  ylim=None,                     

     output='res_his.png', outDir=None,plotq=None,)

'''

import itertools
import matplotlib.pyplot as plt
#from pylab import *
from XSVS_Functions import nbinomPMF,poisson,peval, gammaDist
from numpy import sort,zeros,arange,sqrt,array, linspace

mcolors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
markers = itertools.cycle(list(plt.Line2D.filled_markers))
lstyles = itertools.cycle(['-', '--', '-.','.',':'])




def xsvs_plot_histogram( bins,his, dev, trace_mean, 
                         model=None, title=None, show=True,
            show_poisson = False,show_bn_K = False,show_bn=True,
                         show_gamma=False,  expose_time = 1.0,
                         qlist=None,print_K=False,
                    scalex=False, show_dev=True, xlim=None,  ylim=None,                     
               output='res_his.png', outDir=None,plotq=None,):
    #plot the bins (normalized by tra.mean), his, devation,
    #his.shape = [ timebin_lev, noqs ] 
    
    fig = plt.figure()
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
    if plotq is None:fig.set_size_inches(20,10.5)    
    else:fig.set_size_inches(20./2,12.5/1.6)
    
    timebin_lev, noqs = his.shape

    if plotq is not None:noqs=1
    sx= int( round (sqrt(noqs)) )
    if noqs%sx==0:sy=noqs/sx
    else:sy=noqs/sx+1
    
    time_list = [2**i*expose_time for i in range(timebin_lev)]
    exp_list = time_list
    for qn in range( noqs ):               
        ax = fig.add_subplot(sx,sy,qn+1 )        
        fig.subplots_adjust(wspace = 0.5)
        fig.subplots_adjust(hspace = 0.5) 
        
        for tn in xrange(len(exp_list)):
            marker = markers.next()
            color = mcolors.next()

            if plotq is not None:qn=plotq
            #print qn
            histogramBins = array( bins[tn,qn] )
            hist_data = his[tn,qn]
            hist_err =  dev[tn,qn]
            #if tn==0:print dev[0,0]

            if model is not None:
                x = linspace(histogramBins[0],histogramBins[-1],num=300)
                #print x.shape
                K_= model[tn,qn]['K']['value']
                M= model[tn,qn]['M']['value']
                fit_nbinomPMF = peval(x,[K_,M])

                K = trace_mean['q%i'%qn] * 2**(tn)  #using real count average
                
                fit_poisson = poisson(x,K)
                if show_gamma:
                    K_gm= model[tn,qn]['K_gm']['value']
                    M_gm= model[tn,qn]['M_gm']['value']
                    fit_gamma = gammaDist(x,[K_gm,M_gm])
                
                fit_nbinomPMF_ = peval(x,[K,M])
                if print_K:print 'the fit-K and count K are:  %s--%s'%(K_, K)
            # Plot the histogram for the first q value for all exposures\n",

            # First, get the q value\n",
            #q_val = q_value_list[j]
            if qlist is not None:
                qv = qlist[qn]
                q_txt = r'$q= %.5f$ $\AA^{-1}$'%(qv) 
                #q_txt = r'$q = %.5f $' % (qv)   
             
            else:
                #q_txt = r'$q = Q_%s $' % (qn+1)
                q_txt = r'$q= Q_%s $'%(qn+1) 
            
            #q_txt = r'$q = %.2f \times 10^{-2} \; \AA^{-1}$' % (q_val*1e2)
            # Get the exposure time\n",
            exp_val =  tn
            
            if title is not None:q_txt=r'%s '%title+q_txt            
            ax.set_title(q_txt,fontsize= 20 )
            # Plot the data 
            #x-axis is K/<K>, y-axis is P(K), title = q_txt
            
            if scalex!=False:
##                K_ave=0
##                tsi,tbin,tqs,tfra = tra.shape
##                for si in range(tsi):
##                    K_ave += tra[si][tn,qn].mean() * 2**(tn+1)
                #K_ave= model[tn,qn]['K']['value']
                if model is None:K_= trace_mean['q%i'%qn] * 2**(tn) 
                K_ave=  K_
            else:K_ave=1            
            #print len( histogramBins[:-1]/K_ave)
            #print len(hist_data[:-1])
            #print len(hist_err[:-1])
            if not show_dev:hist_err_ = hist_err *0
            else:hist_err_ = hist_err
            if model is not None:ls=''
            else:ls='--'

            #print (histogramBins[:-1]/K_ave).shape
            #print ls
             
            ax.errorbar(histogramBins[:-1]/K_ave,hist_data,
                     hist_err_, ls=ls, markersize=10, lw=2,                          
                marker=marker,color=color,label=r'%.3f' % (exp_list[tn]))
            
            if model is not None:
                if show_bn:
                    ax.errorbar(x/K_ave,fit_nbinomPMF ,
                    fit_nbinomPMF*0, ls='-',  lw=4,color=color,)
                                #label='nbinomPMF')
                
                if show_gamma:
                    ax.errorbar(x/K_ave,fit_gamma,
                    fit_gamma*0, ls=':',  lw=6, color=color,)
                                #label='gamma')                    
                    
                if show_poisson:
                    #print 'show_poisson as -- line'
                    ax.errorbar(x/K_ave,fit_poisson,
                    fit_poisson*0, ls='--',  lw=6, color=color,)
                                #label='poisson')
                if show_bn_K:
                    ax.errorbar(x/K_ave,fit_nbinomPMF_ ,
                    fit_nbinomPMF*0, ls='-.',  lw=4,color=color,)
                               # label='nbinomPMF_Kcons')
                    

    
        #ax.set_xlim(-1,20)
        if xlim is not None:
            x1,x2=xlim
            ax.set_xlim( x1,x2)
            
        if ylim is not None:
            y1,y2=ylim
            ax.set_ylim( y1,y2)
            
        #ax.set_ylim(1e-4,5)    
        #ax.set_yscale('log')        
        #ax.legend(frameon=0,bbox_to_anchor=(1.2,1),title=r'$t_e$')
        if noqs==1:fontsize=20
        else:fontsize=8
        ax.legend(frameon=0,bbox_to_anchor=(1.0,1),
                  title=r'$t_e$',fontsize=fontsize)
        
        #ax.legend(frameon=0,bbox_to_anchor=(1.1,1),title=r'$t_e$')
        if scalex!=False:ax.set_xlabel(r'$K/ \langle K \rangle$',fontsize=20)
        else:ax.set_xlabel(r'$K$',fontsize=20)
        ax.set_ylabel(r'$P(K)$',fontsize=20)
        
        
    #ax.legend(frameon=0,bbox_to_anchor=(1.0,1),title=r'$t_e$',fontsize=12)
##    legend = ax.legend(loc='upper right', shadow=True)
##    frame = legend.get_frame()
##    frame.set_facecolor('0.90')
##    # Set the fontsize
##    for label in legend.get_texts():
##        label.set_fontsize(10)
##    for label in legend.get_lines():
##        label.set_linewidth(1.5)  # the legend line width

            
    fig.tight_layout()       
    if outDir!=None:output =outDir + output
    #print (output)
    plt.savefig( output )    
    if show:plt.show()



 

