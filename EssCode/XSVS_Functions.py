# -*- coding: utf-8 -*-
######################################################################################
################################################################################
#####################Fitting####################################################
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
1) nbinomPMF(x,K,M):
    Negative Binomial (Poisson-Gamma) distribution function
2)residuals(params,y,x):
    Residuals function used for least squares fitting
3)residuals_err(params,y,x,yerr):
    Residuals function used for least squares fitting
4)residuals_logerr(params,y,x,yerr):
    Residuals function used for least squares fitting
5)residuals2(params,y,x,K):
    Residuals function used for least squares fitting with
    *K* parameter fixed.
6)residuals2_err(params,y,x,yerr,K):
    Residuals function used for least squares fitting with
    *K* parameter fixed
7)residuals2_logerr(params,y,x,yerr,K):
    Residuals function used for least squares fitting with
    *K* parameter fixed.
8)gammaDist(x,params):
    Gamma distribution function
9)residuals_gamma(params,y,x)
10) residuals_gamma_err(params,y,x,yerr)
11)residuals2_gamma(params,y,x,K)
12)residuals2_gamma_err(params,y,x,yerr,K)
13)poisson(x,K)  
    Poisson distribution function.
14)gauss(x, p)
15)residuals_gauss(params,y,x):
16)beta_Time( x, params)
17)residuals_beta_Time(params,y,x)
18)residuals_ploy(params,y,x)
19)residuals_ploy2(params,y,x)
20)g2_tau( x, params)
21)residuals_g2_tau(params,y,x)
22) peval(x,params)
23) HisMtdData( data, bins): 
    histrogram a multi-dimensional data, e.g.,with shape as M*N
'''



from scipy.special import gamma, gammaln
#from sympy import Symbol
import scipy.misc
import sys
from time import time

from numpy import pi,sin,arctan,sqrt,mgrid,where,shape,exp,linspace,std,arange
from numpy import power,log,log10,array,zeros,reshape,mean,histogram,round,int_
from numpy import indices,hypot,digitize,ma,histogramdd,apply_over_axes,sum
from numpy import around,intersect1d, ravel, unique,hstack,vstack,zeros_like
from numpy.linalg import lstsq
from numpy import polyfit,poly1d,ones
import pandas as pd




#all the functions for SVS


########nbinomPMF
def nbinomPMF(x,K,M):
    '''
    Negative Binomial (Poisson-Gamma) distribution function.
    K is  average photon counts,
    M is the number of coherent modes    
    the probability density of photon, P(x), satify this  function.    
    '''
    K = float(K)
    M = float(M)
    coeff = exp(gammaln(x+M)-gammaln(x+1)-gammaln(M))
    Pk = coeff*power(M/(K+M),M)
    coeff2 = power(K/(M+K),x)
    Pk *= coeff2
    return Pk

def peval(x,params):
    K,M = params
    pr = M/(K+M)
    result = nbinomPMF(x,K,M)
    return result

def residuals(params,y,x):
    '''Residuals function used for least squares fitting
    '''    
    K,M = params
    K = float(K)
    M = float(M)
    pr = M/(K+M)
    result = ( y - nbinomPMF(x,K,M) )
    return result
def residuals_err(params,y,x,yerr):
    '''Residuals function used for least squares fitting
    '''    
    K,M = params
    K = float(K)
    M = float(M)
    pr = M/(K+M)
    result = ( y - ( nbinomPMF(x,K,M)))/yerr 
    return result
def residuals_logerr(params,y,x,yerr):
    '''Residuals function used for least squares fitting
    '''
    
    K,M = params
    K = float(K)
    M = float(M)
    pr = M/(K+M)
    result = ( y - log10( nbinomPMF(x,K,M)))/yerr 
    return result

def residuals2(params,y,x,K):
    '''Residuals function used for least squares fitting with
    *K* parameter fixed.
    '''
    M = params
    pr = M/(K+M)
    result = (y - nbinomPMF(x,K,M) )
    return result
def residuals2_err(params,y,x,yerr,K):
    '''Residuals function used for least squares fitting with
    *K* parameter fixed.
    '''
    M = params
    K = float(K)
    M = float(M)
    pr = M/(K+M)
    result = ( y - ( nbinomPMF(x,K,M)))/yerr 
    return result
def residuals2_logerr(params,y,x,yerr,K):
    '''Residuals function used for least squares fitting with
    *K* parameter fixed.
    '''
    M = params
    K = float(K)
    M = float(M)
    pr = M/(K+M)
    result = ( y - log10( nbinomPMF(x,K,M)))/yerr 
    return result

#########gamma
def gammaDist(x,params):
    '''Gamma distribution function
    M,K = params, where K is  average photon counts <x>,
    M is the number of coherent modes,
    In case of high intensity, the beam behavors like wave and
    the probability density of photon, P(x), satify this gamma function.
    '''
    
    K,M = params
    K = float(K)
    M = float(M)
    coeff = exp(M*log(M) + (M-1)*log(x) - gammaln(M) - M*log(K))
    Gd = coeff*exp(-M*x/K)
    return Gd
def residuals_gamma(params,y,x):
    
    K,M = params
    K = float(K)
    M = float(M)
    pr = M/(K+M)
    result = (  y - ( gammaDist(x,[K,M]))  )
    return result

def residuals_gamma_err(params,y,x,yerr):
    
    K,M = params
    K = float(K)
    M = float(M)
    pr = M/(K+M)
    result = (  y - ( gammaDist(x,[K,M]))  )/yerr 
    return result

def residuals2_gamma(params,y,x,K):
    
    M = params
    K = float(K)
    M = float(M)
    pr = M/(K+M)
    result = (  y - ( gammaDist(x,[K,M]))  )
    return result

def residuals2_gamma_err(params,y,x,yerr,K):
    
    M = params
    K = float(K)
    M = float(M)
    pr = M/(K+M)
    result = (  y - ( gammaDist(x,[K,M]))  )/yerr 
    return result


#########poisson
def poisson(x,K):    
    '''Poisson distribution function.
    K is  average photon counts    
    In case of low intensity, the beam behavors like particle and
    the probability density of photon, P(x), satify this poisson function.
    '''
    K = float(K)    
    Pk = exp(-K)*power(K,x)/gamma(x+1)
    return Pk

###################gauss
def gauss(x, p): 
    import numpy as np
    (yo, A, xc, w,) = p
    return yo + A * np.exp(-(x - xc) ** 2 / (2.0 * w ** 2))

def residuals_gauss(params,y,x): 
    result = (  y - gauss( x, params)  )
    return result



def beta_Time( x, params):
    gamma, beta1, betai = params
    bt = beta1 * ( exp(-2*gamma*x) -1. + 2*gamma*x)/(2*gamma**2*x**2) + betai
    return bt

def residuals_beta_Time(params,y,x):
    
    gamma, beta1, betai = params 
    result = (  y - beta_Time( x, params)  )
    return result    


def residuals_ploy(params,y,x):
    f=poly1d(params)
    return y - f(x)

def residuals_ploy2(params,y,x):
    f= params*x**2
    return y - f(x)


def g2_tau( x, params):
    gamma, beta1, gi = params
    g2 = beta1 * exp(-2*gamma*x) + gi
    return g2

def residuals_g2_tau(params,y,x):
    
    gamma, beta1, gi = params
    result = (  y - g2_tau( x, params)  )
    return result


####for histrogram

def HisMtdData( data, bins):
    '''histrogram a multi-dimensional data, e.g.,with shape as M*N
    Input: data --array
           bins --with shape as N
    Output: histrogram, bins_edge as list
    #####current not work for the ring shaped data, because of different length
    '''
    H, edges = histogramdd( data, bins) #H has shape as bins
    axis = len( H.shape )
    his = list( zeros(axis))
    for ax in range(axis):
        all_ax = list( range(axis))
        all_ax.pop( ax )
        his_ax= ( apply_over_axes(sum, H, all_ax) )
        his[ax]=array( his_ax.flatten(), dtype = int )
    return his,edges


































