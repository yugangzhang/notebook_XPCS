from math import pi


def get_Lambda(E,u='SI'):
    """
    calculates X-ray wavelength as a function of Energy [keV] in optional units.
    Syntax: getLambda(E,u), 
    where E=X-ray energy; optional: u= 'A','nm','um','cm','mm','m','SI' (='m'), default in the absence of u: 'SI'  
    
    """
    hPlank=6.62606876e-34;
    cvac=2.99792458e8;
    Qelectron=1.602176463e-19;
    scale=1
    #l=hPlank*cvac/(E*1000*Qelectron)
    l=hPlank*cvac/(E*1000*Qelectron);
    if u is 'A':
            scale=1e10;return l*scale # Angstroem
    elif u is 'nm':
            scale=1e9; return l*scale # nm
    elif u is 'um':
            scale=1e6; return l*scale # um
    elif u is 'mm':
            scale=1e3; return l*scale # mm
    elif u is 'cm':
            scale=1e2; return l*scale # cm
    elif u is 'm' or u is 'SI':
            scale=1; return l*scale
    else:
            print 'invalid option, type "get_Lambda(\'?\')" for available options and syntax'

   

def qpix(x,lambda_ = 1.37760, Xenergy=None, Ldet = 5000. ):

    #9 KeV >>>  1.37760 A
    ''' DOCUMENT qpix(x)
       x is pixel size, 75X75 um for eiger 1/4M  
       retuns q (in 1/A), x in mm   
       q=2kSin(tth/2)=k*tth=2pi/lam*tth=(2 pi/lam) * (x/L)
       lambda=1.37760 A
       L=5000 mm (5 m) have to look in the book for the exact one	
    '''
    if Xenergy is not None:lambda_ = get_Lambda(Xenergy,u='SI')
    lambda_=lambda_
    Ldet=Ldet	
    qq=(2*pi/lambda_)*(x/Ldet);
    return qq










 
