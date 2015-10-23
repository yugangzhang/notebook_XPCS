def load_eiger_mask( sid, inDir):
    '''
    Give the eiger master path and filename to get piel mask
    Input:
        sid: filename
        InDir: filepath
    Output:
        maskps: a NDarray, the mask file,
    
    '''
    
    import h5py
    import numpy as np
    
    fp = inDir + sid
    f= h5py.File( fp, 'r' )
    ins = f['entry']['instrument']
    det = ins['detector']['detectorSpecific']
    pmask = det['pixel_mask']
    ps = np.array( pmask)
    ps[ np.where(ps!=0)]= 100
    ps[ np.where(ps==0)]= 0
    ps[ np.where(ps==100)]= 1
    maskps = np.array( ps, dtype=bool)

    return maskps


