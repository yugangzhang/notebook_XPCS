# -*- coding: utf-8 -*-
######################################################################################
#####################Make Mask##########################################
#########By Dr. Yugang Zhang, yuzhang@bnl.gov, Octa 18, 2015##############################################
#######################################################################################
######################################################################################
############################################################################################################
#######################################################################################
######################################################################################
################################################################################
 
from pylab import *
from matplotlib.path import Path
import matplotlib.patches as patches
#from GetEdf import get_edf
from numpy import indices
import scipy.misc, scipy.ndimage

class ROI(object):
    """ control_button= 'control' for the segment mode;
        control_button= 'shift' for the free draw mode;
    """
    def __init__(self, axes, fig):
        self.previous_point = []
        self.start_point = []
        self.end_point = []
        self.line = None  
        #self.linem = None 
        self.fig =  fig
        self.fig.canvas.draw()
        
        self.seg_point={}
        self.path_point={}
        self.path_point_key = 0
        self.seg_point_key = 0        
        self.path_point['0']=[]
        self.seg_point['0']=[]
        
    def motion_notify_callback(self, event):
        ''' control_button= 'shift' for the free draw mode;'''
        if event.inaxes:
            #print event.key, event.button
            ax = event.inaxes
            x, y = event.xdata, event.ydata
            if event.button == None and self.line != None: # Move line around
                #print 'button=None, here'
                self.line.set_data([self.previous_point[0], x],
                                   [self.previous_point[1], y])      
                self.fig.canvas.draw()
            elif event.button == 1 and event.key == 'shift': # Free Hand Drawing
                if self.line == None: # if there is no line, create a line
                    #print 'No line here'                        
                    self.line = Line2D([x,  x], [y, y],marker='o' )                        
                    self.start_point = [x,y]
                    self.previous_point =  self.start_point
                    ax.add_line(self.line)
                    self.fig.canvas.draw()
                    self.path_point['%s'%self.path_point_key].append( [x,y] )
                    #self.path_point.append( [x,y] )
                else:
                    #print 'free move here'
                    self.line = Line2D([self.previous_point[0], x],
                                  [self.previous_point[1], y])                 
  
                    ax.add_line(self.line)
                    self.previous_point = [x, y]
                    self.fig.canvas.draw()
                    #self.path_point.append( [x,y] )
                    self.path_point['%s'%self.path_point_key].append( [x,y] )
                    
            elif (event.button == 3)  and (self.line != None) and (event.key == 'shift'): # close the loop
                self.line.set_data([self.previous_point[0],  self.start_point[0]],
                                   [self.previous_point[1], self.start_point[1]])
                print 'close the loop'
                ax.add_line(self.line)
                self.fig.canvas.draw()                
                self.line = None
                self.path_point_key+=1
                self.path_point['%s'%self.path_point_key]=[]
                

    def button_press_callback(self, event):
        """ control_button= 'control' for the segment mode;"""
        #key = self.seg_point_key
        #print key
        if event.inaxes:            
            x, y = event.xdata, event.ydata            
            ax = event.inaxes            
            if (event.key == 'control') and (event.button == 1):                
                # If you press the right button
                    if self.line == None: # if there is no line, create a line
                        #print 'No line here'                        
                        self.line = Line2D([x,  x], [y, y],marker='o' )                        
                        self.start_point = [x,y]
                        self.previous_point =  self.start_point
                        ax.add_line(self.line)
                        self.fig.canvas.draw()
                        self.seg_point['%s'%self.seg_point_key].append( [x,y] )
                        #self.seg_point.append( [x,y] )
                    else: # if there is a line, create a segment
                        
                        self.line = Line2D([self.previous_point[0], x],
                                           [self.previous_point[1], y],
                                           marker='o')
                        self.previous_point = [x,y]
                        ax.add_line(self.line)
                        self.fig.canvas.draw()
                        self.seg_point['%s'%self.seg_point_key].append( [x,y] )
                        #self.seg_point.append( [x,y] )
                        #print 'make a line here'

            elif (event.button == 3)  and (self.line != None) and (event.key == 'control'): # close the loop
                self.line.set_data([self.previous_point[0],  self.start_point[0]],
                                   [self.previous_point[1], self.start_point[1]])
                print 'close the loop'
                ax.add_line(self.line)
                self.fig.canvas.draw()                
                self.line = None
                self.seg_point_key+=1
                self.seg_point['%s'%self.seg_point_key]=[]
                #print self.seg_point_key

def get_path2( data, verts ):
    path1 = Path(verts)
    index = path1.contains_points(data[:,:2])
    #print data[index, :2]
    
    plot(data[:,0],data[:,1], 'b.')
    patch = patches.PathPatch(path1, facecolor='orange', lw=2)
    gca().add_patch(patch)
    plot(data[index,0], data[index,1], 'r.')
    show()
    #return index

def get_mask( img, verts):
    #verts is a dict
    dx,dy=img.shape
    mask = zeros( dx*dy, dtype=bool)
    data = zeros( [dx,dy,2])
    #print data.shape,data[:,:,0].shape
    #data[:,:,0]= range(dx)
    data[:,:,0]= range(dy)
    #for i in range(dy):
    for i in range(dx):
        data[i,:,1] = i     
    data=data.reshape( [dx*dy,2])    
    for k in verts.keys():
        path1 = Path(verts[k])
        #print verts[k],data
        mask = mask | path1.contains_points(data)
        #mask = path1.contains_points(data)
        print mask.shape, len(where(mask)[0])
    mask=mask.reshape( [dx,dy])
    return mask

def show_mask( img, mask, only_mask=False, show_point=True ):
    img_ = img.copy()
    if only_mask:img_[:]=0
    img_[mask]=100 *0
    #print img_, mask
    dx,dy=img.shape    
    fig, ax = plt.subplots(nrows=1)
    vmin=img_.min();vmax=img_.max()
    #ax.set_xlim(0,dx-1)
    #ax.set_ylim(0,dy-1)
    if show_point:
        data = zeros( [dx,dy,2])
        data[:,:,0]= range(dy)
        for i in range(dx):
            data[i,:,1] = i 
        #print data
        data=data.reshape( [dx*dy,2])
        mask=mask.reshape( dx*dy,)
        plot(data[mask,0], data[mask,1], 'r.')
    ax.imshow( img_,cmap=plt.get_cmap('winter'), vmin=vmin,vmax=2.0 )    
    
    show()

def show_mask2(mask):
    dimy,dimx=mask.shape
    y,x = indices( [ dimy,dimx] )
    x[:]=0
    fig, ax = plt.subplots(nrows=1)
    #vmin=img_.min();vmax=img_.max()
    ax.set_xlim(0,dimx-1)
    ax.set_ylim(0,dimy-1)
    x[mask]=1
    ax.imshow( x )
    show()

def get_path( img, verts ):
    #verts=array( verts, int)
    path1 = Path(verts)
    print path1
    dx,dy=img.shape
    data = zeros( [dx,dy,2])
    data[:,:,0]= range(dx)
    for i in range(dy):
        data[i,:,1] = i 
    #print data
    data=data.reshape( [dx*dy,2])
    #print data.shape
    #print path1,data
    index = path1.contains_points(data)
    print index.shape, len(where(index)[0])
    #print data[index, :2]
    
    #plot(data[:,0],data[:,1], 'b.')
    fig, ax = plt.subplots(nrows=1)
    vmin=img.min();vmax=img.max()
    ax.imshow( img,cmap=plt.get_cmap('winter'), vmin=vmin,vmax=vmax )
    #ax.set_xlim(0,dx-1)
    #ax.set_ylim(0,dy-1)
    patch = patches.PathPatch(path1, facecolor='orange', lw=2)
    gca().add_patch(patch)
    plot(data[index,0], data[index,1], 'r.')
    show() 

def conver_to_int(dict):
    d={}
    for k in dict.keys():
        if dict[k]!=[]:
            a = array(dict[k],int)
            d[k]=a
    return d










    

