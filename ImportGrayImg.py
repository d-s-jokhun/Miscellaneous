#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import cv2
import math


# In[ ]:



def ImportGrayImg (path, ScaleFactor=1, ImgSize=0):
       
    img = cv2.imread(path, cv2.IMREAD_UNCHANGED).astype(np.float32)
        
    W = math.floor(img.shape[1] * ScaleFactor)
    H = math.floor(img.shape[0] * ScaleFactor)
    
    if ImgSize>0 and (W>ImgSize or H>ImgSize):
        return
    
    img = cv2.resize(img, (W,H))
    img = np.subtract(img, np.amin(img))
    img = np.divide(img, np.amax(img))

    
    width_to_pad=ImgSize-img.shape[1]
    height_to_pad=ImgSize-img.shape[0]
    
    width_start, width_end, height_start, height_end = (0,0,0,0)
    
    if width_to_pad>0:
        width_start = width_to_pad//2
        width_end = width_to_pad - width_start
        
    if height_to_pad>0:
        height_start = height_to_pad//2
        height_end = height_to_pad - height_start
        
    if len(img.shape) == 2:
        img = np.pad(img,((height_start,height_end),(width_start,width_end)))
    elif len(img.shape) == 3:
        img = np.pad(img,((height_start,height_end),(width_start,width_end),(0,0)))
    
    return img

