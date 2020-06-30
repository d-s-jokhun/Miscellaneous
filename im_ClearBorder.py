#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Takes in a labelled image (labelled mask) and outputs a copy of the image but without the objects which are touching the image border.


# In[2]:


import cv2
import numpy as np


# In[3]:


def im_ClearBorder (labelled_img):
    n_objs = np.amax(labelled_img)
    if n_objs==0:
        return (labelled_img)
    
    labelled_ClearedBorder=np.zeros(labelled_img.shape, dtype=bool)
    
    for n in range(n_objs):
        obj_n_mask = labelled_img==n+1
        if np.sum(obj_n_mask)>0:
            BoundRec = cv2.boundingRect(np.uint8(obj_n_mask))
            x1=BoundRec[0]
            x2=BoundRec[0]+BoundRec[2]-1
            y1=BoundRec[1]
            y2=BoundRec[1]+BoundRec[3]-1
            canvas_h,canvas_w = labelled_img.shape
            if x1>0 and x1<canvas_w-1 and x2>0 and x2<canvas_w-1 and y1>0 and y1<canvas_h-1 and y2>0 and y2<canvas_h-1:
                labelled_ClearedBorder += obj_n_mask
            
    labelled_ClearedBorder = labelled_ClearedBorder.astype(labelled_img.dtype) * labelled_img
        
    return labelled_ClearedBorder


# In[ ]:




