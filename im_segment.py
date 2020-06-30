#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Imput = (a mask (or labelled mask) and an image) or (a mask or labelled mask)
# Outputs a crop if the mask is not labelled
# Outputs a list of crops if the mask is labelled


# In[2]:


import cv2
import numpy as np


# In[5]:


def im_segment (mask, image=None):
        
    n_objs = np.amax(mask)
    if n_objs==0:
        return (None) 
    
    if type(image) == type(None):
        image = mask
    
    Segments=[]
    for n in range(n_objs):
        obj_n_mask = mask==n+1
        if np.sum(obj_n_mask)>0:
            BoundRec = cv2.boundingRect(np.uint8(obj_n_mask))
            x1=BoundRec[0]
            x2=BoundRec[0]+BoundRec[2]-1
            y1=BoundRec[1]
            y2=BoundRec[1]+BoundRec[3]-1

            obj_n_img = image * obj_n_mask
            crop = obj_n_img[y1:y2,x1:x2]
            Segments.append(crop)
    
    return Segments


# In[ ]:




