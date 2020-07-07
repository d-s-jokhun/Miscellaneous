#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import math


# In[41]:



def DataPartition(FullSet, Partition=[0.6,0.2,0.2], RanSeed=None):
    assert np.sum(Partition)<=1
    np.random.seed(RanSeed)
    l = len(FullSet)
    new_order = np.random.choice(FullSet, size = l, replace = False)

    Tr_MaxIdx = np.ceil(l*Partition[0]).astype(np.int)
    Val_MaxIdx = np.ceil(l*(Partition[0]+Partition[1])).astype(np.int)
    Ts_MaxIdx = np.ceil(l*(Partition[0]+Partition[1]+Partition[2])).astype(np.int)

    Tr_Set = new_order[0:Tr_MaxIdx]
    Val_Set = new_order[Tr_MaxIdx:Val_MaxIdx]
    Ts_Set = new_order[Val_MaxIdx:Ts_MaxIdx]
    
    return Tr_Set, Val_Set, Ts_Set


# In[ ]:





# In[ ]:





# In[ ]:




