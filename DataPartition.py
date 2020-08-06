#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np


# In[1]:



def DataPartition(Data2Partion_dict, Partition=[0.6,0.2,0.2], RanSeed=None):
    assert np.sum(Partition)==1, 'Partition must be a list of fractions (Train_Frac, Val_Frac, Test_Frac) whose sum is equal to 1'
    np.random.seed(RanSeed)
    
    Ascd_ClassSizes = sorted([len(items[1]) for items in Data2Partion_dict.items()])
    len_min = Ascd_ClassSizes[0]
    len_max = Ascd_ClassSizes[-1]
    
    Test_size = round(Partition[2]*len_min)
    Val_size = round(Partition[1]*len_min)
    Train_size = len_max - (Val_size+Test_size)
    
    PartitionedData_dict = {}
    for key in Data2Partion_dict.keys():
        new_order = np.random.choice(Data2Partion_dict[key], size=len(Data2Partion_dict[key]), replace = False)
        Tr_Set_temp=[]; Val_Set=[]; Ts_Set=[];
        Ts_Set = new_order[:Test_size]
        Val_Set = new_order[Test_size:Test_size+Val_size]
        Tr_Set_temp = new_order[Test_size+Val_size:]
        Tr_Set=[]
        quotient=Train_size//len(Tr_Set_temp) ; remainder=Train_size%len(Tr_Set_temp)
        for count in range(quotient):
            Tr_Set.extend(Tr_Set_temp)
        Tr_Set.extend(np.random.choice(Tr_Set_temp, size=remainder, replace = False))
        PartitionedData_dict[key]={'Tr_Set':Tr_Set,'Val_Set':Val_Set,'Ts_Set':Ts_Set}
        
    return PartitionedData_dict


# In[ ]:





# In[ ]:





# In[ ]:




