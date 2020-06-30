#!/usr/bin/env python
# coding: utf-8

# In[1]:


import tensorflow as tf
from tensorflow.python.client import device_lib


# In[2]:


print(device_lib.list_local_devices())


# In[3]:


print(tf.config.list_physical_devices('GPU'))


# In[4]:


print(tf.test.is_built_with_cuda())


# In[5]:


print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))


# In[ ]:




