#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#  Imports using wget.download, unzips it and delete the downloaded zip file
# dwnld_args = (url,local_destination)


# In[9]:


from wget import download
from zipfile import ZipFile
from os import path, remove


# In[8]:


def wgetDwnld_n_unzip(dwnld_args):
    download(dwnld_args)
    NewFile=path.join(dwnld_args[1],path.basename(dwnld_args[0]))
    with ZipFile(NewFile, 'r') as ZipObj:
        ZipObj.extractall(dwnld_args[1])
    remove(NewFile)
    


# In[10]:





# In[ ]:




