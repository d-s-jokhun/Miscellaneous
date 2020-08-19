#!/usr/bin/env python
# coding: utf-8

# In[1]:


import zipfile


# In[ ]:


def namelist_in_archive (archive):
    with zipfile.ZipFile(archive) as archive:
        namelist = archive.namelist()
    return namelist

