
#%%
import tensorflow as tf
import os

os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]="0"


#%%

model = tf.keras.applications.NASNetLarge(
    input_shape=None,
    include_top=True,
    weights=None,
    input_tensor=None,
    pooling=None,
    classes=1000,
)

#%%



