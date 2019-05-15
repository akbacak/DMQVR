import cv2
import keras
import numpy as np
from keras.applications import VGG16
from keras.models import load_model
import tensorflow as tf
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
import glob

listing = os.listdir("Python/Frames/")

count = 1

for file in listing:
    listing_2  = os.listdir("Python/Frames/" + file + "/" )

    X = []
    for images in listing_2:
        image =  plt.imread("Python/Frames/" + file + "/" + images )
        X.append (image)
    X = np.array(X)
    print(X.shape)

    image_size=224,
    base_model = VGG16(weights='imagenet', include_top=False, input_shape=(X.shape[1:]))

    batch_size = 48
    XX = base_model.predict(X, batch_size=batch_size, verbose=0, steps=None)
    print(XX.shape)
    np.save(open("Python/" + file + ".npy", 'w'), XX)
    count+=1







'''
# https://stackoverflow.com/questions/42389870/how-to-put-many-numpy-files-in-one-big-numpy-file-file-by-file?rq=1
# This code merge multiple npy files in to one

import matplotlib.pyplot as plt 
import numpy as np
import glob
import os, sys
fpath ="/home/ubuntu/keras/enver/dmlvh/bf.npy"
npyfilespath ="/home/ubuntu/keras/enver/dmlvh/preprocessed_videos/"   

os.chdir(npyfilespath)
npfiles= glob.glob("*.npy")
npfiles.sort()
all_arrays = []
for i, npfile in enumerate(npfiles):
    all_arrays.append(np.load(os.path.join(npyfilespath, npfile)))
np.save(fpath, np.concatenate(all_arrays))


'''
