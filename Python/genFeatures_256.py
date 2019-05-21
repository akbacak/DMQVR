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
    np.save(open("Python/Outputs/" + file + ".npy", 'w'), XX)



#coding=utf-8
from keras.models import Sequential,Input,Model,InputLayer
from keras.models import model_from_json
from keras.models import load_model
import numpy as np    # for mathematical operations
import os



np.set_printoptions(linewidth=8192)


json_file = open('/home/ubuntu/keras/enver/dmlvh2/models/dmlvh2_mLSTM_256_model.json', 'r')
model = json_file.read()
json_file.close()
model = model_from_json(model)
model.load_weights("/home/ubuntu/keras/enver/dmlvh2/models/dmlvh2_mLSTM_256_weights.h5")
model = Model(inputs=model.input, outputs=model.get_layer('dense_1').output) # dense_1 for features , dense_2 for predictions



X1 = []
X1 = np.load(open("Python/Outputs/video_1.npy"))
X1.shape
print(X1.shape[:])
X1 = X1.reshape(1, X1.shape[0], X1.shape[1] * X1.shape[2] * X1.shape[3] )
binary_codes_1 = model.predict(X1, batch_size=64, verbose=0, steps=None)
binary_codes_1 = binary_codes_1 
binary_codes_1 = binary_codes_1.astype(float)


X2 = []
X2 = np.load(open("Python/Outputs/video_2.npy"))
X2.shape
print(X2.shape[:])
X2 = X2.reshape(1, X2.shape[0], X2.shape[1] * X2.shape[2] * X2.shape[3] )
binary_codes_2 = model.predict(X2, batch_size=64, verbose=0, steps=None)
binary_codes_2 = binary_codes_2
binary_codes_2 = binary_codes_2.astype(float)

X3 = []
X3 = np.load(open("Python/Outputs/video_3.npy"))
X3.shape
print(X3.shape[:])
X3 = X3.reshape(1, X3.shape[0], X3.shape[1] * X3.shape[2] * X3.shape[3] )
binary_codes_3 = model.predict(X3, batch_size=64, verbose=0, steps=None)
binary_codes_3 = binary_codes_3 
binary_codes_3 = binary_codes_3.astype(float)



np.savetxt('Python/Outputs/features_q1.txt', binary_codes_1,  fmt='%f',  delimiter=',')
np.savetxt('Python/Outputs/features_q2.txt', binary_codes_2,  fmt='%f',  delimiter=',')
np.savetxt('Python/Outputs/features_q3.txt', binary_codes_3,  fmt='%f',  delimiter=',')
