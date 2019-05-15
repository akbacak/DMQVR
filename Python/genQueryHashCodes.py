#coding=utf-8
from keras.models import Sequential,Input,Model,InputLayer
from keras.models import model_from_json
from keras.models import load_model
import numpy as np    # for mathematical operations
import os



np.set_printoptions(linewidth=8192)


json_file = open('models/dmlvh2_mLSTM_512_2_model.json', 'r')
model = json_file.read()
json_file.close()
model = model_from_json(model)
model.load_weights("models/dmlvh2_mLSTM_512_2_weights.h5")
model = Model(inputs=model.input, outputs=model.get_layer('dense_1').output) # dense_1 for features , dense_2 for predictions


X1 = []
X1 = np.load(open("Python/video_1.npy"))
X1.shape
print(X1.shape[:])
X1 = X1.reshape(1, X1.shape[0], X1.shape[1] * X1.shape[2] * X1.shape[3] )
binary_codes_1 = model.predict(X1, batch_size=64, verbose=0, steps=None)
binary_codes_1 = binary_codes_1 > 0.5
binary_codes_1 = binary_codes_1.astype(int)

X2 = []
X2 = np.load(open("NPY_videos2/video_2.mp4.npy"))
X2.shape
print(X2.shape[:])
X2 = X2.reshape(1, X2.shape[0], X2.shape[1] * X2.shape[2] * X2.shape[3] )
binary_codes_2 = model.predict(X2, batch_size=64, verbose=0, steps=None)
binary_codes_2 = binary_codes_2 > 0.5
binary_codes_2 = binary_codes_2.astype(int)



np.savetxt('Python/prd_q1.txt', binary_codes_1,  fmt='%d',  delimiter=',')
np.savetxt('Python/prd_q2.txt', binary_codes_2,  fmt='%d',  delimiter=',')
