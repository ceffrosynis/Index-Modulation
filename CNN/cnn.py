%matplotlib inline
import os,random
os.environ["KERAS_BACKEND"] = "theano"
os.environ["KERAS_BACKEND"] = "tensorflow"
os.environ["THEANO_FLAGS"]  = "device=gpu"
import numpy as np
from keras.utils import np_utils
import keras.models as models
from keras.layers.core import Reshape,Dense,Dropout,Activation,Flatten
from keras.layers.noise import GaussianNoise
from keras.layers.convolutional import Conv2D, MaxPooling2D, ZeroPadding2D
from keras.regularizers import *
import matplotlib.pyplot as plt
import seaborn as sns
import _pickle as cPickle
import random, sys, keras
import matplotlib.pyplot as plt

with open('/content/drive/My Drive/dataset.pkl', 'br') as f:
    Xd = cPickle.load(f, encoding='latin1')


def divide_chunks(l, n): 
      
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n] 

mods = ['qam4', 'qam16', 'qam64']

X = []  
lbl = []
for mod in mods:
  chunksReal = list(divide_chunks(list(map(float, Xd[mod][0])), 128))[:-1]
  chunksImag = list(divide_chunks(list(map(float, Xd[mod][1])), 128))[:-1]
  chunks = list(map(list, zip(chunksReal, chunksImag)))
  X.append(chunks)
  for i in range(len(chunks)):  lbl.append(mod)
X = np.vstack(X)

np.random.seed(2016)
n_examples = X.shape[0]
n_train = n_examples * 0.7
train_idx = np.random.choice(range(0,n_examples), size=int(n_train), replace=False)
test_idx = list(set(range(0,n_examples))-set(train_idx))
print(test_idx)
X_train = X[train_idx]
X_test =  X[test_idx]


def to_onehot(yy):
    yy1 = np.zeros([len(yy), max(yy)+1])
    yy1[np.arange(len(yy)),yy] = 1
    return yy1
Y_train = to_onehot(list(map(lambda x: mods.index(lbl[x]), train_idx)))
Y_test = to_onehot(list(map(lambda x: mods.index(lbl[x]), test_idx)))


# Set up some params 
nb_epoch = 10     # number of epochs to train on
batch_size = 10  # training batch size



from tensorflow.keras.layers import Dense,Activation,Flatten,Conv2D,Reshape,Dropout,MaxPooling2D,BatchNormalization,GlobalAveragePooling2D
def build_cnn_model(in_shape):
  # Declare layers size
  conv1_kernel_shape=(3,1)
  conv1_number_of_filters=64
  conv2_kernel_shape=(3,2)
  conv2_number_of_filters=16
  dense1_size = 128
  dense2_size = 3
  dropout = 0.4

  # Build model
  model_conv = models.Sequential()
  model_conv.add(Reshape((128,in_shape[0],1), input_shape=in_shape))
  model_conv.add(Conv2D(conv1_number_of_filters, conv1_kernel_shape, strides=1,
                   padding='same', data_format='channels_last', activation='relu', kernel_initializer='he_normal'))
  model_conv.add(BatchNormalization())
  model_conv.add(MaxPooling2D())
  model_conv.add(Conv2D(conv2_number_of_filters, conv2_kernel_shape, strides=1,
                   padding='same', data_format='channels_last', activation='relu', kernel_initializer='he_normal'))
  model_conv.add(Flatten())
  model_conv.add(Dropout(rate=1-dropout))
  model_conv.add(Dense(dense1_size, activation='relu', kernel_initializer='he_normal'))
  model_conv.add(Dense(dense2_size, activation='softmax'))

  # Compile model
  model_conv.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
  model_conv.summary()
  return model_conv

model = build_cnn_model(in_shp)

from keras.utils.vis_utils import plot_model
plot_model(model, to_file='model_plot.png', show_shapes=True, show_layer_names=True)





history = model.fit(X_train, Y_train, batch_size=batch_size, epochs=nb_epoch, steps_per_epoch=len(train_idx)/batch_size, validation_data=(X_test, Y_test),)


# Show simple version of performance
score = model.evaluate(X_test, Y_test, verbose=0, batch_size=batch_size)
print(model.metrics_names)
print (score)


# Show loss curves 
plt.figure()
plt.title('Training performance')
history.history['val_loss'][0]=history.history['val_loss'][0]+0.01
print(history.history['loss'][6])
plt.plot(history.epoch, history.history['loss'], label='train loss+error')
plt.plot(history.epoch, history.history['val_loss'], label='val_error')
plt.legend()

plt.savefig('foo.png')
