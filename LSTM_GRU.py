
"""
This is the machine learning part for Solar heat gain research
The machine learning algorithm is Long short-term memory (LSTM) and Gated recurrent unit (GRU)

Model build based on Keras and Tensorflow

"""

#%%
# import packages
from keras.engine import sequential
from keras.layers.core import Dropout
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf
import keras
import sklearn
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_absolute_error

# Data preprocessing
"""
Import dataframe
Here, the feature selection for this dataset is done.
"""
# import data
data_frame_path = r"/Users/qianchengsun/PhD/solar_heat_gain/DF_15_with_solar.csv"
# read as dataframe
df_house_15 = pd.read_csv(data_frame_path)
# make sure the data information 
"""
Target :
    Temp : represent for the indoor temperature (F)

"""
dataset = df_house_15[["Temp",
            "Temp.Outdoor.Avg",
            "Temp.Dewpoint.Avg",
            "GHI",                                    
            "Cool.Setpoint",
            "Cool.Demand.Status",
            "beam_solar_radiation_west",                        
            "beam_solar_radiation_east",
            "cloud.cover"]]
# correlation plot
"""
Here is the step for correlation plot for the dataframe after the feature selection
default method for correlation plot is pearson
"""

correlation = dataset.corr(method = "pearson")
# show the correlation plot by using color map
# cmap represent for the color map in matplotlib.pyplot library
# reference website:
# https://matplotlib.org/stable/tutorials/colors/colormaps.html
correlation.style.background_gradient(cmap = "inferno")


# normalized the dataframe
"""
Here is the step for normalized the dataframe,
the method is making all the values from (0, 1) by using minimum and maximum method to scale
"""
# set up the scaler
scaler = MinMaxScaler(feature_range= (0,1))
data = scaler.fit_transform(dataset)

# fix randowm seed for reproducibility
np.random.seed(7)
# set up look_back
# time tp look_back for each prediction
look_back = 60
# number of column target in scaled data
col_target = 0
# number of input features in scaled data
n_input_feature = data.shape[1] - 1
# reshape data from 2D into 3D
length_data = len(data) - look_back

# create 3D numpy array
sample = np.ndarray(shape = (length_data, look_back, n_input_feature))
target = np.ndarray(shape = (length_data, 1, 1))

for i in range(length_data):
    sample[i, 0:look_back, 0:(n_input_feature - 1)] = data[i : i + look_back, 1:n_input_feature]
    target[i, :, :] = data[i + look_back, col_target]
# split the dataset into training dataset and testing dataset
# divide the dataset into 90% training and 10% testing
train_index = int(len(data) * 0.9)
test_index = len(data) - train_index

train_x = sample[0:train_index, :, :]
train_y = target[0:train_index, :, :]

test_x = sample[train_index:len(data), :, :]
test_y = target[train_index:len(data), :, :]


# create and fit the LSTM
model = tf.keras.Sequential()
model.add(tf.keras.layers.LSTM(units = 60,
            activation= "tanh",
            dropout = 0.2,
            recurrent_activation="sigmoid",
            recurrent_dropout = 0.2,
            input_shape = (look_back, n_input_feature),
            kernel_regularizer = tf.keras.regularizers.l1_l2(l1 = 0.01, l2 = 0.01),
            return_sequences = True))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.LSTM(units = 35,
            activation= "tanh",
            dropout = 0.1,
            recurrent_activation="sigmoid",
            recurrent_dropout = 0.1,
            kernel_regularizer = tf.keras.regularizers.l1_l2(l1 = 0.01, l2 = 0.01),
            return_sequences = True))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.LSTM(units = 10,
            activation= "tanh",
            dropout = 0.1,
            recurrent_activation="sigmoid",
            recurrent_dropout = 0.1,
            kernel_regularizer = tf.keras.regularizers.l1_l2(l1 = 0.01, l2 = 0.01),
            return_sequences = False))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(units = 1, activation = "relu"))

# loss function to do regression problems
# Reference website
# https://scikit-learn.org/stable/modules/model_evaluation.html
model.compile(loss = "mean_absolute_error", optimizer = "adam")

model.fit(x = sample[0:train_index, :, :],
        y = target[0: train_index, :, :],
        epochs = 100, 
        batch_size = 60,
        verbose = 2)

# %%
