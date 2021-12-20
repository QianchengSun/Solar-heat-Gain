
"""
This is the machine learning part for Solar heat gain research
The machine learning algorithm is Long short-term memory (LSTM) and Gated recurrent unit (GRU)

Model build based on Keras and Tensorflow


In this example, the epoch = 10 is used to test if the code can run successfully

Package usage:

    pandas : which is used for the data preprocessing for the dataframe

    numpy, tensorflow, keras, sklearn : is used for machine learning model build and training 
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
from sklearn.metrics import r2_score


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
"""
# show the correlation plot by using color map
# cmap represent for the color map in matplotlib.pyplot library
# reference website:
# https://matplotlib.org/stable/tutorials/colors/colormaps.html
"""
correlation.style.background_gradient(cmap = "inferno")


# normalized the dataframe
"""
Here is the step for normalized the dataframe,
the method is making all the values from (0, 1) by using minimum and maximum method to scale

# set up the scaler
# Reference website:
# https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html
"""
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
"""
# LSTM layer build
# Reference website:
# https://www.tensorflow.org/api_docs/python/tf/keras/layers/LSTM
"""
model.add(tf.keras.layers.LSTM(units = 60,
            activation= "tanh",
            dropout = 0.2,
            recurrent_activation="sigmoid",
            recurrent_dropout = 0.2,
            input_shape = (look_back, n_input_feature),
            kernel_regularizer = tf.keras.regularizers.l1_l2(l1 = 0.01, l2 = 0.01),
            return_sequences = True))
"""
# Batch Normalization
# Reference website:
# https://arxiv.org/abs/1502.03167

# Batch Normalization for LSTM model
# Reference website
# https://www.tensorflow.org/api_docs/python/tf/keras/layers/BatchNormalization
"""
model.add(tf.keras.layers.BatchNormalization())
"""
Build GRU layer
Reference website:
https://www.tensorflow.org/api_docs/python/tf/keras/layers/GRU
"""
model.add(tf.keras.layers.GRU(units = 35,
            activation= "tanh",
            dropout = 0.1,
            recurrent_activation="sigmoid",
            recurrent_dropout = 0.1,
            kernel_regularizer = tf.keras.regularizers.l1_l2(l1 = 0.01, l2 = 0.01),
            return_sequences = True))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.GRU(units = 10,
            activation= "tanh",
            dropout = 0.1,
            recurrent_activation="sigmoid",
            recurrent_dropout = 0.1,
            kernel_regularizer = tf.keras.regularizers.l1_l2(l1 = 0.01, l2 = 0.01),
            return_sequences = False))
model.add(tf.keras.layers.BatchNormalization())
"""
# dense layer for LSTM
# Reference website :
# https://www.tensorflow.org/api_docs/python/tf/keras/layers/Dense
"""
model.add(tf.keras.layers.Dense(units = 1, activation = "relu"))
"""
# loss function to do regression problems
# Reference website
# https://scikit-learn.org/stable/modules/model_evaluation.html
"""
model.compile(loss = "mean_absolute_error", 
            optimizer = "adam", 
            metrics = [tf.keras.metrics.MeanAbsoluteError()])

history = model.fit(x = sample[0:train_index, :, :],
        y = target[0: train_index, :, :],
        epochs = 10, 
        batch_size = 60,
        validation_split = 0.2, # here setup the validation dataset
        verbose = 2,
        # setup the early stopping for machine learning model
        # Reference website:
        # https://keras.io/api/callbacks/early_stopping/
        callbacks = tf.keras.callbacks.EarlyStopping(monitor = "loss",
                                                    min_delta = 0,
                                                    patience = 2,
                                                    verbose = 1,
                                                    mode = "auto",
                                                    restore_best_weights = True))
# visualize the loss function
"""
Data visualization for Deep learning Model using Matplotlit.pyplot

Reference website:

https://www.pluralsight.com/guides/data-visualization-deep-learning-model-using-matplotlib

Here develop the loss function plot in order to decide if the model has overfitting or underfitting
"""
loss_train = history.history["loss"]
loss_validation = history.history["val_loss"]
epochs = range(0, len(loss_train))
plt.plot(epochs, loss_train, color = "g", label = "Training_loss")
plt.plot(epochs, loss_validation, color = "r", label = "Validation_loss")
plt.title("Training and Validation loss")
plt.xlabel("Epochs")
plt.ylabel("Loss value")
plt.legend()
plt.show()

# visualize the mean absolute error 
train_mae = history.history["mean_absolute_error"]
valiedation_mae = history.history["val_mean_absolute_error"]
plt.plot(epochs, train_mae, color = "g", label = "Training_MAE")
plt.plot(epochs, valiedation_mae , color = "r", label = "Validation_MAE")
plt.title("Training and Validation MAE")
plt.xlabel("Epochs")
plt.ylabel("MAE value")
plt.legend()
plt.show()

# make predictions
"""
Apply the LSTM model on the data set
"""
trainPredict = model.predict(train_x)
# scale back the predictions
testPredict = model.predict(test_x)

# Scale back the prediction 
# find the min value and the max value of the dataset
"""
Reference website :
This link is shown about how to scale the data
https://scikit-learn.org/stable/modules/preprocessing.html

From the Mathmatical formulation, the scale back can be calculated
"""
min_value = np.min(dataset)[col_target]
max_value = np.max(dataset)[col_target]
spread = max_value - min_value
# scale back the value for the target prediction
trainPredict_scaleback = trainPredict * spread + min_value
testPredict_scaleback = testPredict * spread + min_value

# scale back to the value for training dataset
trainY = train_y.reshape((len(train_y), 1)) # reshape the data from 3D to 2D
trainY_scaleback = trainY * spread + min_value
# scale back to the value for testing dataset
testY = test_y.reshape((len(test_y), 1))
testY_scaleback = testY * spread + min_value
# calculate Mean squared error
"""
Reference website:
https://scikit-learn.org/stable/modules/generated/sklearn.metrics.mean_absolute_error.html
"""
# obtain the mean absolute error for training dataset
train_score = mean_absolute_error(trainY_scaleback[:,0], trainPredict_scaleback[:,0])
print("Train Score : MAE:", (train_score))
# obtain the mean absolute error for testing dataset
test_score = mean_absolute_error(testY_scaleback[:,0], testPredict_scaleback[:,0])
print("Test score : MAE :", (test_score))

# obtain the R-squared value for training dataset
"""
Reference website :
https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html
"""
train_r2 = r2_score(trainY_scaleback[:,0], trainPredict_scaleback[:,0])
print("Train Score : R-squared value :", (train_r2))
# obtain the R-squared value for testing  dataset
test_r2 = r2_score(testY_scaleback[:,0], testPredict_scaleback[:,0])
print("Test Score : R-squared value :", (test_r2))

# %%
# plot the model result 
# create an empty array for Predicted value
trainPredictPlot = np.empty(shape = (len(data), 1))
trainPredictPlot[:, :] = np.nan
# predicted data is depend on 3 parts : look_back, train, test
trainPredictPlot[0 : look_back,:] = np.array(dataset["Temp"][0 : look_back]).reshape((look_back, 1))
trainPredictPlot[(look_back + 1) : (train_index + look_back + 1), :] = trainPredict_scaleback
trainPredictPlot[(train_index + look_back) : len(data) ] = testPredict_scaleback
# shift test predictions for plotting
originPlot = np.empty(shape = (len(data), 1))
originPlot[:, :] = np.array(dataset["Temp"]).reshape((len(data), 1))
# create plot for the prediction vs actual
plt.plot(trainPredictPlot, color = "green", label = "Predicted_Indoor_Temperature F")
plt.plot(originPlot, color = "red", label = "Actual_Indoor_Temperature F")
plt.title("Predicted Temperature vs Actual Temperature")
plt.xlabel("index value")
plt.ylabel("Temperature F ")
plt.legend()
plt.show()



# %%
