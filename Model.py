import os
import tensorflow as tf

from tensorflow.keras.layers import Conv2D, MaxPooling2D, Dense, Flatten, Dropout
from tensorflow.keras.models import Sequential


def setup_model():
    """
    recreate model 4 of paper
    :return: model with specific layer
    """
    model = Sequential()

    # convolutional layer (filter number, filter size, step/stride)
    model.add(Conv2D(16, (5, 5), 5, activation='relu', input_shape=(64, 64, 3)))

    # max pooling layer (pool size, step/stride)
    model.add(MaxPooling2D((5, 5), 5))

    # four fully connected layer (neurons/units)
    model.add(Flatten())
    model.add(Dense(512, activation='relu'))
    model.add(Dense(256, activation='relu'))
    model.add(Dense(128, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))

    model.compile('adam', loss=tf.losses.BinaryCrossentropy(), metrics=['accuracy'])
    # print(model.summary())

    return model


def save_model(model, name):
    """
    save the model
    :param model: Tensorflow model
    """
    model.save(os.path.join('models', name + '.h5'))


def train_model(model, train, val):
    """
    trains the given model
    :param model: neural network model
    :param train: training data set
    :param val: validation data set
    :return: trained model
    :return: history of training process
    """
    # logging:
    logdir = 'logs'
    tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=logdir)

    # training
    hist = model.fit(train, epochs=20, validation_data=val, callbacks=[tensorboard_callback])

    return model, hist