import os
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf

from tensorflow.keras.layers import Conv2D, MaxPooling2D, Dense, Flatten, Dropout
from tensorflow.keras.metrics import Precision, Recall, BinaryAccuracy
from tensorflow.keras.models import Sequential


def load_data(level):
    '''
    loads the image data for a specific level
    :param level
    :return train: training data set
    :return val: validation data set
    :return test: test data set
    '''
    # TODO: does images need to have a specific format (pixel)
    data = tf.keras.utils.image_dataset_from_directory('data/level' + str(level))

    train_size = int(len(data) * .7)
    val_size = int(len(data) * .2)
    test_size = int(len(data) * .1)

    # TODO: evtl crossvalidation (in each training epoch other validation set)
    train = data.take(train_size)
    val = data.skip(train_size).take(val_size)
    test = data.skip(train_size + val_size).take(test_size)

    return train, val, test


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
    model.add(Dense(2, activation='sigmoid'))

    model.compile('adam', loss=tf.losses.BinaryCrossentropy(), metrics=['accuracy'])
    print(model.summary())

    return model


def save_model(model):
    """
    save the model
    :param model: Tensorflow model
    """
    model.save(os.path.join('models', 'IntronExonClassifier.h5'))


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


def plot_training_hist(hist):
    """
    plots the loss and accuracy performance of the training process
    :param hist: history of training
    """
    # loss
    fig = plt.figure()
    plt.plot(hist.history['loss'], color='teal', label='loss')
    plt.plot(hist.history['val_loss'], color='orange', label='val_loss')
    fig.suptitle('Loss', fontsize=20)
    plt.legend(loc="upper left")
    plt.show()

    # accuracy
    fig = plt.figure()
    plt.plot(hist.history['accuracy'], color='teal', label='accuracy')
    plt.plot(hist.history['val_accuracy'], color='orange', label='val_accuracy')
    fig.suptitle('Accuracy', fontsize=20)
    plt.legend(loc="upper left")
    plt.show()


def evaluate_model(model, test):
    """
    evaluates the trained model
    :param model: trained model
    :param test: test data set
    """
    pre = Precision()
    re = Recall()
    acc = BinaryAccuracy()

    for batch in test.as_numpy_iterator():
        X, y = batch
        yhat = model.predict(X)
        pre.update_state(y, yhat)
        re.update_state(y, yhat)
        acc.update_state(y, yhat)

    print(pre.result(), re.result(), acc.result())


def predict(model, image):
    """
    predicts a class for a given image
    :param model: trained model
    :param image: to be predicted
    :return:
    """
    plt.imshow(image)
    plt.show()

    yhat = model.predict(np.expand_dims(image / 255, 0))

    # TODO: which class is which?
    if yhat > 0.5:
        print(f'Predicted class is Intron')
    else:
        print(f'Predicted class is Exon')


if __name__ == '__main__':
    model = setup_model()
    save_model(model)
    # load_model('imageclassifier.h5')

    train, val, test = load_data(1)
    trained_model, hist = train_model(model, train, val)
    plot_training_hist(hist)
    evaluate_model(trained_model, test)
