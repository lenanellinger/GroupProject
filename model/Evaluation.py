import numpy as np
import matplotlib.pyplot as plt

import tensorflow as tf
from tensorflow.keras.metrics import Precision, Recall, BinaryAccuracy


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
    :param test: test images set
    """
    pre = Precision()
    re = Recall()
    acc = BinaryAccuracy()

    for batch in test.as_numpy_iterator():
        X, y = batch
        yhat = model.predict(X)
        print(yhat)
        pre.update_state(y, yhat)
        re.update_state(y, yhat)
        acc.update_state(y, yhat)

    tf.print("Precision: ", pre.result())
    tf.print("Recall: ", re.result())
    tf.print("Accuracy: ", acc.result())

    return float(pre.result()), float(re.result()), float(acc.result())


def predict(model, image):
    """
    predicts a class for a given image
    :param model: trained model
    :param image: to be predicted
    """
    plt.imshow(image)
    plt.show()

    yhat = model.predict(np.expand_dims(image[:, :, :3] / 255, 0))
    if yhat > 0.5:
        print(f'Predicted class is Intron')
    else:
        print(f'Predicted class is Exon')
