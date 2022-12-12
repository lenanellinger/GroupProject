import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Dense, Flatten, Dropout
import os


def loadData(level):
    '''

    :param level
    :return:
    '''
    data = tf.keras.utils.image_dataset_from_directory('data/level' + str(level))

    train_size = int(len(data) * .7)
    val_size = int(len(data) * .2)
    test_size = int(len(data) * .1)


    # evtl crossvalidation
    train = data.take(train_size)
    val = data.skip(train_size).take(val_size)
    test = data.skip(train_size + val_size).take(test_size)


def setUpModel():
    """
    recreate model 4 of paper
    :return: model with layer
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

def saveModel(model):
    """
    save the model
    :param model: Tensorflow model
    """
    model.save(os.path.join('models', 'IntronExonClassifier.h5'))


if __name__ == '__main__':
    model = setUpModel()
    saveModel(model)
    # load_model('imageclassifier.h5')
    print("hello")