import tensorflow as tf


def load_data(level, cmap):
    '''
    loads the image data for a specific level
    :param level
    :param cmap
    :return data: whole dataset
    '''
    data = tf.keras.utils.image_dataset_from_directory('data/' + cmap + '/level' + str(level),
                                                       image_size=(64, 64),
                                                       batch_size=32)
    return data


def get_train_val_test_dataset(data, train_size, val_size, test_size):
    '''
    creates training, validation and test data set of given size

    :param data: loaded data
    :param train_size: size of training dataset
    :param val_size: size of validation dataset
    :param test_size: size of test dataset

    :return train: training data set
    :return val: validation data set
    :return test: test data set
    '''
    train = data.take(train_size)
    val = data.skip(train_size).take(val_size)
    test = data.skip(train_size + val_size).take(test_size)

    return train, val, test


def get_train_test_dataset(data, train_size, test_size):
    '''
    creates training and test data set of given size for cross validation

    :param data: loaded data
    :param train_size: size of training dataset
    :param test_size: size of test dataset

    :return train: training data set
    :return test: test data set
    '''
    train = data.take(train_size)
    test = data.skip(train_size).take(test_size)

    return train, test
