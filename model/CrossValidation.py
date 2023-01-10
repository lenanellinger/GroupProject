import os
import numpy as np
import tensorflow as tf

from sklearn.model_selection import KFold, StratifiedKFold

from Model import setup_model
from LoadData import load_data, get_train_test_dataset
from Evaluation import evaluate_model


def get_model_name(level, k):
    return 'model_level' + str(level) + "_fold" + str(k) + '.h5'


def cross_validation(data, level):
    """
    validates the model with k-fold cross validation
    :param data: whole dataset
    :param level
    :return average validation accuracy
    :return average validation loss
    """
    validation_accuracy = []
    validation_loss = []

    kf = KFold(n_splits=5)
    skf = StratifiedKFold(n_splits=5, random_state=7, shuffle=True)
    fold = 1

    train_size = int(len(data) * .9)
    test_size = int(len(data) * .1)
    train_dataset, test_dataset = get_train_test_dataset(data, train_size, test_size)
    train_images = np.concatenate(list(train_dataset.map(lambda x, y: x)))
    train_labels = np.concatenate(list(train_dataset.map(lambda x, y: y)))

    for train, val in skf.split(train_images, train_labels):
        model = setup_model()
        checkpoint = tf.keras.callbacks.ModelCheckpoint(os.path.join('trainedModels', get_model_name(level, fold)),
                                                        monitor='val_accuracy', verbose=1,
                                                        save_best_only=True, mode='max')
        callbacks_list = [checkpoint]
        model.fit(train_images[train], train_labels[train], epochs=20, validation_split=0.2,
                         callbacks=callbacks_list)
        model.load_weights(os.path.join('trainedModels', get_model_name(level, fold)))
        results = model.evaluate(train_images[val], train_labels[val])
        results = dict(zip(model.metrics_names, results))

        validation_accuracy.append(results['accuracy'])
        validation_loss.append(results['loss'])

        tf.keras.backend.clear_session()

        fold += 1

    return np.mean(validation_accuracy), np.mean(validation_loss)


if __name__ == '__main__':
    average_val_accuracy = []
    average_val_loss = []

    for level in range(1, 7):
        data = load_data(level)
        acc, loss = cross_validation(data, level)
        average_val_accuracy.append(acc)
        average_val_loss.append(loss)

    print()
    print("\t| Acccuracy\t| Loss")
    for i in range(6):
        print(str(i+1) + "\t| " + str(round(average_val_accuracy[i], 2)) + "\t\t| " + str(round(average_val_loss[i], 2)))
