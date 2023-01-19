import numpy as np
import tensorflow as tf
import csv

from Model import setup_model
from LoadData import load_data, get_train_test_dataset


def get_model_name(level, k):
    return 'model_level' + str(level) + "_fold" + str(k) + '.h5'


def cross_validation(data, level, cmap):
    """
    validates the model with k-fold cross validation
    :param cmap: color map
    :param level
    :param data: whole dataset
    :return average validation accuracy
    :return average validation loss
    """
    k = 5

    train_size = int(len(data) * .9)
    test_size = int(len(data) * .1)
    train_dataset, test_dataset = get_train_test_dataset(data, train_size, test_size)

    validation_accuracy_folds = []
    fold_len = int(len(train_dataset) / k)
    for val_fold in range(k):
        print()
        print("Level: " + str(level))
        print("Cmap: " + cmap)
        print("Fold: " + str(val_fold))
        val = train_dataset.skip(val_fold * fold_len).take(fold_len)

        train_1 = train_dataset.take(val_fold * fold_len)
        train_2 = train_dataset.skip((val_fold + 1) * fold_len).take((k - val_fold - 1) * fold_len)
        train = train_1.concatenate(train_2)

        model = setup_model()
        model.fit(train, epochs=20, validation_data=val)
        results = model.evaluate(val)
        results = dict(zip(model.metrics_names, results))

        validation_accuracy_folds.append(results['accuracy'])

        tf.keras.backend.clear_session()

    return np.mean(validation_accuracy_folds)


def cv_images_colormap():
    average_val_accuracy = {}
    color_maps = ["Greys", "cividis", "viridis", "plasma", "inferno", "magma"]

    for cmap in color_maps:
        average_val_accuracy[cmap] = []
        for level in range(1, 7):
            data = load_data('images/' + cmap + '/level' + str(level))
            acc = cross_validation(data, level, cmap)
            average_val_accuracy[cmap].append(acc)

    with open('cv_colormap_level.csv', 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        csv_writer.writerow(["cmap", "1", "2", "3", "4", "5", "6"])
        for cmap in color_maps:
            csv_writer.writerow([cmap] + average_val_accuracy[cmap])

    print()
    print("\t| 1\t| 2\t| 3\t| 4\t| 5\t| 6")
    print("---------------------------------------")
    for cmap in color_maps:
        print(cmap + "\t| " + str(round(average_val_accuracy[cmap][0], 2)) + "\t|"
              + str(round(average_val_accuracy[cmap][1], 2)) + "\t|"
              + str(round(average_val_accuracy[cmap][2], 2)) + "\t|"
              + str(round(average_val_accuracy[cmap][3], 2)) + "\t|"
              + str(round(average_val_accuracy[cmap][4], 2)) + "\t|"
              + str(round(average_val_accuracy[cmap][5], 2)))


if __name__ == '__main__':
    cv_images_colormap()
