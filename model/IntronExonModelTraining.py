import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import re
import os

from Model import setup_model, train_model, save_model
from LoadData import load_data, get_train_val_test_dataset
from Evaluation import evaluate_model, plot_training_hist


def training(path, model_name):
    """
    trains a model with the data given in the path
    :param model_name: to save the trained model
    :param path: with sub folders intron and exon
    :return precision
    :return recall
    :return accuracy
    """
    model = setup_model()

    data = load_data(path)

    train_size = int(len(data) * .7)
    val_size = int(len(data) * .2)
    test_size = int(len(data) * .1)

    train, val, test = get_train_val_test_dataset(data, train_size, val_size, test_size)
    trained_model, hist = train_model(model, train, val)
    plot_training_hist(hist)
    pre, re, acc = evaluate_model(trained_model, test)

    save_model(trained_model, model_name)

    return pre, re, acc


def mixed_sequence_test():
    trained_model = tf.keras.models.load_model('trainedModels/IntronExonClassifier.h5')
    exon_intron_rates = []
    with open("data/sequenceData/1_exon_5000_length.txt") as file:
        for line in file:
            if line[0] == ">":
                ratio = re.findall(r'\d[.]\d+', line)
                exon_intron_rates.append(float(ratio[0]))

    folder_dir = "../data/images/mixed_inferno_level5"

    for i, image in enumerate(os.listdir(folder_dir)):
        im = plt.imread(os.path.join(folder_dir, image))[:, :, :3]
        print(im.shape)
        yhat = trained_model.predict(np.expand_dims(im, axis=0))
        print("Actual: " + str(exon_intron_rates[i]) + ", Predicted: " + str(yhat[0]))


if __name__ == '__main__':
    for length in ["100", "300", "500", "1000"]:
        for trim in ["rndm", "None"]:
            model_name = length + "_" + length + "_" + trim
            path = "../data/images/train_data_" + model_name
            pre, re, acc = training(path, "model_" + model_name)
