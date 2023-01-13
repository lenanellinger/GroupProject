import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import re
import os
import pickle

from LoadData import load_data
from Evaluation import evaluate_model


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
        yhat = trained_model.predict(np.expand_dims(im, axis=0))
        print("Actual: " + str(exon_intron_rates[i]) + ", Predicted: " + str(yhat[0]))


def evaluate_labeled_test_data(model_name, image_folder):
    trained_model = tf.keras.models.load_model('trainedModels/' + model_name)

    image_dataset = load_data("../data/images/" + image_folder)

    pre, re, acc = evaluate_model(trained_model, image_dataset)

    return pre, re, acc


if __name__ == '__main__':
    statistics = {}
    for length in ["100", "300", "500", "1000", "5000"]:
        for trim in ["rndm", "None"]:
            for level in range(1, 7):
                model_name = "level" + str(level) + "/model_" + length + "_" + length + "_" + trim + ".h5"
                image_folder = "test_data_realish_" + length + "_" + length + "_" + trim + "/level" + str(level)
                pre, re, acc = evaluate_labeled_test_data(model_name, image_folder)
                statistics[model_name] = [pre, re, acc]

    with open('../evaluation/statistics_test.pkl', 'wb') as f:
        pickle.dump(statistics, f)
    print(statistics)
