import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import re
import os

from Model import setup_model, train_model, save_model
from LoadData import load_data, get_train_val_test_dataset
from Evaluation import evaluate_model, plot_training_hist


def train_process():
    level = 5
    cmap = "inferno"
    model = setup_model()

    data = load_data(level, cmap)

    train_size = int(len(data) * .7)
    val_size = int(len(data) * .2)
    test_size = int(len(data) * .1)

    train, val, test = get_train_val_test_dataset(data, train_size, val_size, test_size)
    trained_model, hist = train_model(model, train, val)
    plot_training_hist(hist)
    pre, re, acc = evaluate_model(trained_model, test)


def test_process():
    trained_model = tf.keras.models.load_model('trainedModels/IntronExonClassifier.h5')
    exon_intron_rates = []
    with open("data/sequenceData/generator/1_exon_5000_length.txt") as file:
        for line in file:
            if line[0] == ">":
                ratio = re.findall(r'\d[.]\d+', line)
                exon_intron_rates.append(float(ratio[0]))

    folder_dir = "../data/images/mixed_inferno_level5/inferno/level5/mixed_inferno_level5"

    for i, image in enumerate(os.listdir(folder_dir)):
        im = plt.imread(os.path.join(folder_dir, image))[:, :, :3]
        print(im.shape)
        yhat = trained_model.predict(np.expand_dims(im, axis=0))
        print("Actual: " + str(exon_intron_rates[i]) + ", Predicted: " + str(yhat[0]))


if __name__ == '__main__':
    # train_process()
    test_process()
