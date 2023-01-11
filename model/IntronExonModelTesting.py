import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import re
import os


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
    mixed_sequence_test()
