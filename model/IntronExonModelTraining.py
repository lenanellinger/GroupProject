import pickle
import numpy as np

from Model import setup_model, train_model, save_model
from LoadData import load_data, get_train_val_test_dataset
from Evaluation import evaluate_model, plot_training_hist


def training(path, model_name=None):
    """
    trains a model with the data given in the path
    :param model_name: to save the trained model, if None, not saved
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
    # plot_training_hist(hist)
    pre, re, acc = evaluate_model(trained_model, test)

    if model_name is not None:
        save_model(trained_model, model_name)

    return pre, re, acc


if __name__ == '__main__':
    statistics = {}
    for length in ["100", "300", "500", "1000", "5000"]:
        for trim in ["rndm", "None"]:
            for level in range(1, 7):
                model_name = length + "_" + length + "_" + trim
                path = "../data/images/train_data_" + model_name + "/level" + str(level)
                pre, re, acc = training(path, "level" + str(level) + "/model_" + model_name)
        print(statistics)

    with open('../evaluation/statistics_train.pkl', 'wb') as f:
        pickle.dump(statistics, f)
    print(statistics)
