import pickle

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
    # plot_training_hist(hist)
    pre, re, acc = evaluate_model(trained_model, test)

    save_model(trained_model, model_name)

    return pre, re, acc


if __name__ == '__main__':
    statistics = {}
    for length in ["100", "300", "500", "1000"]:
        for trim in ["rndm", "None"]:
            model_name = length + "_" + length + "_" + trim
            path = "../data/images/train_data_" + model_name
            pre, re, acc = training(path, "model_" + model_name)
            statistics[model_name] = [pre, re, acc]

    with open('../evaluation/statistics_train_specific_length.pkl', 'wb') as f:
        pickle.dump(statistics, f)
    print(statistics)
