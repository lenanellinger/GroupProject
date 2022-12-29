from Model import setup_model, train_model
from LoadData import load_data, get_train_val_test_dataset
from Evaluation import evaluate_model, plot_training_hist


if __name__ == '__main__':
    level = 5
    model = setup_model()

    data = load_data(level)

    train_size = int(len(data) * .7)
    val_size = int(len(data) * .2)
    test_size = int(len(data) * .1)

    train, val, test = get_train_val_test_dataset(data, train_size, val_size, test_size)
    trained_model, hist = train_model(model, train, val)
    # plot_training_hist(hist)
    evaluate_model(trained_model, test)
