import cv2
import os
import numpy as np
import imutils
import tensorflow as tf

from tensorflow.keras.preprocessing.image import img_to_array
from tensorflow.keras.preprocessing.image import load_img
from tensorflow.keras.models import Model


"""
Class Activation Map
inspired by https://pyimagesearch.com/2020/03/09/grad-cam-visualize-class-activation-maps-with-keras-tensorflow-and-deep-learning/
"""


class GradCAM:
    def __init__(self, model, layer_name):
        self.model = model
        self.layerName = layer_name

    def compute_heatmap(self, image, eps=1e-8):
        # compute gradients
        grad_model = Model(
            inputs=[self.model.inputs],
            outputs=[self.model.get_layer(self.layerName).output,
                     self.model.output])
        with tf.GradientTape() as tape:
            inputs = tf.cast(image, tf.float32)
            (conv_outputs, predictions) = grad_model(inputs)
            loss = predictions[:, 0]
        grads = tape.gradient(loss, conv_outputs)
        cast_conv_outputs = tf.cast(conv_outputs > 0, "float32")
        cast_grads = tf.cast(grads > 0, "float32")
        guided_grads = cast_conv_outputs * cast_grads * grads
        conv_outputs = conv_outputs[0]
        guided_grads = guided_grads[0]

        # create heatmap
        weights = tf.reduce_mean(guided_grads, axis=(0, 1))
        cam = tf.reduce_sum(tf.multiply(weights, conv_outputs), axis=-1)
        (w, h) = (image.shape[2], image.shape[1])
        heatmap = cv2.resize(cam.numpy(), (w, h))
        numer = heatmap - np.min(heatmap)
        denom = (heatmap.max() - heatmap.min()) + eps
        heatmap = numer / denom
        heatmap = (heatmap * 255).astype("uint8")
        return heatmap

    def overlay_heatmap(self, heatmap, image, alpha=0.5,
                        colormap=cv2.COLORMAP_VIRIDIS):
        # create overlay heatmap
        heatmap = cv2.applyColorMap(heatmap, colormap)
        output = cv2.addWeighted(image, alpha, heatmap, 1 - alpha, 0)
        return heatmap, output


def class_activation_map(level, model_name, img_folder, image_name, layer_name):
    # load model and image
    model = tf.keras.models.load_model('trainedModels/' + "level" + str(level) + "/" + model_name)
    image = os.path.join(img_folder, image_name)
    orig = cv2.imread(image)
    image = load_img(image, target_size=(64, 64))
    image = img_to_array(image)
    image = np.expand_dims(image, axis=0)
    prediction = model.predict(image)
    if prediction > 0.5:
        label = "intron"
    else:
        label = "exon"

    # create class activation map
    cam = GradCAM(model, layer_name)
    heatmap = cam.compute_heatmap(image)
    heatmap = cv2.resize(heatmap, (orig.shape[1], orig.shape[0]))
    heatmap, output = cam.overlay_heatmap(heatmap, orig, alpha=0.5)
    cv2.putText(output, label, (1, 6), cv2.FONT_HERSHEY_SIMPLEX,
                0.3, (255, 255, 255), 1)
    output = np.vstack([orig, heatmap, output])
    output = imutils.resize(output, height=700)
    folder = img_folder + "/CAM_" + model_name.replace(".h5", "")
    is_exist = os.path.exists(folder)
    if not is_exist:
        os.makedirs(folder)
    if not cv2.imwrite(os.path.join(folder, image_name), output):
        raise Exception("Could not write image")
    # cv2.imshow("Output", output)
    # cv2.waitKey(0)


if __name__ == '__main__':
    for length in ["100", "300"]:
        for type in ["/exon", "/intron"]:
            for level in range(3, 5):
                folder_dir = "../data/images/train_data_" + length + "_" + length + "_rndm/level" + str(level) + type
                for image in os.listdir(folder_dir):
                    class_activation_map(level, "model_" + length + "_" + length + "_rndm.h5", folder_dir, image,
                                         "conv2d_2")
