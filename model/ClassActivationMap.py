from tensorflow.keras.models import Model
import tensorflow as tf
import numpy as np
import cv2
import os
from tensorflow.keras.applications import ResNet50
from tensorflow.keras.applications import VGG16
from tensorflow.keras.preprocessing.image import img_to_array
from tensorflow.keras.preprocessing.image import load_img
from tensorflow.keras.applications import imagenet_utils
import numpy as np
import argparse
import imutils

class GradCAM:
    def __init__(self, model, layerName):
        self.model = model
        self.layerName = layerName

    def compute_heatmap(self, image, eps=1e-8):
        # construct our gradient model by supplying (1) the inputs
        # to our pre-trained model, (2) the output of the (presumably)
        # final 4D layer in the network, and (3) the output of the
        # softmax activations from the model
        gradModel = Model(
            inputs=[self.model.inputs],
            outputs=[self.model.get_layer(self.layerName).output,
                     self.model.output])
        with tf.GradientTape() as tape:
            inputs = tf.cast(image, tf.float32)
            (convOutputs, predictions) = gradModel(inputs)
            loss = predictions[:, 0]
        grads = tape.gradient(loss, convOutputs)
        castConvOutputs = tf.cast(convOutputs > 0, "float32")
        castGrads = tf.cast(grads > 0, "float32")
        guidedGrads = castConvOutputs * castGrads * grads
        # the convolution and guided gradients have a batch dimension
        # (which we don't need) so let's grab the volume itself and
        # discard the batch
        convOutputs = convOutputs[0]
        guidedGrads = guidedGrads[0]
        # compute the average of the gradient values, and using them
        # as weights, compute the ponderation of the filters with
        # respect to the weights

        weights = tf.reduce_mean(guidedGrads, axis=(0, 1))
        cam = tf.reduce_sum(tf.multiply(weights, convOutputs), axis=-1)
        # grab the spatial dimensions of the input image and resize
        # the output class activation map to match the input image
        # dimensions
        (w, h) = (image.shape[2], image.shape[1])
        heatmap = cv2.resize(cam.numpy(), (w, h))
        # normalize the heatmap such that all values lie in the range
        # [0, 1], scale the resulting values to the range [0, 255],
        # and then convert to an unsigned 8-bit integer
        numer = heatmap - np.min(heatmap)
        denom = (heatmap.max() - heatmap.min()) + eps
        heatmap = numer / denom
        heatmap = (heatmap * 255).astype("uint8")
        # return the resulting heatmap to the calling function
        return heatmap

    def overlay_heatmap(self, heatmap, image, alpha=0.5,
                        colormap=cv2.COLORMAP_VIRIDIS):
        # apply the supplied color map to the heatmap and then
        # overlay the heatmap on the input image
        heatmap = cv2.applyColorMap(heatmap, colormap)
        output = cv2.addWeighted(image, alpha, heatmap, 1 - alpha, 0)
        # return a 2-tuple of the color mapped heatmap and the output,
        # overlaid image
        return (heatmap, output)


def class_activation_map(level, model_name, img_folder, image_name):
    model = tf.keras.models.load_model('trainedModels/' + "level" + str(level) + "/" + model_name)
    image = os.path.join(img_folder, image_name)
    orig = cv2.imread(image)

    image = load_img(image, target_size=(64, 64))
    image = img_to_array(image)
    image = np.expand_dims(image, axis=0)
    pred = model.predict(image)
    if pred > 0.5:
        label = "intron"
    else:
        label = "exon"
    cam = GradCAM(model, "conv2d_2")
    heatmap = cam.compute_heatmap(image)
    # resize the resulting heatmap to the original input image dimensions
    # and then overlay heatmap on top of the image
    heatmap = cv2.resize(heatmap, (orig.shape[1], orig.shape[0]))
    (heatmap, output) = cam.overlay_heatmap(heatmap, orig, alpha=0.5)
    # draw the predicted label on the output image
    cv2.putText(output, label, (1, 6), cv2.FONT_HERSHEY_SIMPLEX,
                0.3, (255, 255, 255), 1)
    output = np.vstack([orig, heatmap, output])
    output = imutils.resize(output, height=700)
    folder = img_folder + "/CAM_" + model_name.replace(".h5", "")
    isExist = os.path.exists(folder)
    if not isExist:
        os.makedirs(folder)
    if not cv2.imwrite(os.path.join(folder, image_name), output):
        raise Exception("Could not write image")
    #cv2.imshow("Output", output)
    #cv2.waitKey(0)

if __name__ == '__main__':
    for length in ["300"]:
        for type in ["/exon", "/intron"]:
            for level in [3]:
                folder_dir = "../data/images/train_data_" + length + "_" + length + "_rndm/level" + str(level) + type
                for image in os.listdir(folder_dir):
                    class_activation_map(level, "model_" + length + "_" + length + "_rndm.h5", folder_dir, image)

