import ensembl_rest
import random
import tensorflow as tf
import numpy as np
import sys
import os
import re
from tensorflow.keras.preprocessing.image import img_to_array
from tensorflow.keras.preprocessing.image import load_img
import matplotlib.pyplot as plt
sys.path.insert(1, 'data\\images\\')
import imageCreation
sys.path.insert(2, 'data\\sequenceData\\generator\\')
from mixed_sequence_generator import mixedSequenceGenerator
random.seed(666)
class Windows:
    '''
    an object that contains a set of windows of one certain window size. Also can store the ground truth
    and prediction for each single window
    '''
    def __init__(self, sequence, window_size: int):
        self.sequences = []
        self.size = window_size
        self.groundTruth = []
        for baseIdx in range(len(sequence.sequence) - window_size):
            self.sequences.append(sequence.sequence[baseIdx: baseIdx + window_size])
            self.groundTruth.append(sum(sequence.groundTruth[baseIdx: baseIdx + window_size])/window_size)
        print(sum(self.groundTruth))
        self.predictions = []
        self.predictionsIdx = []

class Sequence:
    '''
    an object that holds a sequence, together with its real exon informations and predicted exon informations
    '''
    def __init__(self, sequence):
        self.sequence = sequence
        self.groundTruth = np.zeros(len(sequence))
        self.prediction = np.zeros(len(sequence))

class SlidingWindow:
    '''
    object that contains exactly one gene, the sequence and all windows
    '''
    def __init__(self, gene: str, geneDebug, debugger: bool = False):
        '''
        defines the sequences used. Can input a sequence name and create the whole gene sequence together with the exon info
        or can take a mixedSequence-tuple
        :param gene: gene Name
        :param geneDebug: the mixed sequence
        :param debugger: weather to use the gene Name or the mixedSequence
        '''
        if not debugger:
            self.gene = gene
            self.windowSets = []
            try:     # symbol_lookup is error prone!
                entry = ensembl_rest.symbol_lookup('human', gene.split('\n', 1)[0],
                                                   params={'expand': True})
                self.sequence = Sequence(ensembl_rest.sequence_id(entry['id'])['seq'])
                self.start = entry['start']
                self.end = entry['end']
                for trans in entry['Transcript']:
                    for exon in trans['Exon']:
                        self.sequence.groundTruth[exon['start'] - self.start: exon['end'] - self.start + 1] += 1.0
                cut_off_left = int(round(len(self.sequence.sequence) * 0.2))
                cut_off_right = int(round(len(self.sequence.sequence)* 0.4))
                self.sequence.sequence = self.sequence.sequence[cut_off_left:cut_off_right]
                self.sequence.groundTruth = self.sequence.groundTruth[cut_off_left:cut_off_right]


            except:
                print("error occured")
                self.sequence = None
        else:
            self.gene = "debugger"
            self.windowSets = []
            self.sequence = Sequence(geneDebug[0])
            self.start = 0
            self.end = len(self.sequence.sequence)
            self.sequence.groundTruth[geneDebug[2][0]: geneDebug[2][1]+1] += 1.0


    def createWindows(self, sizes: [int] = [100, 300, 500]):
        '''
        create all windows
        :param sizes: an array of all sizes of windows
        :return:
        '''
        for s in sizes:
            self.windowSets.append(Windows(self.sequence, s))

    def predict(self, levels: [int], all_models):
        '''
        For each set of window, for each window predict the prob of being an exonic sequence. The prediction
        is saved in the windows object. For each window a figure of each input level is created and predicted
        with a model trained with said level
        :param levels: an array containing all levels on which the windows should be predicted
        :param all_models: a dictionary containing all models
        :return:
        '''

        for windows in self.windowSets:
            for level in levels:
                models_level = all_models[level]
                for modelId, model in models_level.items():
                    currentPrediction = np.zeros(len(windows.sequences))
                    # Generate all figures of a certain level for the given window
                    imageCreation.generate_test_images_direct_input(windows.sequences, level = level, output_name = str(windows.size))
                    for i, imageId in enumerate(os.listdir(str(windows.size))):
                        image = load_img(os.path.join(str(windows.size), imageId), target_size=(64, 64))
                        image = img_to_array(image)
                        image = np.expand_dims(image, axis=0)#
                        idx = int(re.findall(r'\d+', imageId)[0])
                        print(idx)
                        print(imageId)
                        print(os.path.join(str(windows.size), imageId))
                        currentPrediction[idx] = 1- model.predict(image)

                    windows.predictions.append(currentPrediction)
                    windows.predictionsIdx.append("%s_%s" %(str(level), modelId))

    def plotResults(self):
        '''
        visualize the results of all Models used and all windows used
        :return:
        '''
        for wind in self.windowSets:
            i = 0

            for model in wind.predictions:
                plt.plot(range(len(wind.sequences)), wind.groundTruth, label=str(wind.size) + " ground Truth")
                plt.plot(range(len(wind.sequences)), model, label = "%s_%s" %(str(wind.size), str(i)))
                plt.title(wind.predictionsIdx[i] +"_"+ str(wind.size))
                plt.ylabel("Probability of being an exon")
                plt.xlabel("Window along the sequence")
                plt.legend()
                plt.show()
                plt.savefig(wind.predictionsIdx[i]+ str(wind.size)+'.png')
                i += 1

if __name__ == '__main__':
    # Import all models
    all_models = {4: {"100_rndm_": tf.keras.models.load_model("model\\trainedModels\\level4\\model_100_100_rndm.h5"),
                    "300_rndm_": tf.keras.models.load_model("model\\trainedModels\\level4\\model_300_300_rndm.h5"),
                    "500_rndm_": tf.keras.models.load_model("model\\trainedModels\\level4\\model_500_500_rndm.h5"),
                    "1000_rndm_": tf.keras.models.load_model("model\\trainedModels\\level4\\model_1000_1000_rndm.h5"),
                    "5000_rndm_": tf.keras.models.load_model("model\\trainedModels\\level4\\model_5000_5000_rndm.h5"),
                    # "intronexon" : tf.keras.models.load_model("model\\trainedModels\\IntronExonClassifier.h5"),
                }, 3: {
                    "100_rndm_3": tf.keras.models.load_model("model\\trainedModels\\level3\\model_100_100_rndm.h5"),
                    "300_rndm_3": tf.keras.models.load_model("model\\trainedModels\\level3\\model_300_300_rndm.h5"),
                    "500_rndm_3" : tf.keras.models.load_model("model\\trainedModels\\level3\\model_500_500_rndm.h5"),
                    "1000_rndm_3" : tf.keras.models.load_model("model\\trainedModels\\level3\\model_1000_1000_rndm.h5"),
                    "5000_rndm_3" : tf.keras.models.load_model("model\\trainedModels\\level3\\model_5000_5000_rndm.h5")
                      }}
    # Generate exactly one mixed sequence
    mixedSequence = mixedSequenceGenerator("data\\sequenceData\\generator\\prot-cod_genes.txt")
    random.seed(42)
    sequences = mixedSequence.mixedSequences(amount=1,length=2000,minlength=500, maxlength=600 , exonNum=1)
    sequence = sequences[0]
    # init the slider with the mixed sequence as input
    slider = SlidingWindow(debugger=True, geneDebug= sequence, gene= "TP53\n")
    # let the slider slide over the sequence to generate windows
    slider.createWindows(sizes=[100, 300, 500])
    # predict the sequences within the windows
    slider.predict(levels=[3], all_models=all_models)
    print(slider.gene)
    # visualize the result
    slider.plotResults()