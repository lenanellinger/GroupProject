import ensembl_rest
import random
import tensorflow as tf
import numpy as np
import sys
import os
import re
import matplotlib.pyplot as plt
sys.path.insert(1, 'data\\images\\')
import imageCreation
sys.path.insert(2, 'data\\sequenceData\\generator\\')
from mixed_sequence_generator import mixedSequenceGenerator
random.seed(666)
class Windows:
    '''
    an object that contains the windows of one certain window size
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

class Sequence:
    '''
    an object that holds a sequence, together with its real exon informations and predicted exon informations
    '''
    def __init__(self, sequence):
        self.sequence = sequence
        self.groundTruth = np.zeros(len(sequence))
        self.prediction = np.zeros(len(sequence))
        self.timesPredicted = np.zeros(len(sequence))  # needed for the edge-correction

    def increaseGroundTruth(self, idx):
        self.groundTruth[idx] +=1.0
    def increasePrediction(self, idx):
        self.prediction[idx] +=1.0

class SlidingWindow:
    '''
    object that contains exactly one gene, the sequence and all windows
    '''
    def __init__(self, gene: str, geneDebug, debugger: bool = False):
        if not debugger:
            self.gene = gene
            self.windowSets = []
            try:
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
        for s in sizes:
            self.windowSets.append(Windows(self.sequence, s))

    def predict(self, levels: [int], all_models):

        for windows in self.windowSets:
            for level in levels:
                models_level = all_models[level]
                for model in models_level.values():
                    currentPrediction = np.zeros(len(windows.sequences))
                    # Generate all figures of a certain level for the given window
                    imageCreation.generate_test_images_direct_input(windows.sequences, level = level, output_name = str(windows.size))
                    for i, image in enumerate(os.listdir(str(windows.size))):
                        im = plt.imread(os.path.join(str(windows.size), image))[:, :, :3]
                        idx = int(re.findall(r'\d+', image)[0])
                        print(idx)
                        print(image)
                        print(os.path.join(str(windows.size), image))
                        currentPrediction[idx] = 1- model.predict(np.expand_dims(im, axis = 0))
                        #currentPrediction[idx] = round(currentPrediction[idx])
                        #self.sequence.prediction[idx : idx + windows.size] += 1 - all_models[level][windows.size - 50](np.expand_dims(im, axis = 0))
                        #self.sequence.timesPredicted[idx : idx + windows.size] += 1.0
                    windows.predictions.append(currentPrediction)

    def adjustEdgesAndGroundTruth(self):
        maxPredictions = np.max(self.sequence.timesPredicted)
        for i in range(len(self.sequence.timesPredicted)):
            self.sequence.prediction[i] *= (maxPredictions / self.sequence.timesPredicted[i])
        self.sequence.groundTruth *= maxPredictions

    def plotResults(self):
        for wind in self.windowSets:
            i = 0
            plt.plot(range(len(wind.sequences)), wind.groundTruth, label=str(wind.size) + " ground Truth")
            for model in wind.predictions:
                plt.plot(range(len(wind.sequences)), model, label = "%s_%s" %(str(wind.size), str(i)))
                i += 1
            plt.legend()
            plt.show()





if __name__ == '__main__':
    # TODO: import the models
    all_models = {4: {#"100_rndm_": tf.keras.models.load_model("model\\trainedModels\\model_100_100_rndm.h5"),
                  #"300_rndm_": tf.keras.models.load_model("model\\trainedModels\\model_300_300_rndm.h5"),
                  #"500_rndm_": tf.keras.models.load_model("model\\trainedModels\\model_500_500_rndm.h5"),
                        #"1000_rndm_": tf.keras.models.load_model("model\\trainedModels\\level4\\model_1000_1000_rndm.h5"),
                     # "5000_rndm_": tf.keras.models.load_model("model\\trainedModels\\level4\\model_5000_5000_rndm.h5"),
                      #"100_no_": tf.keras.models.load_model("model\\trainedModels\\model_100_100_None.h5."),
                      #"300_no_": tf.keras.models.load_model("model\\trainedModels\\model_300_300_None.h5"),
                      #"500_no_": tf.keras.models.load_model("model\\trainedModels\\model_500_500_None.h5"),
                      #"1000_no_": tf.keras.models.load_model("model\\trainedModels\\level4\\model_1000_1000_None.h5"),
                      #"5000_no_": tf.keras.models.load_model("model\\trainedModels\\level4\\model_5000_5000_None.h5"),
                       # "intronexon" : tf.keras.models.load_model("model\\trainedModels\\IntronExonClassifier.h5"),
                        "100_rndm_3": tf.keras.models.load_model("model\\trainedModels\\level3\\model_100_100_rndm.h5"),
        "300_rndm_3": tf.keras.models.load_model("model\\trainedModels\\level3\\model_300_300_rndm.h5"),
        "500_rndm_3" : tf.keras.models.load_model("model\\trainedModels\\level3\\model_500_500_rndm.h5"),
        "1000_rndm_3" : tf.keras.models.load_model("model\\trainedModels\\level3\\model_1000_1000_rndm.h5"),
        "5000_rndm_3" : tf.keras.models.load_model("model\\trainedModels\\level3\\model_5000_5000_rndm.h5")
                      }}
    print("test 1: done")
    mixedSequence = mixedSequenceGenerator("data\\sequenceData\\generator\\prot-cod_genes.txt")
    random.seed(42)
    sequences = mixedSequence.mixedSequences(amount=1,length=2000,minlength=550, maxlength=600 , exonNum=1)
    sequence = sequences[0]

    print("test 2: done")
    slider = SlidingWindow(debugger=True, geneDebug= sequence, gene= "TP53\n")
    slider.createWindows(sizes=[100, 300, 500])
    slider.predict(levels=[3], all_models=all_models)
   # slider.adjustEdgesAndGroundTruth()
    slider.plotResults()