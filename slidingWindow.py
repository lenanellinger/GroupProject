import ensembl_rest
import json
import random
from random import randrange
import numpy as np
from

random.seed(666)
class Windows:
    '''
    an object that contains the windows of one certain window size
    '''
    def __init__(self, sequence: str, window_size: int):
        self.windows = []
        self.size = window_size
        for baseIdx in range(len(sequence) - window_size):
            self.windows.append(sequence[baseIdx: window_size])

class Sequence:
    '''
    an object that holds a sequence, together with its real exon informations and predicted exon informations
    '''
    def __init__(self, sequence):
        self.sequence = sequence
        self.groundTruth = np.zeros(len(sequence))
        self.prediction = np.zeros(len(sequence))

    def increaseGroundTruth(self, idx):
        self.groundTruth[idx] +=1
    def increasePrediction(self, idx):
        self.prediction[idx] +=1

class SlidingWindow:
    '''
    object that contains exactly one gene, the sequence and all windows
    '''
    def __init__(self, gene: str):
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
                    self.sequence.groundTruth[exon['start'] - self.start: exon['end'] - self.start + 1] += 1

        except:
            self.sequence = None

    def createWindows(self, sizes: [int] = [100, 300, 500]):
        for s in sizes:
            self.windowSets.append(Windows(self.sequence, s))

    def predict(self):
        for windows in self.windowSets:
            idx = 0
            for wind in windows.windows:
                # TODO: Generate picture
                # TODO: predict with model (two models even) if intron or exon
                prediction = True # TODO: get an output for the figure that looks like this
                if prediction:
                    self.sequence.prediction[idx: idx + windows.size + 1] +=1



if __name__ == '__main__':

    # TODO: import the models
    slider = SlidingWindow("TNFRSF18\n")
    slider.createWindows()
    print("hdsafad")