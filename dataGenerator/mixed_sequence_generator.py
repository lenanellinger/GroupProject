import ensembl_rest
import numpy as np
import json
import random
from random import randrange
import sys
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import os


class mixedSequenceGenerator:

    def __init__(self, proteinList):
        f = open(proteinList)
        self.allGenes = []
        for line in f:
            self.allGenes.append(line)
        f.close()

    def geneAnnotator(self, gene):
        try:
            entry = ensembl_rest.symbol_lookup('human', gene.split('\n', 1)[0],
                                               params={'expand': True})
        except:
            return ("",[],[],[])
        splittedCanonical = entry['canonical_transcript'].split('.', 1)[0]
        id = 0
        for i in range(len(entry['Transcript'])):
            if entry['Transcript'][i]['id'] == splittedCanonical:
                id = i
        Canonical = entry['Transcript'][id]

        sequence = ensembl_rest.sequence_id(Canonical['id'])['seq']
        startPos = Canonical['start']
        ExonAnnotator = np.full((len(sequence)), False)
        for exon in Canonical['Exon']:
            exonStart = exon['start'] - startPos
            exonEnd = exon['end'] - startPos
            ExonAnnotator[exonStart:exonEnd+1] = True
        ExIntLen = []
        ExIntType = []
        prior = 2
        segmentLength = 0
        for anot in ExonAnnotator:
            if int(anot) != prior:
                ExIntLen.append(segmentLength)
                ExIntType.append(prior)
                segmentLength = 0
            segmentLength += 1
            prior = int(anot)
        ExIntLen.append(segmentLength)
        ExIntType.append(prior)

        return (sequence, ExonAnnotator, ExIntLen, ExIntType)

    def generateSegment(self,left, middle, right, length):
        '''
        takes three segment lengths and returns segment a given length with random lengths on the sides
        :param left: length of intron on the left
        :param middle: length of exon (or mulitple exons in middle)
        :param right: length of intron on the right
        :param length: output lenght of the total segment
        :return: tuple containing the idx of all borders
        '''
        

    def mixedSequences(self, amount, length, exonNum ):
        sequences = []
        while len(sequences) < amount:
            rndmGene = self.geneAnnotator(random.sample(self.allGenes, 1))
            if exonNum == 1:
                idxStart = 0
                for segment in range(len(rndmGene[3])):
                    if rndmGene[3][segment] == 1:
                        segmentLength = sum(rndmGene[2][segment - 1:segment + 2])
                        if segmentLength > length:
                            sequences.append((rndmGene[0][idxStart: idxStart + segmentLength]))
                        idxStart += sum(rndmGene[2][segment-1: segment+1])


test = mixedSequenceGenerator("prot-cod_genes.txt")
test2 = test.geneAnnotator("A1CF\n")
print(test2[2])
print(test2[3])




