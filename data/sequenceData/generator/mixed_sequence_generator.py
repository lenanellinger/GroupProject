import ensembl_rest
import numpy as np
import random
from random import randrange


class mixedSequenceGenerator:

    def __init__(self, proteinList):
        '''
        saves all genes in an array
        :param proteinList: file path of the genes
        '''
        f = open(proteinList)
        self.allGenes = []
        for line in f:
            self.allGenes.append(line)
        f.close()

    def geneAnnotator(self, gene):
        '''
        gets an name of a gene an fetches it from the ensembl database. Information about the exons and introns are
        gathered and the sequence of the gene is saved in the object together with all exon and intron positions
        :param gene:
        :return: A tuple with [0] the sequence [1] the corresponding exon information [2] a list with information about segments in the sequence
        [3] The length of said segments [4] the start pos of said segments
        '''
        print(gene)
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
        ExIntPos = []
        prior = 2
        segmentLength = 0
        segmentPosition = 0
        for anot in ExonAnnotator:
            if int(anot) != prior:
                ExIntLen.append(segmentLength)
                ExIntType.append(prior)
                ExIntPos.append(segmentPosition)
                segmentLength = 0
            segmentLength += 1
            segmentPosition +=1
            prior = int(anot)
        ExIntLen.append(segmentLength)
        ExIntType.append(prior)
        ExIntPos.append(segmentPosition)

        return (sequence, ExonAnnotator, ExIntType, ExIntLen, ExIntPos)

    def generateSegment(self,left, middle, right, length, currentPos):
        '''
        takes three segment lengths and returns segment a given length with random lengths on the sides
        :param left: length of intron on the left
        :param middle: length of exon (or mulitple exons in middle)
        :param right: length of intron on the right
        :param length: output lenght of the total segment
        :return: tuple containing the idx of all borders
        '''
        if left < right:
            maxThreshold = min(left, length - middle)
            minThreshold = max(0, length - middle - right)
            leftAdd = randrange(minThreshold, maxThreshold + 1)
            rightAdd = length - leftAdd - middle
        else:
            maxThreshold = min(right, length - middle)
            minThreshold = max(0, length - middle -left)
            rightAdd = randrange(minThreshold, maxThreshold)
            leftAdd = length - middle - rightAdd
        return (left - leftAdd + currentPos, left + middle + rightAdd + currentPos)


    def mixedSequences(self, amount, length, minlength, maxlength, exonNum ):
        '''
        generates a input amount of mixedSequences(exon and introns), which all have a certain length,
        The exon is searched in random sampled genes and has a minimal length and a maximal length.
        exonNum must be set to 1 (is not yet implemented!!!!!)
        :param amount: The amount of returned mixed sequences
        :param length: The length of the sequences
        :param minlength: minimal length of the exon
        :param maxlength: maximal length of the exon
        :param exonNum: the number of exons within the mixed sequence (NOT IMPLEMENTED)
        :return:
        '''
        sequences = []
        while len(sequences) < amount:
            rndmGene = self.geneAnnotator(random.sample(self.allGenes, 1)[0])

            if exonNum == 1:
                for segment in range(len(rndmGene[2])):
                    if rndmGene[2][segment] == 1 and segment != len(rndmGene[2])-1:  # this solution is supoptimal
                        segmentLength = sum(rndmGene[3][segment - 1:segment + 2])
                        if segmentLength > length and rndmGene[3][segment] > minlength and rndmGene[3][segment] < maxlength:
                            segmentBorders = self.generateSegment(rndmGene[3][segment-1], rndmGene[3][segment],
                                                                  rndmGene[3][segment +1], length,
                                                                  rndmGene[4][segment-1]-rndmGene[3][segment-1]
                                                                  )
                            sequences.append([rndmGene[0][segmentBorders[0]:segmentBorders[1] +1],
                                        segmentBorders,
                                        (rndmGene[4][segment] - rndmGene[3][segment] - segmentBorders[0], rndmGene[4][segment] - segmentBorders[0]),
                                        rndmGene[3][segment] / length])
            print(str(len(sequences)) + " Sequences created")
        return sequences



def main():
    mixedSequGenerator = mixedSequenceGenerator("prot-cod_genes.txt")
    sequences = mixedSequGenerator.mixedSequences(1, 5000,300, 400, 1)
    with open("1_exon_5000_length.txt", "w") as e:
        for i in sequences:
            e.write("> Exon Position: %s. Exon-Intron-Ratio: %s\n" %(str(i[2]), str(i[3])))
            e.write("%s \n" %(i[0]))

if __name__ == '__main__':
    main()



