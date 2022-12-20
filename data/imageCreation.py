import numpy as np
import itertools
from matplotlib import pyplot as plt
from Bio import SeqIO

def getSequences(filePath):
    """

    :return: list of sequences
    """
    sequences = []
    with open(filePath) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(record.seq)
    return sequences

def createPercentageMatrix(sequence, level, keywords):
    """
    creates the percentage matrix of the given level

    :param sequence: string
    :param level: int
    :param keywords: keyword substrings
    :return: percentage matrix
    """
    if len(sequence) < level:
        print("Sequence is shorten than level")
        return

    dimension = int(np.sqrt(len(keywords)))
    pMatrix = np.zeros((dimension, dimension))

    count_x = 0
    count_y = 0
    for key in keywords:
        count = sequence.count(key)
        pMatrix[count_y][count_x] = count/len(sequence)

        if count_x == dimension-1:
            count_y += 1
            count_x = 0
        else:
            count_x += 1

    return pMatrix

def convertMatrixToImage(pMatrix, intronExon, index, level):
    """
    Converts Percentage Matrix To Image Of Size 875x656x3

    :param pMatrix: percentage Matrix
    :return: image (matrix size 875x656x3)
    """
    n = int(64/pMatrix.shape[0])
    image = np.kron(pMatrix, np.ones((n, n)))
    path = "level" + str(level) + "/" + str(intronExon) + "/test" + str(index) + ".png"
    plt.imsave(path, image, cmap='cividis')


if __name__ == '__main__':
    sequencesIntron = getSequences("introns.txt")
    sequencesExon = getSequences("exon.txt")

    for level in range(1, 7):
        alphabet = ['A', 'C', 'G', 'T']
        keywords = [''.join(i) for i in itertools.product(alphabet, repeat=level)]
        counter = 0
        for i in keywords:
            print(i + " ", end = "")
            counter += 1
            if counter == 16:
                print()
                counter = 0



        for i, seq in enumerate(sequencesIntron):
            pMatrix = createPercentageMatrix(seq, level, keywords)
            convertMatrixToImage(pMatrix, "intron", i, level)
        for i, seq in enumerate(sequencesExon):
            pMatrix = createPercentageMatrix(seq, level, keywords)
            convertMatrixToImage(pMatrix, "exon", i, level)