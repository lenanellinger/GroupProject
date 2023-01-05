import numpy as np
import os
import itertools
from matplotlib import pyplot as plt
from Bio import SeqIO

def getSequences(filePath):
    """
    returns a list of sequences from specific file path
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
    :param intronExon: if intron or exon
    :param index: index of image
    :param level
    :return: image
    """
    n = int(64/pMatrix.shape[0])
    image = np.kron(pMatrix, np.ones((n, n)))

    # change colormap:
    # - Greys
    # perceptually uniform colormap (best choice in many applications)
    # - cividis
    # - viridis
    # - plasma
    # - inferno
    # - magma
    color_map = "inferno"
    if intronExon == "mixed":
        path = "mixed/" + color_map + "/level" + str(level) + "/" + str(intronExon)
    else:
        path = color_map + "/level" + str(level) + "/" + str(intronExon)
    file_name = "image" + str(index) + ".png"

    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)

    plt.imsave(os.path.join(path, file_name), image, cmap=color_map)


def generate_train_images():
    sequencesIntron = getSequences("../dataGenerator/introns.txt")
    sequencesExon = getSequences("../dataGenerator/exon.txt")

    if len(sequencesExon) > len(sequencesIntron):
        sequencesExon = sequencesExon[:len(sequencesIntron)]
    else:
        sequencesIntron = sequencesIntron[:len(sequencesExon)]

    for level in range(1, 7):
        alphabet = ['A', 'C', 'G', 'T']
        keywords = [''.join(i) for i in itertools.product(alphabet, repeat=level)]
        counter = 0
        for i in keywords:
            print(i + " ", end="")
            counter += 1
            if counter == 16:
                print()
                counter = 0
        print()

        for i, seq in enumerate(sequencesIntron):
            pMatrix = createPercentageMatrix(seq, level, keywords)
            convertMatrixToImage(pMatrix, "intron", i, level)
        for i, seq in enumerate(sequencesExon):
            pMatrix = createPercentageMatrix(seq, level, keywords)
            convertMatrixToImage(pMatrix, "exon", i, level)


def generate_test_mixed_images():
    mixedSequences = getSequences("../dataGenerator/1_exon_5000_length.txt")
    level = 5

    alphabet = ['A', 'C', 'G', 'T']
    keywords = [''.join(i) for i in itertools.product(alphabet, repeat=level)]

    for i, seq in enumerate(mixedSequences):
        pMatrix = createPercentageMatrix(seq, level, keywords)
        convertMatrixToImage(pMatrix, "mixed", i, level)


if __name__ == '__main__':
    generate_train_images()
    # generate_test_mixed_images()