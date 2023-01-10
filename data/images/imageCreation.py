import numpy as np
import os
import itertools
from matplotlib import pyplot as plt
from Bio import SeqIO


def get_sequences(file_path):
    """
    returns a list of sequences from specific file path
    :return: list of sequences
    """
    sequences = []
    with open(file_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(record.seq)
    return sequences


def create_percentage_matrix(sequence, level, keywords):
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
    p_matrix = np.zeros((dimension, dimension))

    count_x = 0
    count_y = 0
    for key in keywords:
        count = sequence.count(key)
        p_matrix[count_y][count_x] = count / len(sequence)

        if count_x == dimension - 1:
            count_y += 1
            count_x = 0
        else:
            count_x += 1

    return p_matrix


def convert_matrix_to_image_and_save(p_matrix, intron_exon, index, level, folder_name):
    """
    Converts Percentage Matrix To Image Of Size 64x64x3
    and saves image

    :param p_matrix: percentage Matrix
    :param intron_exon: if intron or exon or mixed_inferno_level5
    :param index: index of image
    :param level
    :param folder_name: in which folder the images should be saved, if empty folder name is created
    :return: image
    """
    n = int(64 / p_matrix.shape[0])
    image = np.kron(p_matrix, np.ones((n, n)))

    color_map = "inferno"
    if folder_name == "":
        if intron_exon != "intron" and intron_exon != "exon":
            path = intron_exon + "_" + color_map + "_level" + str(level) + "/" + str(intron_exon)
        else:
            path = color_map + "/level" + str(level) + "/" + str(intron_exon)
    else:
        path = folder_name + "/" + str(intron_exon)
    file_name = "image" + str(index) + ".png"

    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)

    plt.imsave(os.path.join(path, file_name), image, cmap=color_map)


def generate_images_for_level(level, sequences_intron, sequences_exon, folder_name=""):
    """
    generate images for a given level
    :param level
    :param sequences_intron: list of intron sequences
    :param sequences_exon: list of exon sequences
    :param folder_name: name of folder to be saved, if empty folder name is created
    """
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

    for i, seq in enumerate(sequences_intron):
        pMatrix = create_percentage_matrix(seq, level, keywords)
        convert_matrix_to_image_and_save(pMatrix, "intron", i, level, folder_name)
    for i, seq in enumerate(sequences_exon):
        pMatrix = create_percentage_matrix(seq, level, keywords)
        convert_matrix_to_image_and_save(pMatrix, "exon", i, level, folder_name)


def generate_train_images(introns_file, exons_file, level=None):
    """
    creates images for given intron and exon sequences
    :param level: if None all levels are created
    :param introns_file: fasta format with intron sequences
    :param exons_file: fasta format with exon sequences
    """
    sequences_intron = get_sequences("../sequenceData/" + introns_file)
    sequences_exon = get_sequences("../sequenceData/" + exons_file)

    if len(sequences_exon) > len(sequences_intron):
        sequences_exon = sequences_exon[:len(sequences_intron)]
    else:
        sequences_intron = sequences_intron[:len(sequences_exon)]

    if level is not None:
        mode = exons_file.replace("Exon_", "").replace(".txt", "")
        generate_images_for_level(level, sequences_intron, sequences_exon, "train_data_" + mode)
    else:
        # generate all levels
        for level in range(1, 7):
            generate_images_for_level(level, sequences_intron, sequences_exon)


def generate_mixed_images():
    """
    generates mixed sequences (an intron with an exon contained)
    """
    mixed_sequences = get_sequences("../sequenceData/1_exon_5000_length.txt")
    level = 5

    alphabet = ['A', 'C', 'G', 'T']
    keywords = [''.join(i) for i in itertools.product(alphabet, repeat=level)]

    for i, seq in enumerate(mixed_sequences):
        p_matrix = create_percentage_matrix(seq, level, keywords)
        convert_matrix_to_image_and_save(p_matrix, "mixed", i, level, folder_name="")


if __name__ == '__main__':
    for length in ["100", "300", "500", "1000"]:
        for trim in ["rndm", "None"]:
            introns = "Intron_" + length + "_" + length + "_" + trim + ".txt"
            exons = "Exon_" + length + "_" + length + "_" + trim + ".txt"
            generate_train_images(introns, exons, level=4)
