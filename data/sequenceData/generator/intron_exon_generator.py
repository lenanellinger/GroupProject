import ensembl_rest
import json
import random
from random import randrange

random.seed(42)




class dataGenerator:

    def __init__(self, proteinList):
        '''
        inits the class and save all genes in an array.
        :param proteinList: the file path where all prots are stored
        '''
        f = open(proteinList)
        self.allGenes = []
        for line in f:
            self.allGenes.append(line)
        f.close()
        self.exonArrays = []
        self.intronArrays = []
        self.normalExons = []
        self.normalIntrons = []

    def createExonIntronsInit(self, sequLen=5000, sequNum=2000):
        '''
        Creates a certain number of sequences that have at least the length of sequLen.
        The generated exons and introns are saved in the object
        :param sequLen: the min length of the sequences
        :param sequNum: the number of sequences returned
        :return:
        '''
        longExons = []
        concatExon = ""
        longIntrons = []
        concatIntron = ""
        genes = []
        while len(longExons) < sequNum:
            exons = []
            introns = []
            randomGenes = random.sample(self.allGenes, 15)
            for gene in randomGenes:
                genes.append(gene)
                self.allGenes.remove(gene)

            for gene in randomGenes:
                print("test")

                try:
                    entry = ensembl_rest.symbol_lookup('human', gene.split('\n', 1)[0],
                                                       params={'expand': True})
                except:
                    continue
                splittedCanonical = entry['canonical_transcript'].split('.', 1)[0]
                id = 0
                for i in range(len(entry['Transcript'])):
                    if entry['Transcript'][i]['id'] == splittedCanonical:
                        id = i
                Canonical = entry['Transcript'][id]
                exon_borders = []
                for exon in Canonical['Exon']:
                    exons.append(ensembl_rest.sequence_id(exon['id'])['seq'])

                    exon_border = (exon['start'], exon['end'])
                    exon_borders.append(exon_border)

                '''
                Get all introns
                '''
                exon_borders.sort(key=lambda x: x[0])
                sequence = ensembl_rest.sequence_id(Canonical['id'])['seq']
                startPos = Canonical['start']
                for t in range(len(exon_borders) - 1):
                    introns.append(sequence[exon_borders[t][1] - startPos:exon_borders[t + 1][0] - startPos])

            for ex in exons:
                if len(ex) >= sequLen:
                    longExons.append(ex)
                    continue
                concatExon += ex
                if len(concatExon) >= sequLen:
                    longExons.append(concatExon)
                    concatExon = ""

            if len(longIntrons) <= 2000:
                for intron in introns:
                    if len(intron) >= sequLen:
                        longIntrons.append(intron)
                        continue
                    concatIntron += intron
                    if len(concatIntron) >= sequLen:
                        longIntrons.append(concatIntron)
                        concatIntron = ""
            print(str(len(longExons)) + " exons are created, " + str(len(longIntrons)) + " Introns are created")
            print(len(genes))

        self.intronArrays.append(
            (longIntrons, sequLen, sequLen, "None"))  # the first sequLen being the real length the second the original
        self.exonArrays.append((longExons, sequLen, sequLen, "None"))
        print(len(genes))

    def Trimmer(self, length, sequences, method="rndm"):
        '''
        Select gets a sequence and trims it to the input length. Three different trimming methods are useable
        rnmd gets a random substring, left the prefix and right the suffix of the sequence
        :param length: The length the sequence should be trimmed to
        :param sequences: the input sequence as string
        :param method: how to trim the sequence
        :return: the trimmed sequence
        '''
        result = []
        for sequ in sequences:
            if len(sequ) <= length:
                continue
            match method:
                case "rndm":
                    cutOffLimit = len(sequ) - length
                    cutOff = randrange(cutOffLimit)
                case "left":
                    cutOff = len(sequ) - length
                case "right":
                    cutOff = 0
            result.append(sequ[cutOff: cutOff + length])
        return result

    def sequenceTrimmer(self, length, source, method="rndm"):
        '''
        Trims all sequences with a source-length (from createExonIntronsInit) to the input length
        Different methods are usable (see Trimmer)
        :param length: The length the sequences should be trimmed to
        :param source: The original source-length (or minimal length)
        :param method: how to trim the sequences
        :return: None, saves the new sequences in the object
        '''
        addingTrimmedE = []
        addingTrimmedI = []
        for i in range(len(self.intronArrays)):
            if self.intronArrays[i][3] == "None" and self.intronArrays[i][2] == source:
                addingTrimmedI.append((self.Trimmer(length, self.intronArrays[i][0], method), length, source, method))
                addingTrimmedE.append((self.Trimmer(length, self.exonArrays[i][0], method), length, source, method))
        self.exonArrays = self.exonArrays + addingTrimmedE
        self.intronArrays = self.intronArrays + addingTrimmedI

    def saveSequ(self, IntEx, IntExString):
        '''
        saves all Exons arrays created within the object or all intron arrays
        :param IntEx: Save introns or exons?
        :param IntExString: Information how to store the file
        :return:
        '''
        if IntEx == "Exon":
            savedObject = self.exonArrays
        else:
            savedObject = self.intronArrays
        for array in savedObject:
            with open("../%s_%s_%s_%s.txt" % (IntExString, str(array[1]), str(array[2]), array[3]), "w") as e:
                for i in range(len(array[0])):
                    e.write("> %s; %s, %s, %s, %s \n" % (i, IntExString, str(array[1]), str(array[2]), array[3]))
                    e.write(array[0][i] + "\n")
                e.close()

def main():
    Generator = dataGenerator("prot-cod_genes.txt")
    for l in [500,1000]:
        Generator.createExonIntronsInit(sequLen=l, sequNum=2000)
        print("Data generation done")
    for l in [500,1000]:
        Generator.sequenceTrimmer(length = l,source = l, method="rndm")
        print("Trimming done")
    Generator.saveSequ("Exon", "Exon")
    Generator.saveSequ("Intron", "Intron")

if __name__ == '__main__':
    main()
