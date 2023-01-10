import ensembl_rest
import json
import random
from random import randrange

random.seed(42)




class dataGenerator:

    def __init__(self, proteinList):
        f = open(proteinList)
        self.allGenes = []
        for line in f:
            self.allGenes.append(line)
        f.close()
        self.exonArrays = []
        self.intronArrays = []

    def createExonIntronsInit(self, sequLen=5000, sequNum=2000):
        longExons = []
        concatExon = ""
        longIntrons = []
        concatIntron = ""
        while len(longExons) < sequNum:
            exons = []
            introns = []
            randomGenes = random.sample(self.allGenes, 5)
            for gene in randomGenes:
                self.allGenes.remove(gene)

            for gene in randomGenes:

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

            for intron in introns:
                if len(intron) >= sequLen:
                    longIntrons.append(intron)
                    continue
                concatIntron += intron
                if len(concatIntron) >= sequLen:
                    longIntrons.append(concatIntron)
                    concatIntron = ""
            print(str(len(longExons)) + " exons are created, " + str(len(longIntrons)) + " Introns are created")

        self.intronArrays.append(
            (longIntrons, sequLen, sequLen, "None"))  # the first sequLen being the real length the second the original
        self.exonArrays.append((longExons, sequLen, sequLen, "None"))

    def Trimmer(self, length, sequences, method="rndm"):
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
        addingTrimmedE = []
        addingTrimmedI = []
        for i in range(len(self.intronArrays)):
            if self.intronArrays[i][3] == "None" and self.intronArrays[i][2] == source:
                addingTrimmedI.append((self.Trimmer(length, self.intronArrays[i][0], method), length, source, method))
                addingTrimmedE.append((self.Trimmer(length, self.exonArrays[i][0], method), length, source, method))
        self.exonArrays = self.exonArrays + addingTrimmedE
        self.intronArrays = self.intronArrays + addingTrimmedI

    def saveSequ(self, IntEx, IntExString):
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
    for l in [1000]:
        Generator.createExonIntronsInit(sequLen=l, sequNum=2000)
    for l in [1000]:
        Generator.sequenceTrimmer(length = l,source = l, method="rndm")
    Generator.saveSequ("Exon", "Exon")

if __name__ == '__main__':
    main()
