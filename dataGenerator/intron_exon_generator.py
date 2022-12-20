import ensembl_rest
import json
import random

f = open("prot-cod_genes.txt", "r")
allGenes = []
for line in f:
    allGenes.append(line)
f.close()

random.seed(42)
exons = []
longExons = []
concatExon = ""
introns = []
longIntrons = []
concatIntron = ""
while len(longExons) < 50:
    randomGenes = random.sample(allGenes, 15)
    idx = []
    for gene in randomGenes:
        allGenes.remove(gene)
    idx.sort(reverse=True)



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
        if len(ex) >= 5000:
            longExons.append(ex)
            continue
        concatExon += ex
        if len(concatExon) >= 5000:
            longExons.append(concatExon)
            concatExon = ""

    for intron in introns:
        if len(intron) >= 5000:
            longIntrons.append(intron)
            continue
        concatIntron += intron
        if len(concatIntron) >= 5000:
            longIntrons.append(concatIntron)
            concatIntron = ""
for exon in longExons:
    print(len(exon))

with open("exon.txt", "w") as e:
    for i in range(len(longExons)):
        e.write(">" + str(i) + "\n")
        e.write(longExons[i] + "\n")
with open("introns.txt", "w") as i:
    for j in range(len(longIntrons)):
        i.write(">" + str(j) + "\n")
        i.write(longIntrons[j] + "\n")