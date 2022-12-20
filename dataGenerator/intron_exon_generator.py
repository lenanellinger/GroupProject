import ensembl_rest
import json
import random

f = open("prot-cod_genes.txt", "r")
allGenes = []
for line in f:
    allGenes.append(line)
f.close()

randomGenes = random.sample(allGenes, 100)

exons = []
introns = []
gene_exons = []
gene_introns = []

for gene in randomGenes:
    entry = ensembl_rest.symbol_lookup('human', gene.split('\n', 1)[0],
                                       params={'expand': True})
    splittedCanonical = entry['canonical_transcript'].split('.', 1)[0]
    id = 0
    for i in range(len(entry['Transcript'])):
        if entry['Transcript'][i]['id'] == splittedCanonical:
            id = i
    Canonical = entry['Transcript'][id]
    exon_borders = []
    for exon in Canonical['Exon']:
        exons.append(ensembl_rest.sequence_id(exon['id'])['seq'])
        gene_exons.append(gene)
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
        gene_introns.append(gene)

with open("exon.txt", "w") as e:
    for i in range(len(exons)):
        e.write(">" + gene_exons[i] )
        e.write(exons[i] + "\n")
with open("introns.txt", "w") as i:
    for j in range(len(introns)):
        i.write(">" + gene_introns[j] )
        i.write(introns[j] + "\n")