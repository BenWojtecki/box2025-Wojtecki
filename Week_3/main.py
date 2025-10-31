from xopen import xopen
from readfa import readfq
import matplotlib.pyplot as plt
import numpy as np
# de Bruijn graph (dBG)

#Question 1: Les informations minimales sont les k-1 mers présent dans les noeuds ainsi que les k-mers 
#sur les arrêtes. Pour faire cela, on peut utiliser un ensemblre pour les k-mers et une liste/dictionnaire pour représenter les arrêtes 

#Question 2:

def DNA_complementary(seq):
    new_dna = ''.join(reversed(seq))
    new_dna = new_dna.upper()
    strg = ""
    for letter in new_dna: 
        if (letter == 'A'): 
            strg += "T"
        elif(letter == 'T'):
            strg += "A"
        elif(letter == 'G'):
            strg += "C"
        elif(letter == 'C'):
            strg += "G"
    return strg


def subword(seq, k):
    subword = set()
    for i in range(len(seq) - k + 1):
        strg = seq[i:i+k]
        subword.add(strg)
    return subword

def canonical_kmer(seq):
    r = DNA_complementary(seq)
    return min(seq, r)

def create_dbg(file_path, k, t):
    counts = {}
    compteur = 0
    with xopen (file_path) as fasta:
        for _, seq, _ in readfq(fasta):
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                can_kmer = canonical_kmer(kmer)
                counts[can_kmer] = counts.get(can_kmer, 0) + 1
                compteur += 1
    kmers = []
    for kmer, c in counts.items():
        if c >= t: 
            kmers.append(kmer)

    vertices = set()
    edges = []
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        vertices.add(prefix)
        vertices.add(suffix)
        edges.append((prefix, suffix))
    
    total_stored = len(kmers)

    print("k-mers stored :", compteur)
    print("k-mers in the graph :", total_stored)
    return vertices, edges

def histogram_frequency(file_path, k):
    freq = [] 
    with xopen(file_path) as fasta:
        for _, seq, _ in readfq(fasta):
                counts = {}
    compteur = 0
    with xopen (file_path) as fasta:
        for _, seq, _ in readfq(fasta):
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                can_kmer = canonical_kmer(kmer)
                counts[can_kmer] = counts.get(can_kmer, 0) + 1
                compteur += 1
    for t in range(1, 40):
        kmers = []
        value = 0
        for kmer, c in counts.items():
            if c >= t: 
                kmers.append(kmer)
                value += 1
        freq.append(value)
    ypoints = np.arange(1, 40, 1)
    plt.plot(ypoints, freq)
    plt.show()

    vertices = set()
    edges = []
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        vertices.add(prefix)
        vertices.add(suffix)
        edges.append((prefix, suffix))
    
    total_stored = len(kmers)

    print("k-mers stored :", compteur)
    print("k-mers in the graph :", total_stored)
    return vertices, edges



def main():
    print("=============Perfect reads, forward strand only=============")
    #vertices_ecoli_genome, edges_ecoli_genome = create_dbg("ecoli_genome_150k.fa", 31, 1)
    #For k = 31 and t = 1, we have k-mers stored 153563, k-mers in the graph 153469

    print("=============Perfect reads, both forward/reverse strands=============")
    #vertices_sample_perfect_forward, edges_sample_perfect_forward = create_dbg("reads-20250930T134338Z-1-001/reads/ecoli_sample_perfect_reads_forward.fasta.gz", 31, 1)
    #For k = 31 and t = 1, we have k-mers stored 2100000, k-mers in the graph 149866

    print("=============0.1% error, both strands=============")
    vertices_sample_perfect, edges_sample_perfect = histogram_frequency("reads-20250930T134338Z-1-001/reads/ecoli_sample_perfect_reads.fasta.gz", 31)
    #For k = 31 and t = 1, we have k-mers stored 2100000, k-mers in the graph 149869

    print("=============1% error, only both strands=============")
    vertices_sample, edges_sample = create_dbg("reads-20250930T134338Z-1-001/reads/ecoli_sample_reads_01.fasta.gz", 31, 1)
    #For k = 31 and t = 1, we have k-mers stored 2100000, k-mers in the graph 702101

    print("=============END=============")

if __name__ == "__main__":
    main()

#Unitigs
