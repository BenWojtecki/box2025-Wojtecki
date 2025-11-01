from xopen import xopen
from readfa import readfq
import matplotlib.pyplot as plt
import numpy as np

files = [
        "ecoli_genome_150k.fa",
        "reads-20250930T134338Z-1-001/reads/ecoli_sample_perfect_reads_forward.fasta.gz",
        "reads-20250930T134338Z-1-001/reads/ecoli_sample_perfect_reads.fasta.gz",
        "reads-20250930T134338Z-1-001/reads/ecoli_sample_reads_01.fasta.gz",
    ]


# de Bruijn graph (dBG)

#Question 1: Les informations minimales sont les k mers présent dans les noeuds ainsi que les k-1 mers
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
    with xopen (file_path) as fasta:
        for _, seq, _ in readfq(fasta):
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                can_kmer = canonical_kmer(kmer)
                counts[can_kmer] = counts.get(can_kmer, 0) + 1

    kmers = []
    for kmer, c in counts.items():
        if c >= t: 
            kmers.append(kmer)

    edges = set()
    prefix_dict = {}
    for kmer in kmers:
        prefix = kmer[:-1]
        prefix_dict.setdefault(prefix, []).append(kmer)

    for kmer in kmers:
        suffix = kmer[1:]
        for next_kmer in prefix_dict.get(suffix, []):
            edges.add((kmer, next_kmer))

                        
    print("k-mers stored :", sum(counts.values()))
    print("k-mers in the graph :", len(kmers))
    print("edges in the graph :", len(edges))
    return kmers, edges

def histogram_frequency(file_path, k):
    counts = {}
    with xopen (file_path) as fasta:
        for _, seq, _ in readfq(fasta):
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                can_kmer = canonical_kmer(kmer)
                counts[can_kmer] = counts.get(can_kmer, 0) + 1
    
    max_freq = max(counts.values())
    freq = []   

    for t in range(1, 40):
        compteur = 0
        for c in counts.values():
            if c >= t:
                compteur += 1   
        freq.append(compteur)

    plt.plot(range(1, 40), freq)
    plt.xlabel("Abundance threshold (t)")
    plt.ylabel("Number of kmers >= t")
    plt.title(f"K-mer frequency histogram (k={k})")
    plt.show()
    return counts, freq

#Unitigs

def build_adjacency(kmers):
    adj = {}
    indeg = {}
    prefix_map = {}

    for kmer in kmers:
        indeg[kmer] = 0
        adj[kmer] = []
        prefix = kmer[:-1]
        if prefix not in prefix_map:
            prefix_map[prefix] = []
        prefix_map[prefix].append(kmer)

    for kmer in kmers:
        suffix = kmer[1:]
        if suffix in prefix_map:
            for next_kmer in prefix_map[suffix]:
                adj[kmer].append(next_kmer)
                indeg[next_kmer] += 1

    return adj, indeg


def unitig_from(kmers, edges, k, start_kmer):
    adj, indeg = build_adjacency(kmers)
    unitig = start_kmer
    right_kmer = start_kmer

    while right_kmer in adj and len(adj[right_kmer]) == 1:
        next_kmer = adj[right_kmer][0]
        incoming_edges = 0
        for a, b in edges:
            if b == next_kmer:
                incoming_edges += 1
        
        if incoming_edges == 1:
            unitig += next_kmer[-1]
            right_kmer = next_kmer
        else:
            break
    
    left_kmer = start_kmer
    while True:
        anterior = []
        for a, b in edges:
            if b == left_kmer:
                anterior.append(a)
        
        if len(anterior) == 1:
            anterior_kmer = anterior[0]
            if len(adj[anterior_kmer]) == 1:
                unitig = anterior_kmer[0] + unitig
                left_kmer = anterior_kmer
            else:
                break
        else:
            break

    return unitig


def create_unitigs(kmers, edges, k):
    adj, indeg = build_adjacency(kmers)
    visited = set()
    unitigs = []

    for kmer in kmers:
        if (indeg[kmer] != 1 or len(adj[kmer]) != 1) and kmer not in visited:
            seq = kmer
            current = kmer
            visited.add(current)

            while current in adj and len(adj[current]) == 1:
                next_kmer = adj[current][0]
                if indeg[next_kmer] == 1:
                    seq += next_kmer[-1]
                    current = next_kmer
                    visited.add(current)
                else:
                    break

            unitigs.append(seq)

    return unitigs

def main():

    print("=============Question 3=============")
    create_dbg("ecoli_genome_150k.fa", 31, 1)
    #For k = 31 and t = 1, we have k-mers stored 153563 also there is 153469 nodes in the graph and 92796 edges
    print("=============Question 4=============\n")

    print("=============Perfect reads, forward strand only=============")
    create_dbg("ecoli_genome_150k.fa", 31, 1)
    #For k = 31 and t = 1, we have k-mers stored 153563 also there is 153469 nodes in the graph and 92796 edges

    print("=============Perfect reads, both forward/reverse strands=============")
    create_dbg("reads-20250930T134338Z-1-001/reads/ecoli_sample_perfect_reads_forward.fasta.gz", 31, 1)
    #For k = 31 and t = 1, we have k-mers stored 2100000, also there is 149869 nodes in the graph and 90618 edges

    print("=============0.1% error, both strands=============")
    create_dbg("reads-20250930T134338Z-1-001/reads/ecoli_sample_perfect_reads.fasta.gz", 31, 1)
    #For k = 31 and t = 1, we have k-mers stored 2100000, also there is 149869 nodes in the graph and 90621 edges

    print("=============1% error, only both strands=============")
    create_dbg("reads-20250930T134338Z-1-001/reads/ecoli_sample_reads_01.fasta.gz", 31, 1)
    #For k = 31 and t = 1, we have k-mers stored 2100000, also there is 213031 nodes in the graph and 129297 edges

    print("============= Question 5 =============")
    print("Histogramme de fréquences pour différents seuils d’abondance t")
    histogram_frequency("reads-20250930T134338Z-1-001/reads/ecoli_sample_perfect_reads.fasta.gz", 31)
    #optimal t is aroung 20 for reads with 0.1% error, both strands

    histogram_frequency("reads-20250930T134338Z-1-001/reads/ecoli_sample_reads_01.fasta.gz", 31)
    #optimal t is aroung 5 for reads with 0.1% error, only both strands

    print("============= Question 7=============")
    target_kmer = "CGCTCTGTGTGACAAGCCGGAAACCGCCCAG"
    for file in files:
        print(f"\n Dataset: {file}")
        print("--------------------------------------------------")
        for t in [1, 2, 3, 4]:
            kmers, edges = create_dbg(file, 31, t)
            if target_kmer in kmers:
                u = unitig_from(kmers, edges, 31, target_kmer)
                print(f"Length of the unitig containing the target k-mer : {len(u)}")
            else:
                print("Target k-mer not present in the graph")
        print("--------------------------------------------------")
    #For the E. coli genome, the unitig containing the target k-mer only exists at t = 1 (33 bp) and disappears as t increases, because low-abundance k-mers are filtered out.
    #For the read datasets, the unitig size remains constant at 33 bp, since the high coverage preserves the true genomic k-mers across all abundance thresholds.
    print("=============END=============")

if __name__ == "__main__":
    main()

#Unitigs
