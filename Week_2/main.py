from readfa import readfq
import matplotlib.pyplot as plt

files = ["Files/file1.fa", "Files/file2.fa", "Files/file3.fa", "Files/file4.fa", "Files/file5.fa", "Files/file6.fa"]
base_map = {"A" : 0, "C": 1, "G" : 2, "T" : 3}

def subword(seq, k):
    subword = set()
    for i in range(len(seq) - k + 1):
        strg = seq[i:i+k]
        subword.add(canonical_kmer(strg))
    return subword


def count_kmers(seq, k):
    counts_kmer = {}
    for i in range(len(seq) - k + 1):
        strg = seq[i:i+k]
        canonical = canonical_kmer(strg)
        if canonical in counts_kmer:
            counts_kmer[canonical] += 1
        else:
            counts_kmer[canonical] = 1
    return counts_kmer

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

def canonical_kmer(seq):
    r = DNA_complementary(seq)
    return min(seq, r)


def jaccard_index(seq1, seq2):
    union_seq = seq1.union(seq2)
    intersect_seq = seq1.intersection(seq2)
    if not union_seq:
        return 0 
    return (len(intersect_seq) / len(union_seq))

def kmer_to_int(kmer):
    num = 0
    for letter in kmer:
        num = num * 4 + base_map[letter]
    return num
 

def plot_kmer_histogram(file, k):
    counts = {}
    with open(file) as f:
        for _, seq, _ in readfq(f):
            seq_counts = count_kmers(seq, k)
            for kmer, c in seq_counts.items():
                if kmer in counts:
                    counts[kmer] += c
                else:
                    counts[kmer] = c

    frequencies = list(counts.values())
    plt.figure(figsize=(8, 5))
    plt.hist(frequencies, bins = 50, edgecolor="black")
    plt.xlabel("k-mer count")
    plt.ylabel("Number of distinct k-mers")
    plt.tight_layout()
    plt.show()


def main():
    kmers_file = {}
    for file in files:
        kmers = set()
        count = 0
        cumulative_count = 0

        with open(file) as f:
            for name, seq, _ in readfq(f):
                count += 1
                cumulative_count += len(seq)
                seq_counts = count_kmers(seq, 20)
                
                result = subword(seq, 20)
                kmers.update(result)

                canonical_seq = DNA_complementary(seq)
                result_canonical = subword(canonical_seq, 20)
                kmers.update(result_canonical)


        kmers_file[file] = kmers
        print(f"{file}: {count} sequences, {cumulative_count} bases, {len(kmers)} distinct canonical 20-mers")

    for i in range(len(files)):
        for j in range(i + 1, len(files)):
            file1, file2 = files[i], files[j]
            ji = jaccard_index(kmers_file[file1], kmers_file[file2])
            print(f"Jaccard index between {file1} and {file2}: {ji:.3f}")
            
#File4 and File6 have a jaccard index near to 1, like File2 and File5. We can say that the sequences in theses files correspond to almost the same genome. 
#This suggest that these theses genomes come from the same species or strain  
    plot_kmer_histogram("Files/file1.fa", 20)

if __name__ == "__main__":
    main()