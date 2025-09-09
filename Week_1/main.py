import matplotlib.pyplot as plt
import numpy as np
import random

#QUESTION 2:
# The species is Streptococcus pneumoniae strain 5141

#QUESTION 6:
# The real genome and the random genome have many different k-mers, but the real genome has more repeats, while the Fibonacci sequence is very simple and makes only a small number of k-mers

def read_file(file):
    seq = []
    with open(file, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def subword(seq, k):
    subword = set()
    for i in range(len(seq) - k + 1):
        strg = seq[i:i+k]
        subword.add(strg)
    return subword

def random_dna(n):
    bases = "ATCG"
    return "".join(random.choices(bases, k=n))

def fibonacci_word(n):
    a0 = "a"
    curr = "ab"
    while(len(curr) < n):
        next = curr + a0
        a0 = curr
        curr = next

    return curr[:n]

def main():
    seq = read_file("genome_hw1.fa")
    arr, arr1, arr2 = [], [], []
    x = range(2, 31)
    for k in range(2, 31):
        result = subword(seq, k)
        arr.append(len(result))

    seq_dna = random_dna(len(seq))
    for i in range(2, 31):
        result = subword(seq_dna, i)
        arr1.append(len(result))

    fibo_seq = fibonacci_word(len(seq))
    for j in range(2, 31):
        result = subword(fibo_seq, j)
        arr2.append(len(result))


    print(len(seq))
    print(len(seq_dna))

    fig, ax = plt.subplots()
    ax.plot(x, arr, marker='o', label="genome_hw1", color="blue")
    ax.plot(x, arr1, marker='o', label="random_genome", color="green")
    ax.plot(x, arr2, marker='o', label="fibonacci_sequence", color="red")

    ax.set_xlabel("k")
    ax.set_ylabel("Nombre de k-mers unique")
    ax.legend()
    plt.show()

    fig2, ax2 = plt.subplots()
    ax2.plot(x, arr2, marker='o', label="fibonacci_sequence", color="red")
    ax2.set_xlabel("k")
    ax2.set_ylabel("Nombre de k-mers unique")
    ax2.legend()
    plt.show()



if __name__ == "__main__":
    main()  
