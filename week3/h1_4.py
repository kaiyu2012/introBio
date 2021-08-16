import numpy as np
import copy

GENS = "ACGT"

INT_GEN = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T'  
    }

GEN_INT = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}


# input: integer value of a pattern and k of k-mer
# output: a pattern string
def int_to_pattern(p: int, k: int) -> str:
    x = ""
    for _ in range(k):
        x += INT_GEN[np.bitwise_and(p, 3)]
        p = p >> 2    
    return x[::-1]


# input: k of k-mir
# output: all possible patterns of k-mir
def k_mir_patterns(k: int) -> dict:
    x = {}
    for i in range(1<<(2*k)):
        x[i] = int_to_pattern(i, k)
    return x


# Input: A piece of Dna (one String only), pattern
# Output: the minimized distance between pattern(K-mer) and dna slice
def minimized_distance(dna_slice: str, pattern: str):
    k = len(pattern)
    m = copy.deepcopy(k)
    for i in range(len(dna_slice) - k + 1):
        count = sum(1 for a, b in zip(dna_slice[i: i+k], pattern) if a != b)
        if count < m:
            m = count    
    return m


# Input:a collection of strings Dna and a string of k-mir
# Output: the distance of pattern (k-mer)
def distance_dna(dna: list, pattern: str) -> int:
    d = 0
    for i in range(len(dna)):
        d += minimized_distance(dna[i], pattern)    
    return d


# Input: An integer k, followed by a collection of strings Dna.
# Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers. 
# (If there are multiple such strings Pattern, then you may return any one.)
def median_string(k: int, dna: list) -> str:
    # x: get all possible k-mirs
    patterns = k_mir_patterns(k)
    print(patterns)
    results = {}
    for i in patterns:
        results[i] = distance_dna(dna, patterns[i]) 
    k = min(results, key=results.get)
    return patterns[k]
        


def main():

    dna = [
        "GGCCTACATGATTGTTACCGTTTGTCGCGCGCCGCACCAACA",
        "TCCCACTTTACTGTAGGCCGCCTAACTATATGATGGTGAATT",
        "CAGTTCCGCCACGCTGATCGCCTAAGTGAGCACATTGTTAAT",
        "CGGCGCTAGCCGGGGGAGACCACATACACCAGACCTAGCCTA",
        "CCCCCTTGATGGGGGACGGTTTGGCGCCTAACCGCGTTACAG",
        "TCTTTGGCTAGTAGCCTACGGTCGTGCCCTACGGAGAGATAT",
        "GGCCTAAATCTAAGTGCCCATTCTGCCAAGTCGACTAGCGTT",
        "CAGTCCTGCCTAGCCCTAACTTCTGAACCCGCCTCAGATATT",
        "TGGATGCGATGCGAGCTGCGCCTACCAGAGAACCATCCCTAT",
        "AACCACGCATAACCTCCCTGCCTAAAGCACCTTAGTGGGAAC"
    ]
    k = 6
    x = median_string(k, dna)
    print("\n" + "The median string: " + x)


if __name__ == "__main__":
    main()