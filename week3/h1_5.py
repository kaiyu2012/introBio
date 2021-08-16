import numpy as np
import copy

GENES = "ACGT"

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def compute_profile(motifs):
    count = {}
    profile = {}
    k = len(motifs[0])
    for symbol in GENES:
        count[symbol] = np.zeros(k, int).tolist()
        profile[symbol] = np.zeros(k, float).tolist()

    t = len(motifs)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1
    
    for symbol in GENES:
        for i in range(k):
            profile[symbol][i] = round(count[symbol][i]/t, 1)
    
    return profile


# Inputs:
#   profile: dna profile matrix
#   dna_slice: a string of dna
# Output:
#   the best k_mer pattern in the dna string
def greedy_kmer_search(profile, dna_slice: str, k: int) -> str:
    results = {}
    for i in range(len(dna_slice) - k + 1):
        chars = dna_slice[i: i + k]
        score = 1
        for j in range(k):
            score *= profile[chars[j]][i+j]
        results[chars] = score

    return max(results, key=results.get)




# Inputs:
#   dna: a collection of Strings, all of the same length
#   k: the length of the motifs/kmers to find
#   t: the number of Strings in Dna
# Output:
#   A collection bestMotifs of k-mers, one from each string in dna, minimizing Score(bestMotifs) from all motif collections of length t visited by the algorithm    
def greedy_motif_search(dna, k, t):
    p = compute_profile(dna)
    bests = []
    for x in dna:
        bests.append(greedy_kmer_search(p, x, k))

    return bests


def main():

    dna = [
        "GGCGTTCAGGCA",
        "AAGAATCAGTCA",
        "CAAGGAGTTCGC",
        "CACGTCAATCAC",
        "CAATAATATTCG"
    ]
    k = 3
    t = 5
    bests = greedy_motif_search(dna, k, t)
    print(bests)

    file = "dna1.txt"
    f = open(file, "r")
    text = f.read().split('\n')
    k = 12
    t = 25

    bests = greedy_motif_search(text, k, t)
    print(bests)


if __name__ == "__main__":
    main()