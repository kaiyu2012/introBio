import numpy as np
import copy



# Input: a piece of pattern and distance d
# Output: neighbors 
chars = "ACGT"
def neighbors(pattern, d):
    assert(d <= len(pattern))

    if d == 0:
        return [pattern]

    r2 = neighbors(pattern[1:], d-1)
    r = [c + r3 for r3 in r2 for c in chars if c != pattern[0]]

    if (d < len(pattern)):
        r2 = neighbors(pattern[1:], d)
        r += [pattern[0] + r3 for r3 in r2]

    return r



# Input: Integers k and d, followed by a string of Dna.
# Output: all the neighbors of k-mir in a pattern dna.
def PatterNeighbors(Pattern, k, d):
    freqMap = []
    freq = []
    for i in range(len(Pattern)-k+1):
        p = Pattern[i:i+k]
        if p not in freq:
            freq.append(p)
    
    freqMap = copy.deepcopy(freq)
            
    for i in range(len(freq)):
        freqMap += neighbors(freq[i], d)

    # remove duplicated
    return(list(dict.fromkeys(freqMap)))


# Input: Integers k and d, followed by a collection of strings Dna.
# Output: All (k, d)-motifs in Dna.
def MotifEnumeration(Dna, k, d):
    Patterns = []
    freqMap = []    
    
    for x in Dna:
        freqMap.append(PatterNeighbors(x, k, d))

    print(freqMap[0])
 
    for x in freqMap[0]:  
        t = True
        for i in range(1, len(freqMap)):
            if x not in freqMap[i]:
                t = False
                break                
        if t:
            Patterns.append(x)

    return Patterns
    





def main():

    # dna = ["ATTTGGC","TGCCTTA","CGGTATC","GAAAATT"]
    dna = ["ACATTCAAATTCCTCTACACATCAT",
            "TGTTTTGACATCGTCCCATCGTGGA",
            "GAGTTTCCTTTCATCTAACGAGCCA",
            "ACAGGTAGCTCGGTATCGTCCATCG",
            "TCTTCGATATAAGGGTATCTGTAGA",
            "ATTTTACGACAACCTTCTTCCCATA"]
    print(*MotifEnumeration(dna, 5, 1))
    exit(0)


    text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
    k = 5
    p = {
        'A': [0.2, 0.2, 0.3, 0.2, 0.3],
        'C':  [0.4, 0.3, 0.1, 0.5, 0.1],
        'G':  [0.3, 0.3, 0.5, 0.2, 0.4],
        'T':  [0.1, 0.2, 0.1, 0.1, 0.2]
    }

    print(ProfileMostProbableKmer(text, k, p))
    
    x = 0.25 ** 9 * 500 * (1000 - 9 + 1)
    print(x)
    # t = "ACGGGGATTACC"
    # p = {
    #     'A':[0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    #     'C':[0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    #     'G':[0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    #     'T':[0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
    #     }
    
    # print(Pr(t, p))
    
    
    

if __name__ == "__main__":
    main()