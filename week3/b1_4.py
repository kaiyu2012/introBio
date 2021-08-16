import numpy as np


# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    c  = 1.0
    for i in range(len(Text)):
        c *= Profile[Text[i]][i]
        if not c:
            break
    
    return c


# Write your ProfileMostProbableKmer() function here.
# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats
def ProfileMostProbableKmer(text, k, profile):
    freqMap = {}
    for i in range(len(text)-k+1):
        p = text[i: i+k]
        if p not in freqMap:
            freqMap[p] = Pr(p, profile)
    # m = max(freqMap.values()) 
    m = max(freqMap, key=freqMap.get)
    return m


def main():


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

    text = "GACAGCCGCTGGTGGATATGTTAACAACCACTGAGGTTATCACCCCAGAAGGAGTACGTAACCCTCGGACTACCTGAATGCTCTGGCATTAGGTGAGAATAACGTGGGCGATCATGCAAATAAAGGGCAGGCCCAAAGCTGCATCCGCAACGTGCTGTGATAGAGTCTGAAGTACAGGGCACGTCAGAATCCGGTCGGGTCCAGTGGTAGAACAAGACGTGTTCGCCTGTCCGTCAGGGCATGCACAAATCGGCCGGGACCGCATGTAGATCTCGACGGCCAATTGCGGTTTAGGGCTAAAATCACGAGCTCGCCGAGATTCGCTGGCAGAGCATCGGCTAGTGACTGCCGTCGTGTTGCAGACAAGTCGGCAAGCTTATGTCTTATGCTAACACCCGAACATTAGCGGAGTCGCTTTATGTACGTGGGGGGGCGCTGCTGTGGGACCTCTCCACCGGCTTAGAATCCACTGTCCCCGCTATCTTGGCTTCGTACGCATTTGTGGAAGGGTCGCTAGTCATAGCCCTGACTGAGATGATCCCGACGAGAAAGGGCCGCCCGGTGCCTACGGGGTGTCAGCGACATGTGACCAAAAGCTGGAAGCACGGCCCAGGCCATTTGCCCACTCGACTCAAGGTAAGACCCATGACCTAGCTTGCTCCTAGCGCGTCTCCCTTGGTGGGTGGCCGTCCGCTGGCCTCAAAACGGATTATTGGAAGGCAGCACTAACAGGCCTCGATGAACCGCCTTTTCTAAGCGTGGTTGGCACCTGCTAAATTCTAGCGCCCGCCTAATTCCGTGTCACACTTTAAGATATGTAATGGCTAATATCAAAGCTTGAACGAGCTGTGCGGTTAATTGCGCCACGAATTAAATATTTGTCAAGCCGACCCCGCCGAATCATTCGCTCCCTTTTCTTCAATGGGCATTATCGGTTATAAGGTGCAACAGGGTCGATGTATTAGCGGCGTGTAGAGGGCTACCTCTCTTTTTTCTCGTG"
    k = 14
    p = {
        'A': [0.197, 0.225, 0.211, 0.268, 0.211, 0.254, 0.324, 0.155, 0.352, 0.239, 0.211, 0.211, 0.183, 0.211],
        'C': [0.197, 0.197, 0.254, 0.296, 0.254, 0.254, 0.169, 0.239, 0.197, 0.169, 0.296, 0.268, 0.225, 0.197],
        'G': [0.338, 0.268, 0.338, 0.169, 0.282, 0.268, 0.296, 0.38, 0.268, 0.352, 0.239, 0.225, 0.282, 0.254],
        'T': [0.268, 0.31, 0.197, 0.268, 0.254, 0.225, 0.211, 0.225, 0.183, 0.239, 0.254, 0.296, 0.31, 0.338]
    }
    print(ProfileMostProbableKmer(text, k, p))
    
    

if __name__ == "__main__":
    main()