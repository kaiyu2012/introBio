import numpy as np

# Input:  Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array


# Reproduce the PatternCount function here.
def PatternCount(Symbol, ExtendedGenome):
    # type your code here
    count = 0
    for i in range(0,len(ExtendedGenome)-len(Symbol)+1):
        if ExtendedGenome[i:i+len(Symbol)] == Symbol:
            count = count+1
    return count



def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

# Input:  A String Genome
# Output: The skew array of Genome as a list.
def SkewArray(Genome):
    skew = [0]
    score = {"A":0, "T":0, "C":-1, "G":1}
    for i in range(0,len(Genome)):
            skew.append(score[Genome[i]] + skew[i])
    return skew


# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Gnome):
    # position = []
    skew = np.array(SkewArray(Gnome))
    position = np.argwhere(skew == np.amin(skew)).ravel()
 
    return position

# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q) -> int:
    # your code here
    l = [len(p), len(q)]
    count = 0
    for i in range(min(l)):
        if p[i] != q[i]:
            count += 1

    count += abs(l[0] - l[1])
    return count
        

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    
    lt = len(Text)
    lp = len(Pattern)
    for i in range(lt-lp+1):
        count = 0
        
        for j in range(lp):
            if Text[i+j] != Pattern[j]:
                count += 1
            if count > d:
                break
        
        if count <= d:
            positions.append(i)
              
    return positions



# Input: String text, k-mer, and d
# out put: A list of k-mer Patterns which has Hamming distance at most d from Pattern
def FrequentWordsWithMismatches(Text, k, d):
    l = len(Text)
    patterns = {}
    for i in range(l-k+1):
        x = Text[i:i+k]
        if x not in patterns:
            patterns[x] = 1
        else:
            patterns[x] += 1
      
    freq = []  
    for p in patterns:
        if patterns[p] == d:
            freq.append(p)
    
    print(freq)
            

def immediateneighbor(pattern):
    neighborhood=[pattern]
    for i in range(len(pattern)):
        symbol=pattern[i]
        for nucleotide in 'ATGC':
            if nucleotide!=symbol:
                neighbor=pattern[:i]+nucleotide+pattern[i+1:]
                neighborhood.append(neighbor)
    return neighborhood




def FrequentWordsWithMismatches_Sys(Text, k, d):
    Patterns = []
    freqMap = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i: i+k]
        neighborhood = Neighbors(Pattern, d)
        for j in range(len(neighborhood)):
            neighbor = neighborhood[j]
            if not freqMap[neighbor]:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] += 1
                
    m = max(freqMap)
    for Pattern in freqMap:
        if freqMap[Pattern] == m:
            Patterns.append(Pattern)
    return Patterns


def Neighbors(Pattern, d):
        if d == 0:
            return {Pattern}
        if len(Pattern) == 1: 
            return {'A', 'C', 'G', 'T'}
        Neighborhood = {}
        SuffixNeighbors = Neighbors(Suffix(Pattern), d)
        # for each string Text from SuffixNeighbors
        #     if HammingDistance(Suffix(Pattern), Text) < d
        #         for each nucleotide x
        #             add x • Text to Neighborhood
        #     else
        #         add FirstSymbol(Pattern) • Text to Neighborhood
        # return Neighborhood
    

def frequent_words_mismatches(seq, k, d):
    from itertools import product
    bases = ['A', 'C', 'G', 'T']
    bases_comb = [''.join(base) for base in product(bases, repeat = k)]
    freq_words = {}

    for pattern in bases_comb:
        count = ApproximatePatternMatching(pattern, seq, d)
        if count not in freq_words:
            freq_words[count] = [pattern]
        else:
            freq_words[count].append(pattern)

    return freq_words[max(freq_words)]     
        
  

def main():
    
    # t1 = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    # x = frequent_words_mismatches(t1, 4, 1)
    # print(x)
    
    # exit(0)
    t1 = "CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT"
    t2 = "CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG"
    
    print(HammingDistance(t1, t2))
    
    t1 = "GCATACACTTCCCAGTAGGTACTG"
    print(SkewArray(t1))
    
    t1 = "CGTGACAGTGTATGGGCATCTTT"
    t2 = "TGT"
    print(len(ApproximatePatternMatching(t1, t2, 1)))
    exit(0)
    
    # Hidden 1.4
    
    text = "GGATGTACCGAAGATGGCAGGCGAGATGTTACATGGTGTGATCGGGGAGATAATGCAAAGTGCAGTCCGAGAAGAGTCTTCAACATAATTTGCAGGCTATGTATACGGTATTGATATAAGTTCGTCTTGGGGCGGTAGTGAAGTCGGTGTTGGCCACTCCTGAGATAAGAATTGTACGTGTGTCTGGTCGCTTATTGGTAGGGCCGTGCCGAAGAAACTGTCGGTTCTCATCCTCGTCTGCGTTGATCCTACCCCTGCCAGATTGGCTCTCTATGGTAGAGGCCCTGCACTGTTGAATCATTTCGGTGGCTACAGGGCACTTCAAAGACGTCCTAAAGTCCCTTGAAGCACCT"
    pattern = "TCCTCGT"
    d = 2
    print(len(ApproximatePatternMatching(text, pattern, d)))
    
    exit(0)
    
    
    # Hidden 1.4
    f = open("dataset_9_4.txt", "r")
    text = f.read().split("\n")

    print(*ApproximatePatternMatching(text[1], text[0], int(text[2])))
    
    
    
    # Hidden 1.4
    text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
    pattern = "ATTCTGGA"
    d = 3
    print(ApproximatePatternMatching(text, pattern, d))
    exit(0)

    
    # Hidden 1.4
    p = "GAAACCTTATTGGGATGTCTGGTGTCAGCTTAAAACGCCGAGCTTACCCGAATAGAGAGAGGGTCACCTAAGATTCGTGTGAACTCCATTGGATCGCACTCTCGCGACATGAGGACCCGGCCGTTCCGCCCGCCCGCTCGGACTGGTATAAGCGGCGCCTGGCCCCGTCTTCACAGGCATTCCATCCCATGATCCCTATTAGACCAAGCCTAGTTTCATGCTAGCATATCACGTGGGCTGCTGAGAAATGTGGTGATCACCGAGCCCTTCAAGCGGTCTGGTACCACTTGGAATGCACCTCGCACCAGCGATAAGCCCATGCCCTCTGCGTGGCATTGTAACAGTAGAGACGGGCGCGACCAGAGGAACTGGGACTTGAGGCCGCAATGAGCAGTACCCCTGATTGCATCTGGCCGAAAGCCAAAAGCTCTGTTCTTTCTGCATCAGGTCTGTATTGTAGTCCAGCAGTTGCCATACATCTTTAAGTATCATCTATTACCCCCTTGCCTCATTACCGTAACTTTAAGCTCCAGGAGGCTCTGGGAGCACACGAGCTCTAAATATCGCCTGGAGGTGTTGCGTACGCAGTCGAGGGCATGGCCGGTCCCGAAATTAGTGATAAACTCTTCCAATTGTCCGACTTATCCAAGCGGAGGTAGCAATACGACTAACCTACGGCAAGCGAGGGCCTCAACGGTGAGAGACTATATCGTTCTCCCCCGGGTGACATCGTCAACAGGGCGCAGTTTCCGTGAGAATCGTAAAGACAGAACCACCTGCTGTTCAGCAAGAGTTTTCTCAACAGCGAAAAATAGATGTGAGTAAGCCGTATTGTCTCTACCTAGCAAAGTAGGCAGCGTGATGCCGCGGGTGGACAGAACTCCCCATGTTCTTTAGGCATCTTCGAAAAAAGTGGCTAGCCACATGCGCCTAAGCAAGGGTACATGTATCCAATGACGTAGGAATAATGACAGGTCCTAATTTGAGTGGTATGTGACACTCGCACAGTCTATATATTCCCGACTTCTTCTTTGGATAGAATCGGTACTGACTGCA"
    q = "ACATAATAATTTGATTGAGCGGCTCAATGCCGTCGATGATCCAATGGTCTACCATACATCCACGTACTGGAATTGTCGTTAATGACATTGGGCCATCTAAGCAGACGCTGCCCGGGTTCTTGTGCCGTTAATAGCTCATCATTACTGGGAACAGCATCCTGAGCAAGATTTTAGGACGTGTGGATAACCATAGCAGCGTTGGCTCAGGGTTCCGATTGAAGAGACCAAGGTAGCTCCACTCCATGATCTAGCGTAACTATGTCCCGCCGTAATAGGTTCCTGGGGGCCTTCCGTATCTTGTACCCAGAAGGGCGCATTGTGAGTCGGTGATGACCACCGGGGCACCGCGTTGACTCCGAAAATCACCTGTATACTCGATTAGGAGACGTACGGGGATGGTCTAGGTGAGCGCGGAACCCGGTCAAATTCACCAACTACTATCCCCACAACCTTCTAAGTTTTCACACTCGGAACCAACGTGGTGGCGCCCTAAACTGCGCTCTTTACGAGGCAGGGCTATGATAAGGCCCCCCGGTAGCATATATAAGGACTAGTTTCGGAGCCGATCATTGGGCAACAGTATGGCCTCACCTCAATACTCATTCATGTGACCCTTTGGCCATCGCTACTAACGGAAAGATATGCAAAAGCCAGGATGCCTACGGGCCTGATCCGTCCTAAATGGTTCTTCTTACAGTGGTTTCCAAAGATTTGAGTAGACGATATCTTTTAGGTTAGCACCCCCGGTCAGTGACAGGTCGCACTTCTGGCTAGAATCTATAGAATAAAATTGGATTATCTCACAGCATGACAGATAAAACTTCCATGGTCCTACTGAACGTCGAGTCCACTCGGGATTATGGCAAGGATTCGTCACTGTTAACTTTCAAGAGGAGGTAATTCCAGACAGTCTTCGATTTATAACGGCTCATAAATGCGGGGAGCGGCGATTTGTCTTCAAGAACGCTGAAGCCATAAAGTTAGCAGGTTGAGAACGCATGTGCGGGCACCATGGTACAAGAACCCCAAGGGGGTTGGTATTACTGTTTATAGTAA"
    print(HammingDistance(p, q))

    
    # 1.3 at 9/14
    text = "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT"
    symbol = "CC"
    print(FasterSymbolArray(text, symbol))
    
    # 1.4 
    text = "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT"
    print(SkewArray(text))
    
    text = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
    print(MinimumSkew(text))
    
    
    # Hidden  1.3
    f = open("dataset_7_10.txt", "r")
    text = f.read()
    print(MinimumSkew(text))
    
    
    
    
    
    

if __name__ == "__main__":
    main()