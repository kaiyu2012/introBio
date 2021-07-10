
from typing import Pattern
import numpy as np

# count the frequent words
def count_fre(text:str, pattern:str) -> int:
    x = text
    y = pattern
    lx = len(x)
    ly = len(y)

    count = 0

    for i in range(lx - ly + 1):
        if x[i:i+ly] == y:
            count += 1
        
    print(count)
    return count


# find the most frequent words in text
def freq_words(text:str, k:int) -> dict:
    freq_table = []
    freq_count = []
    for i in range(len(text) - k + 1):
        x = text[i:i+k]
        new_patt = True
        for j in range(len(freq_table)):
            if x == freq_table[j]:
                freq_count[j] += 1
                new_patt = False
                break
        if new_patt:
            freq_table.append(x)
            freq_count.append(1)
    
    freq = np.array(freq_count)
     
    
    w_index = np.argwhere(freq == np.amax(freq))
 
    most_freq = []

    for i in w_index:
        most_freq.append(freq_table[i[0]])
        
    d = {'freq': freq_count[w_index[0][0]],
         'list': most_freq}
    
    return d


# find the specific frequence words in text
def spec_freq_words(text:str, k:int, t:int) -> list:
    freq_table = []
    freq_count = []
    for i in range(len(text) - k + 1):
        x = text[i:i+k]
        new_patt = True
        for j in range(len(freq_table)):
            if x == freq_table[j]:
                freq_count[j] += 1
                new_patt = False
                break
        if new_patt:
            freq_table.append(x)
            freq_count.append(1)
    
    freq = np.array(freq_count)
    d = []
    for i in range(len(freq)):
        if freq[i] >= t:
            d.append(freq_table[i])    
    
    return d
            
# get the reverse complement
def rev_comp(text:str) -> str:
    r = ''
    l = len(text)
    for i in range(l):
        char = text[l - i - 1]
        if char == 'A':
            r += 'T'
        elif char == 'T':
            r += 'A'
        elif char == 'G':
            r += 'C'
        elif char == 'C':
            r += 'G'
    
    return r


# match the pattern in text
def match_patt(patt:str, text:str) -> list:
    lp = len(patt)
    lt = len(text)
    ind = []
    for i in range(lt - lp + 1):
        if text[i:i+lp] == patt:
            ind.append(i)
    return ind   


def match_patt_size(patt:str, text:str, t: int) -> bool:
    lp = len(patt)
    lt = len(text)
    ind = 0
    for i in range(lt - lp + 1):
        if text[i:i+lp] == patt:
            ind += 1
            if ind >= t:
                break
    if ind >= t:
        return True
    else:
        return False
  
  
  
def match_patt_file(patt:str, file:str) -> list:
    f = open(file, "r")
    text = f.read()
    print(text)
    
    return(match_patt(patt, text))


# find a clump
def clump_finding(text:str, k:int, L:int, t:int) -> list:
    
    lt = len(text)
    x = spec_freq_words(text[0:L], k, t)
    
    for i in range(L, lt + 1):
        p = text[i-k: i]
        r = text[i-L: i]
        if p not in x and match_patt_size(p, r, t):
            x.append(p)
    
    return x
    
    
    
def main():
    # f = open("E_coli.txt", "r")
    # text = f.read()
    # k = 9
    # L = 500
    # t = 3
    
    # f = clump_finding(text, k, L, t)
    
    # print(len(f))
    
    
    # text = "TGTCTAACGGAGTTTGAGCGAGCAAGTGCCAATGAAGTGCCAATTAGCCAAGCGCCTGGGGACCACATGGCATTCATTCATTACATGCACATGCCGTTATCCGCTACTCCGATAAGCCGTCCATCATGCGCACGAACGCGTCGCTCAGTGTATGACGTACACAACTTTATTATGACTACGGAACTTGGGAGAGAGAAGAAATAAGTTACTAATGAGCCGGCCTATGAAGTCTACGTATCTTCGCAGTTATCAATTAGCGGACAGAAGTTCAAATACCGCTATGCTTATTAATATGTTTTTTGTCCTAACTCCTTGTCCTAACTTATTGGTTCCCACAAAGCATAGAAACCCCATTATAGTTAGGTTTGGAGCGTATATTACCAGAAATAGTACCTTGAGGCAATCATGGCCTAATACTAAGGACGTTTAAATAATAGAGTCGTCCGCTGGGGCTTATAACACAGTCGGTTGCAGACTACGGACGTGGAACAGCGATCCGCCGAGCTGGCTAGGGGCAGGCCTAGGGGACATAACGATTCACTTCACTCAACACTTCACTCAACAAGTCATCGAGCATCGAGAGAATCAGCAGTAAGTTAGACCACCCGCCGCATCTTACATCTTAGCATCTTAGGGGTACGAGGGGTATATCACGCGTGTGGGATCATACCGCACTAGACAACGACTACCTAGTAGGATCGCTAAGTCCCGGTTCAGTTAAGTTCAGTTAAAAGCCAAGCACACTACTCCTCTGGGCGTAGTTGTTTCATCAATCCTATTATAGCCCCGTAGTGCTGTGAAGTTCGCTGCTGACAGCCCTGACAGCCCATAAGCGCGGCGACAAGATAACTCCAACTCCCTCCTCGGGCGGGGTACGGGGTATCGGCAACAGCTCCAGAGTCCCTTTAGGGGCCCGCATAGGCACCAGGCATTGTTGATGTAACGAGATCAACCCATTCCATCGCACTTCCTCTCTGATAAACTAAAGAGTATTGGGTGAAATCAAGACATTTGCATAGTTCCCACCAGACAGACCGAGCGACCGAGCCACCAGCTCAAACAGTTTCTAAAGTTCTTGCGATGATAAGTGAGATTTGCGGATTTGAGATTAGCGCGTATGTAAATCTAGTGATCTACTCAGCGAGGCCTCAGACGCGTTCTCATCCCAAAGTTTATATTTGATGGGACCCAGCGTTGCGACGGATATAGGGATCCAGCGCAGAACATAGAATACCCTTTATTCTCTATGAAACATTAGGCAGTACAAAATTTCCAAAATTTCCAAAATTTCCGTGCGCCGCTCTCCGGTATGATGTCTGCGCTAGAGTCTAGCAATGAACGCTCCCCACTCTGAACCACTCTGAAGAACCCAACGGGTTCCGGCAGAACTCTCTAGGCGGCGTCTGATTAAATCTGACGTCCCAAGGCGCCGGTGAAGGAGGCATCATCACACCAGTGCTGGAGGTCGTAGACGAATATATAATCTTTCTTGCATCGCAGCCGCAGCAATGTTGTGTAGCACACCACGTTCCATTCACCGCAGGGGCTACAAGACGTGAATTGTTAGATGCTACTGGCAGTCCTATACCTAGACGTTCGAGCAGGCATCTGTTCACGCGTTAGGTGCTCCCCCAAGTAAAACCTCGGCTCATACATATTTAGTCGTCTCTCGTCTCGTCCGCTTTCCGCTTAGGGTTTACATATTAAACTGCCCCTACGAGCCCCTACGAGTGAGCAGCGCGCTGAAAGTAGACGTAGAGTAAATCGGCGCCTTCGCCCACACCCACAGTGTGTAGCCGACCGCTGTATCATTATATTCGTTCTTCATGTTCCGGTCATGTTCCGGGGTGAAAGTGCATTTAAAAAGTGCATT"
    # k = 10
    # L = 29
    # t = 3
    # f = clump_finding(text, k, L, t)
    # print(f)
    
    text = "CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT"
    q = freq_words(text, k=3)
    print(q)



if __name__ == "__main__":
    main()
    