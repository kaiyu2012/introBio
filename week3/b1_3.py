import numpy as np

# input: a list of strings Motifs
# output: returns the count matrix of  Motifs (as a dictionary of lists)
def count_motifs(motifs) -> dict:
    count = {}
    k = len(motifs[0])
    for symbol in "ACGT":
        count[symbol] = np.zeros(k, int).tolist()

    t = len(motifs)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1
    return count



# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(motifs):
    count = {}
    profile = {}
    k = len(motifs[0])
    for symbol in "ACGT":
        count[symbol] = np.zeros(k, int).tolist()
        profile[symbol] = np.zeros(k, float).tolist()


    t = len(motifs)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1

    
    for symbol in "ACGT":
        for i in range(k):
            profile[symbol][i] = round(count[symbol][i]/t, 1)
    
    return profile



# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    ct = count_motifs(Motifs)
    consensus = ""
    symbol = "ACGT"

    
    for i in range(len(Motifs[0])):        
        x = {symbol[0]:ct[symbol[0]][i], 
                symbol[1]:ct[symbol[1]][i], 
                symbol[2]:ct[symbol[2]][i], 
                symbol[3]:ct[symbol[3]][i]}
        consensus += max(x, key=x.get)
    
    return consensus


# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    ct = count_motifs(Motifs)
    symbol = "ACGT"
    score = 0    
    for i in range(len(Motifs[0])):        
        x = {symbol[0]:ct[symbol[0]][i], 
                symbol[1]:ct[symbol[1]][i], 
                symbol[2]:ct[symbol[2]][i], 
                symbol[3]:ct[symbol[3]][i]}        
        score += max(x.values())
    
    return(len(Motifs[0])*len(Motifs) - score)



def main():
    m = ["AACGTA","CCCGTT","CACCTT","GGATTA","TTCCGG"]
    print(count_motifs(m))
    print(Profile(m))
    print(Consensus(m))
    Score(m)
    
    
    

if __name__ == "__main__":
    main()