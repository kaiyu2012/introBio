import numpy as np
import math

CHARS = "ACGT"

# Input: motifs
# Output: profiles
def computer_profile(motifs) -> dict:    
    counts = {}
    profiles = {}
    columns = len(motifs[0])
    rows = len(motifs)

    for c in CHARS:
        counts[c] = np.zeros(columns)
        profiles[c] =  list(np.zeros(columns))

    for i in range(columns):
        for j in range(rows):
            counts[motifs[j][i]][i] += 1
    
    for i in range(columns):
        s = counts[CHARS[0]][i] + counts[CHARS[1]][i] + counts[CHARS[2]][i] + counts[CHARS[3]][i]
        for x in CHARS:
            profiles[x][i] = counts[x][i] / s
        
    return profiles
    


# Input: motifs
# Output: entropy
def compute_entropy(motifs) -> float:
    p = computer_profile(motifs)
    e = {}
    columns = len(motifs[0])
    for char in CHARS:
        for i in range(columns):
            x = p[char][i]
            if x:
                if x in e:
                    e[x] += 1
                else:
                    e[x] = 1

    ent = 0
    for x in e:        
        ent +=  e[x] * x * math.log2(x)
    
    return ent




def main():

    motifs = [
        "TCGGGGGTTTTT",
        "CCGGTGACTTAC",
        "ACGGGGATTTTC",
        "TTGGGGACTTTT",
        "AAGGGGACTTCC",
        "TTGGGGACTTCC",
        "TCGGGGATTCAT",
        "TCGGGGATTCCT",
        "TAGGGGAACTAC",
        "TCGGGTATAACC"
        ]

    ent = compute_entropy(motifs)
    print(ent)

    
    
    

if __name__ == "__main__":
    main()