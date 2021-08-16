# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)
def Motifs(Profile, Dna):
    motifs = []
    t = len(Dna)
    k = 3
    for i in range(t):
        motif = ProfileMostProbablePattern(Dna[i], k, Profile)
        motifs.append(motif)
    return motifs
# Insert your ProfileMostProbablePattern(Text, k, Profile) and Pr(Pattern, Profile) functions here.
def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p

def ProfileMostProbablePattern(Text, k, Profile):
    p_dict = {}
    for i in range(len(Text)- k +1):
        p = Pr(Text[i: i+k], Profile)
        p_dict[i] = p
    m = max(p_dict.values())
    keys = [k for k,v in p_dict.items() if v == m]
    ind = keys[0]
    return Text[ind: ind +k]


text = ["AAGCCAAA","AATCCTGG","GCTACTTG","ATGTTTTG"]
k = 4
Profile = {
    "A": [0, 0 , 0.25],
    "C": [0.75, 0.5, 0],
    "G": [0, 0, 0.25],
    "T": [0.25, 0.5, 0.5]
}

x = Motifs(Profile, text)
print(x)


text = ["TGACGTTC", "TAAGAGTT", "GGACGAAA", "CTGTTCGC"]
k = 3
Profile = {
    "A": [0, 0.25, 0.5],
    "C": [0, 0, 0],
    "G": [0.5, 0.5, 0],
    "T": [0.5, 0.25, 0.5]
}

x = Motifs(Profile, text)
print(*x)