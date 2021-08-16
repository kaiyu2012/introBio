import random
import math
############## DATA ##############

#k = 8
#t = 5
#N = 100

#dna = ['CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA',
#'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
#'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
#'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
#'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']

############## STDIN ##################

k = 15
t = 20
N = 2000

filename = 'dataset_163_4.txt'
with open(filename, "r") as dataset:
    dna = []
    for line in dataset:
        dna.append(line.strip())
    Text = dna[0]

############## FUNCIONES PARALELAS ################

def _randomkmers(dna,k,t):
    randomkmer = []
    for i in range(t):
        z = random.choice(range(0,len(dna[i])-k+1))
        randomkmer.append(dna[i][z:z+k])
    return randomkmer

def _consensus(motifs):
    k = len(motifs[0])
    count = _count(motifs)
    consensus = ""
    for j in range(k):
        M = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > M:
                M = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def hamming_distance(p, q):
    count = 0
    L = len(p)
    for i in range(L):
        if p[i] != q[i]:
            count += 1
    return count

def _profile(motifs):
    profile = {}
    t = len(motifs)
    k = len(motifs[0])
    countMotifs = _count(motifs)
    for symbol in "ACGT":
        profile[symbol] = []
    for x in countMotifs:
        for y in countMotifs[x]:
            z = y / float(t)*2
            profile[x].append(z)
    return profile

def _count(motifs):
    count = {}
    k = len(motifs[0])
    for symbol in "ACGT":
        count[symbol] = [] #Genero una lista para cada nucleotido en el set count
        for j in range(k):
            count[symbol].append(1) #a cada uno le pongo un 0
    t = len(motifs)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j] #para el simbolo de esa posicion del motivo
            count[symbol][j] += 1 #sumarle un 1 al set count en ese lugar
    return count

def _motifs(profile,dna,k):

    motifs = []
    def maxProb(pattern,profile):
        num = []
        for i in range(len(pattern)):
            proflist = profile[pattern[i]][i]
            num.append(proflist)
        product = math.prod(num)
        return product
    for x in range(len(dna)):      #for i in range t
        prob = 0
        maxpattern = ''
        for i in range(len(dna[x])-k+1):
            pattern = dna[x][i:i+k]  #kmer text
            if maxProb(pattern,profile) > prob:
                prob = maxProb(pattern,profile)
                maxpattern = pattern
        motifs.append(maxpattern)
    return motifs

def _score(motifs):
    score = 0
    consensus = _consensus(motifs)
    for i in range(len(motifs)):
        score += hamming_distance(consensus,motifs[i])
    return score

def profile_most_probable_kmer(text, k, profile):
    mostProbVal = -1
    mostProbKmer = ''

    for i in range(0, 1 + len(text) - k):
        kmer = text[i:i+k]
        probKmerVal = _pr(kmer, profile)
        if probKmerVal > mostProbVal:
            mostProbVal = probKmerVal
            mostProbKmer = kmer

    return mostProbKmer

def _pr(text, profile):
    P = 1

    for i in range(len(text)):
        P = P * profile[text[i]][i]

    return P

def mostprobprof(Text,k,profile):
    prob = 0
    maxpattern = ''
    def maxProb(pattern,profile):
        num = []
        for i in range(len(pattern)):
            proflist = profile[pattern[i]][i]
            num.append(proflist)
        product = math.prod(num)
        return product
    for i in range(len(Text)-k+1):
        pattern = Text[i:i+k]
        if prob < maxProb(pattern,profile):
            prob = maxProb(pattern,profile)
            maxpattern = pattern
    return maxpattern
############## FUNCION PRINCIPAL ##################

def GibbsSampler(dna,k,t,N):
    motifs = _randomkmers(dna, k, t)
    bestMotifs = motifs
    for i in range(N):
        i = random.choice(range(t))
        del motifs[i]
        profile = _profile(motifs)
        motifs.insert(i,mostprobprof(dna[i],k,profile))
        if _score(motifs) < _score(bestMotifs):
            bestMotifs = motifs
    return bestMotifs

#print(GibbsSampler(dna,k,t,N))

bestScore = 100
for start in range(20):
    bMotifs = GibbsSampler(dna,k,t,N)
    if _score(bMotifs) < bestScore:
        bestScore = _score(bMotifs)
        bestMotifs = bMotifs

print('\n'.join(bestMotifs))