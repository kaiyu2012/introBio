def Normalize(P):
    d = {}
    for k,v in P.items():
        d[k] = P[k]/sum(P.values())
    return d


def normalize(P):
    d = {}
    for i in range(len(P)):
        d[i] = P[i]/sum(P)
    return d


x = [0.15, 0.6, 0.225, 0.225, 0.3]

nz = normalize(x)
print(nz)