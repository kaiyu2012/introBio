p1 = float ((600.-15) / (600.-15+1))
p2 = 1 - p1
from itertools import combinations
counter = 0
for seq in combinations(range(10),2):
    counter +=1
# counter
import scipy.special as sc
sc.comb(10, 2, exact=True)

output = pow(p2,2) * pow(p1,8) * counter
print(output)