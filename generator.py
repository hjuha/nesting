import random
from math import log2

group, n, J = input().split()[:3]
n = int(n)
J = min(int(n), int(J))
random.seed(n + 1337)

matrix = [[0 for _ in range(n)] for _ in range(n)]

print(0.1, 0.05, J, 120)
print(n)

if group == "block_diagonal":
        k = 5
        for y in range(n):
                for x in range(n):
                        if (y // k == x // k):
                                matrix[y][x] = random.random()
                print(" ".join([str(z) for z in matrix[y]]))
elif group == "permutation":
        matrix = [[0 for _ in range(n)] for _ in range(n)]
        for i in range(1 + int(log2(n))):
                p = [i for i in range(n)]
                random.shuffle(p)
                for j in range(n):
                        matrix[j][p[j]] = random.random()
        for y in range(n):
                print(" ".join([str(z) for z in matrix[y]]))