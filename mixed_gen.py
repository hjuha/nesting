import random
from math import log

n = 50
s = 40
d = 20
t = 10000
print(n, d, t)

import numpy as np

A = np.zeros((n, n))
for i in range(s // 3):
	y1 = random.randint(0, s - 1)
	x1 = random.randint(0, s - 1)
	x2 = random.randint(0, s - 1)
	if x1 > x2:
		x1, x2 = x2, x1
	for i in range(x2 - x1 + 1):
		if y1 + i < s:
			A[x1 + i, y1 + i] = A[y1 + i, x1 + i] = random.random() / 4

for i in range(s):
	A[i, i] = 0.5 + random.random() / 2

for i in range(s, n):
        for j in range(1 + int(log(n - s)/log(5))):
                A[i, random.randint(s, n - 1)] = .5 + random.random() / 2

for i in range(n):
	for j in range(n):
		print(A[i, j], end=" ")
	print("")
