import random
from math import log

n = int(input())
d = 15
t = 120
print(n, d, t)

"""
for i in range(n):
	a = [random.random() * (2 - min(2, max(1, abs(i - j)))) for j in range(n)]
	for j in range(n):
		print(a[j] / max(a), end=" ")
	print("")
"""

"""
for i in range(n):
	a = [random.random() * (i // 6 == j // 6) for j in range(n)]
	for j in range(n):
		print(a[j] / max(a), end=" ")
	print("")
"""

import numpy as np

A = np.zeros((n, n))
for i in range(n):
	for j in range(1 + int(log(n)/log(2))):
		A[i, random.randint(0, n - 1)] = .5 + random.random() / 2
for i in range(n):
	for j in range(n):
		print(A[i, j], end=" ")
	print("")

"""A = np.zeros((n, n))

for i in range(n // 3):
	y1 = random.randint(0, n - 1)
	x1 = random.randint(0, n - 1)
	x2 = random.randint(0, n - 1)
	if x1 > x2:
		x1, x2 = x2, x1
	for i in range(x2 - x1 + 1):
		if y1 + i < n:
			A[x1 + i, y1 + i] = A[y1 + i, x1 + i] = random.random() / 4

for i in range(n):
	A[i, i] = 0.5 + random.random() / 2

for i in range(n):
	for j in range(n):
		print(A[i, j], end=" ")
	print("")"""
