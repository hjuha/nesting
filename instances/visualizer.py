import matplotlib.pyplot as plt
import numpy as np

s = input().split()
if len(s) == 3:
	n, _, _ = s
	n = int(n)
	a = np.zeros((n, n))
	for i in range(n):
		s = input().split()
		for j in range(n):
			a[i, j] = float(s[j])
	plt.matshow(a)
	plt.show()
else:
        n = int(input())
        a = np.zeros((n, n))
        for i in range(n):
                s = input().split()
                for j in range(n):
                        a[i, j] = float(s[j])
        plt.matshow(a)
        plt.show()
