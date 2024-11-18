import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson as poisson_dist
from tqdm import tqdm

class Q2:
    def __init__(self, lambda_param=25, C=32):
        self.x0 = 12345
        self.x1 = 12345
        self.x2 = 12345
        self.y0 = 12345
        self.y1= 12345
        self.y2 = 12345
        self.lambda_param = lambda_param
        self.C = C
        self.F25 = np.array([poisson_dist.cdf(i, self.lambda_param) for i in range(50)],dtype=float)
        self.F25 = np.append(self.F25, 1)
        self.inf_i = np.searchsorted(self.F25, np.linspace(0, 1, C, endpoint=False))

    def MRG32K3a(self):
        xn = int((1403580 * self.x1 - 810728 * self.x0) % 4294967087)
        yn = int((527612 * self.y2 - 1370589 * self.y0) % 4294944443)
        un = float(((xn - yn) % 4294967087) / 4294967087)

        self.x0 = self.x1
        self.x1 = self.x2
        self.x2 = xn
        self.y0 = self.y1
        self.y1 = self.y2
        self.y2 = yn

        return un

    def recherche_seq(self, U):
        i = -1
        while i<50:
            i += 1
            if poisson_dist.cdf(i, self.lambda_param) >= U:
                return i
        return i

    def recherche_index(self, U):
        s = int(self.C * U)
        i = self.inf_i[s]
        while i < 50 and self.F25[i] < U:
            i += 1
        return i

    def generator_seq(self, n):
        cdf = {i: 0 for i in range(51)}
        for _ in tqdm(range(n), leave=True):
            u = self.MRG32K3a()
            cdf[self.recherche_seq(u)] += 1
        return cdf

    def generator_index(self, n):
        cdf = {i: 0 for i in range(51)}
        for _ in range(n):
            u = self.MRG32K3a()
            cdf[self.recherche_index(u)] += 1
        return cdf


lambda_param = 25
q2 = Q2(lambda_param = lambda_param)


#cdf = q2.generator_index(100_000)
cdf = q2.generator_seq(100_000)


i = list(cdf.keys())
pmf = list(cdf.values())

plt.bar(i, pmf,  color='skyblue', edgecolor='black')
plt.xlabel('i')
plt.ylabel('P(X = i)')
plt.title(f' pmf (k = {lambda_param})')
plt.show()
