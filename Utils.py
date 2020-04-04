# Frechet distance implementation is borrowed from
# https://gist.github.com/MaxBareiss/ba2f9441d9455b56fbc9
import math
import numpy
from scipy.spatial.distance import cdist

def euc_dist(p1, p2):
    return math.sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]))

def frechet_dist(P,Q):
    n = len(P)
    m = len(Q)
    ca = numpy.full((n, m), -1)

    ca[0, 0] = euc_dist(P[0], Q[0])

    for i in range(1, n):
        ca[i, 0] = max(ca[i - 1, 0], euc_dist(P[i], Q[0]))
    for j in range(1, m):
        ca[0, j] = max(ca[0, j - 1], euc_dist(P[0], Q[j]))
    for i in range(1, n):
        for j in range(1, m):
            ca[i, j] = max(min(ca[i - 1, j],
                               ca[i, j - 1],
                               ca[i - 1, j - 1]),
                           euc_dist(P[i], Q[j]))
    return ca[n-1, m-1]

def euc_matrix(P, Q):
    mdist = cdist(P, Q, 'euclidean')
    return mdist