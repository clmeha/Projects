# spec.py
"""Volume 1, Lab 16: Importance Sampling and Monte Carlo Simulations.
<Carol>
<Math347>
<1/26/16>
"""

import numpy as np
import scipy.stats as stats
from matplotlib import pyplot as plt

# Problem 1 
def prob1(n):
    """Approximate the probability that a random draw from the standard
    normal distribution will be greater than 3.
    Returns: your estimated probability.
    """
    h = []
    k = np.random.normal(loc = 0, scale = 1, size = n)
    for i in range(len(k)):
        if k[i] > 3:
            h.append(1)
        else:
            h.append(0)
    return 1./n*np.sum(h)
        
    

# Problem 2
def prob2():
    """Answer the following question using importance sampling: 
            A tech support hotline receives an average of 2 calls per 
            minute. What is the probability that they will have to wait 
            at least 10 minutes to receive 9 calls?
    Returns:
        IS (array) - an array of estimates using 
            [5000, 10000, 15000, ..., 500000] as number of 
            sample points."""
    N = np.arange(5000, 500001, 5000)
    def importance(N):
        h = lambda x: x > 10
        f = lambda x: stats.gamma(9, scale = 0.5).pdf(x)
        g = lambda x: stats.norm(loc = 10).pdf(x)
        X = np.random.normal(10, size = N)
        return 1./N * np.sum(h(X)*f(X)/g(X))
    return np.vectorize(importance)(N)
    

# Problem 3
def prob3():
    """Plot the errors of Monte Carlo Simulation vs Importance Sampling
    for the prob2()."""
    h = lambda x: x > 10
    MC_estimates = []
    for N in xrange(5000, 505000, 5000):
        X = np.random.gamma(9, scale = 0.5, size = N)
        MC = 1./N*np.sum(h(X))
        MC_estimates.append(MC)
    MC_estimates = np.array(MC_estimates)
    m = 1 - stats.gamma(a = 9, scale = 0.5).cdf(10)
    x = np.arange(5000, 500001, 5000)
    y = abs(MC_estimates - m)
    Y = abs(prob2() - m)
    plt.plot(x, y)
    plt.plot(x, Y)
    plt.show()
    
# Problem 4
def prob4():
    """Approximate the probability that a random draw from the
    multivariate standard normal distribution will be less than -1 in 
    the x-direction and greater than 1 in the y-direction.
    Returns: your estimated probability"""
    N = 500000
    h = lambda x: x[0] < -1 and x[1] > 1
    f = lambda x: stats.multivariate_normal(mean = [0, 0]).pdf(x)
    g = lambda x: stats.multivariate_normal(mean = [-1, 1]).pdf(x)
    X = np.random.multivariate_normal(mean = [-1, 1], cov = [[1, 0], [0, 1]], size = N)
    tot = 0
    for i in range(N):
        tot += (h(X[i])*f(X[i])/g(X[i]))
    return 1./N * tot
    