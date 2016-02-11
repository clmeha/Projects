# spec.py
"""Volume I: Monte Carlo Integration
<Carol>
<Math347>
<1/19/16>
"""

import numpy as np

# Problem 1
def prob1(numPoints = 10.**5):
    """Return an estimate of the volume of the unit sphere using Monte
    Carlo Integration.
    
    Inputs:
        numPoints (int, optional) - Number of points to sample. Defaults
            to 10^5.
    Returns:
        volume (int) - Approximate value of the area of the unit sphere.
    """
    points = np.random.rand(3, numPoints)
    points = points*2. - 1
    circleMask = np.sqrt(points[0,:]**2 + points[1,:]**2 + points[2,:]**2) <= 1
    numInCircle = np.count_nonzero(circleMask)
    return 8.*numInCircle/numPoints

# Problem 2
def prob2(numPoints = 10.**5):
    """Return an estimate of the area under the curve,
    f(x) = |sin(10x)cos(10x) + sqrt(x)*sin(3x)| from 1 to 5.
    
    Inputs:
        numPoints (int, optional) - Number of points to sample. Defautls
            to 10^5.
    Returns:
        area (int) - Apprimate value of the area under the 
            specified curve.
    """
    x_val = np.random.rand(int(numPoints))
    f = lambda x: abs(np.sin(10.*x)*np.cos(10.*x) + np.sqrt(x)*np.sin(3*x))
    x_val = x_val*4 + 1
    return (5-1)*sum(f(x_val))/(numPoints)
    
    

# Problem 3
def mc_int(f, mins, maxs, numPoints=500, numIters=100):
    """Use Monte-Carlo integration to approximate the integral of f
    on the box defined by mins and maxs.
    
    Inputs:
        f (function) - The function to integrate. This function should 
            accept a 1-D NumPy array as input.
        mins (1-D np.ndarray) - Minimum bounds on integration.
        maxs (1-D np.ndarray) - Maximum bounds on integration.
        numPoints (int, optional) - The number of points to sample in 
            the Monte-Carlo method. Defaults to 500.
        numIters (int, optional) - An integer specifying the number of 
            times to run the Monte Carlo algorithm. Defaults to 100.
        
    Returns:
        estimate (int) - The average of 'numIters' runs of the 
            Monte-Carlo algorithm.
                
    Example:
        >>> f = lambda x: np.hypot(x[0], x[1]) <= 1
        >>> # Integral over the square [-1,1] x [-1,1]. Should be pi.
        >>> mc_int(f, np.array([-1,-1]), np.array([1,1]))
        3.1290400000000007
    """
    l = []
    dim = len(mins)
    #points = points * (maxs - mins) + mins
    scale = np.diag(maxs - mins)
    product = np.prod(maxs - mins)
    for i in xrange(numIters):
        points  = np.random.rand(dim, numPoints)
        points = scale.dot(points) + np.diag(mins).dot(np.ones_like(points))
        estimate = product*sum(np.apply_along_axis(f, 0, points))/float(numPoints)
        l.append(estimate)
    return np.mean(l)
    
def test():
    f = lambda x: np.hypot(x[0], x[1]) <= 1
    print mc_int(f, np.array([-1, -1]), np.array([1, 1]), numPoints = 10, numIters = 10)
    
# Problem 4
def prob4(numPoints=[500]):
    """Calculates an estimate of the integral of 
    f(x,y,z,w) = sin(x)y^5 - y^5 + zw + yz^3
    
    Inputs:
        numPoints (list, optional) - a list of the number of points to 
            use in the approximation. Defaults to [500].
    Returns:
        errors (list) - a list of the errors when calculating the 
            approximation using 'numPoints' points.
    Example:
    >>> prob4([100,200,300])
    [-0.061492011289160729, 0.016174426377108819, -0.0014292910207835802]
    """
    sample_points = [100, 1000, 10000]
    f = lambda x: np.sin(x[0])*x[1]**5 - x[1]**3 + x[2]*x[3] + x[1]*x[2]**3
    mins = np.array([-1, -1, -1, -1])
    maxs = np.array([1, 1, 1, 1])
    l = []
    for i in sample_points:
        k = mc_int(f, mins, maxs, numPoints = i)
        l.append(k)
    return [abs(i) for i in l]
        