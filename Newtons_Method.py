'''
Lab 14 - Newton's Method.
'''

import numpy as np
from matplotlib import pyplot as plt

def Newtons_method(f, x0, Df, iters=15, tol=.002):
    '''Use Newton's method to approximate a zero of a function.
    
    INPUTS:
    f     - A function handle. Should represent a function from 
            R to R.
    x0    - Initial guess. Should be a float.
    Df    - A function handle. Should represent the derivative 
            of `f`.
    iters - Maximum number of iterations before the function 
            returns. Defaults to 15.
    tol   - The function returns when the difference between 
            successive approximations is less than `tol`.
    
    RETURN:
    A tuple (x, converged, numiters) with
    x           - the approximation to a zero of `f`
    converged   - a Boolean telling whether Newton's method 
                converged
    numiters    - the number of iterations the method computed
    '''
    x0 = float(x0)
    for i in range(iters):
        xk = x0 - f(x0)/Df(x0)
        #print xk
        if abs(x0 - xk) < tol:
            return xk, True, i + 1
        else:
            x0 = xk
    return xk, False, iters
    
def test():
    f = lambda x: x**2 - 1
    Df = lambda x: 2*x
    x0 = 1.5
    print Newtons_method(f, x0, Df)
    
def prob2():
    '''
    Print the answers to the questions in problem 2.
    '''
    f = lambda x: np.cos(x)
    Df = lambda x: -np.sin(x)
    Newtons_method(f, 1, Df, 15, 1e-5)
    Newtons_method(f, 2, Df, 15, 1e-5)
    print '1. For x0 = 1, 4 iterations are required for 5 digits of accuracy. For x0 = 2, 3 iterations are required for 5 digits of accuracy.'
    x = np.linspace(-4, 4, 500)
    y = np.sin(x)/x - x
    #plt.plot(x, y)
    plt.show()
    g = lambda x: np.sin(x)/x - x
    Dg = lambda x: -(x**2 + np.sin(x) - x*np.cos(x))/x**2
    n = Newtons_method(g, 1, Dg, 15, 1e-7)
    print '2.' + str(n[0])
    h = lambda x: x**9
    Dh = lambda x: 9*x**8
    Newtons_method(h, 1, Dh, 100, 1e-5)
    print '3. It takes 81 iterations to get five digits of accuracy. It is slow, because x^9 is flat from about -0.5 to 0.5.'
    k = lambda x: np.sign(x)*np.power(np.abs(x), 1./3)
    Dk = lambda x: (1/3.)*np.sign(x)*np.power(np.abs(x), -2./3)
    #print Newtons_method(k, .01, Dk, 50)
    print '4. There is oscillation happening here and Newton\'s method doesn\'t work because this case diverges.'
    #print '5.'
    
def Newtons_2(f, x0, iters=15, tol=.002):
    '''
    Optional problem.
    Re-implement Newtons method, but without a derivative.
    Instead, use the centered difference method to estimate the derivative.
    '''
    raise NotImplementedError('Newtons Method 2 not implemented')

def plot_basins(f, Df, roots, xmin, xmax, ymin, ymax, numpoints=100, iters=15, colormap='brg'):
    '''Plot the basins of attraction of f.
    
    INPUTS:
    f       - A function handle. Should represent a function 
            from C to C.
    Df      - A function handle. Should be the derivative of f.
    roots   - An array of the zeros of f.
    xmin, xmax, ymin, ymax - Scalars that define the domain 
            for the plot.
    numpoints - A scalar that determines the resolution of 
            the plot. Defaults to 100.
    iters   - Number of times to iterate Newton's method. 
            Defaults to 15.
    colormap - A colormap to use in the plot. Defaults to 'brg'. 
    
    RETURN:
    Returns nothing, but should display a plot of the basins of attraction.
    '''
    xreal = np.linspace(xmin, xmax, numpoints)
    ximag = np.linspace(ymin, ymax, numpoints)
    Xreal, Ximag = np.meshgrid(xreal, ximag)
    Xold = Xreal + 1j * Ximag
    diff = 10
    for i in xrange(iters):
        Xnew = Xold - f(Xold)/Df(Xold)
        diff = np.max(np.absolute(Xnew - Xold))
        Xold = Xnew
    
    
    def match(x):
        return np.argmin(abs(roots - [x]*len(roots)))
        
    Xnewnew = np.vectorize(match)(Xnew)
    
    plt.pcolormesh(Xreal, Ximag, Xnewnew, cmap=colormap)
    plt.show()
    
def test1():
    f = lambda x: x**3-x
    Df = lambda x: 3*x**2 - 1
    plot_basins(f, Df, np.array([-1, 0, 1]), -1.5, 1.5, -1.5, 1.5, 1000, 50)
    
def prob5():
    '''
    Using the function you wrote in the previous problem, plot the basins of
    attraction of the function x^3 - 1 on the interval [-1.5,1.5]X[-1.5,1.5]
    (in the complex plane).
    '''
    f = lambda x: x**3 - 1
    Df = lambda x: 3*x**2
    roots = np.array([1, -1j**(1./3), 1j**(2./3)])
    xmin, xmax, ymin, ymax = -1.5, 1.5, -1.5, 1.5
    plot_basins(f, Df, roots, xmin, xmax, ymin, ymax, 1000, 50)