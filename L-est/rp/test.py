import bayes
from scipy.stats.distributions import norm as n
import numpy as np

e = Monitor('x', 'y')
e.update({
    'x' : np.random.rand(),
    'y' : np.random.rand(),
    'z' : np.random.rand()
    })

o = MCMC(e)

@o.gibbs('x')
def xUp():
    np.random.rand()

@o.gibbs('y')
def yUp():
    np.random.rand()

@o.metropolis('z')
def zLogLik(x, y, z):
    return log(n.pdf(x*y, loc = z, scale = 0.2))
