import numpy as np

class Sphere():
    """
    Creates a sphere class for uniform sampling at a constant radius.
    """
    def __init__(self, r = 1, shift = np.zeros(3)):
        if 'shape' in dir(r) or not shift.shape == (3,):
            raise ValueError("r must be a number and shift must be a (3,) array")
        self.r = r
        self.shift = np.array(shift)

    def sample(self, n = 1):
        """
        Uses the 'trigonometry trick' to sample randomly from the sphere.

        >>> r = np.random.uniform(low = 1, high = 20)
        >>> s1 = Sphere(r)
        >>> dists = np.sqrt((s1.sample(100)**2).sum(axis = 1))
        >>> all(abs(dists - r) < 0.0001)
        True

        >>> shift = np.random.uniform(size = (3,), low = -10, high = 10)
        >>> s2 = Sphere(r, shift = shift)
        >>> centers = np.mean(s2.sample(5000), axis = 0)
        >>> np.linalg.norm(centers - shift) < r/20
        True
        """
        zs = np.random.uniform(low = -1, high = 1, size = n)
        ts = np.random.uniform(low = 0, high = 2*np.pi, size = n)
        out = np.ndarray(shape = (n, 3))
        for i in xrange(n):
            z = zs[i]
            t = ts[i]
            r = np.sqrt(1 - z**2)
            out[i, :] = [r * np.cos(t),
                         r * np.sin(t),
                         z]
        return (out * self.r) + self.shift

    def sampleGen(self):
        """
        Generator method for creating samples.
        """
        while True:
            z = np.random.uniform(low = -1, high = 1)
            t = np.random.uniform(low = 0, high = 2*np.pi)
            r = np.sqrt(1 - z**2)
            x, y = r * cos(t), r * sin(t)

            yield (np.array([x, y, z]) * self.r) + self.shift
    

class Ellipsoid():
    """
    Shapes holding the constraint:
    
    x.t A x = 1
    for some positive definite matrix A.
    
    """
    def __init__(self, shapeMat = np.diag([1,1,1]), shift = np.zeros(3)):
        """
        shapeMat is a 3x3 positive definite array analogue to a
        MVNormal covariance matrix.
        """
        if (not shapeMat.shape == (3,3)) or (not shift.shape == (3,)):
            raise ValueError("shapeMat must be 3x3 and shift must be 3x1")
        self.A = shapeMat
        self.Ai = np.linalg.inv(self.A)
        self.Ac = np.linalg.cholesky(self.A)
        self.Aic = np.linalg.cholesky(self.Ai)

        # TODO: This doesn't actually find the true maxScaleF
        w, v = np.linalg.eig(self.Aic)
        self.maxScaleF = np.max(w ** 2)
        
        self.shift = shift

    def sample(self, n = 1):
        g = self.sampleGen()
        out = np.ndarray(shape = (n, 3))
        for i in xrange(n):
            out[i, :] = g.next()
        return out

    def sampleGen(self):
        g = Sphere().sampleGen()
        while True:
            prop = np.dot(g.next(), self.Aic)
            sf = self.scaleF(prop)
            if np.random.rand() < sf:
                yield (prop + self.shift)
            else:
                continue

    def scaleF(self, x):
        """
        The differential ratio of a projected surface area to its area
        on a unit sphere is:

        R = r1^2/cos(phi)

        where r1 is the projected radius and phi is the angle between
        the normal vector on the unit sphere and its projected point
        on the ellipsoid.
        """
        x = np.matrix(x)
        r1_sq = np.linalg.norm(np.dot(x, self.Aic)) ** 2
        elipnorm = 2 * np.dot(x, self.Ai)
        cosphi = np.linalg.norm(elipnorm)/np.dot(x, elipnorm.transpose())
        
        return ((r1_sq/cosphi)/self.maxScaleF)[0, 0]

def randomA(low = 1, high = 10):
    vals = np.random.uniform(low = low, high = high, size = (3,))
    q, r = np.linalg.qr(np.matrix(np.random.uniform(size = (3,3))))
    return q*np.matrix(diag(vals))*q.transpose()

n = 1000
g = Ellipsoid(diag([.004, .004, .001])).sampleGen()
out = np.ndarray(shape = (n, 3))
for i in xrange(n):
    out[i, :] = g.next()
np.savetxt('/Users/tel/Temp/test.csv', out, delimiter = ',')
        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
