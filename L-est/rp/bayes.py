import numpy as np
import scipy.optimize
import collections
import inspect

class Monitor(dict):
    """
    Monitor presents a Dictionary interface that abstracts over a
    history storing functionality. On creation a (possibly empty) list
    of monitored keys is passed and each time the value of one of
    those keys is updated, it is added to a history queue. This queue
    can be accessed for any of the monitored keys via the history()
    method. If history() is called on an unmonitored key an error is
    raised.

    >>> m = Monitor('x', 'y', 'z')
    >>> m['x'] = 3
    >>> m['x'] = 4
    >>> m['x'] = 5
    >>> m['x'] = 6
    >>> m['y'] = 0
    >>> m['x']
    6
    >>> m['y']
    0
    >>> m.history('x')
    deque([3, 4, 5, 6])
    >>> m.history('m')
    Traceback (most recent call last):
    ...
    KeyError: 'm'
    >>> 'z' in m.keys()
    True
    >>> m['q'] = 3
    >>> m['q']
    3
    """

    def __init__(self, *watched):
        dict.__init__(self)
        self.watched = watched
        self._histories = dict()
        self.clearHistories()
        for k in self.watched:
            # Let the dictionary know each monitored key will be used
            # later.
            dict.__setitem__(self, k, None)

    def __setitem__(self, k, v):
        if k in self.watched:
            self._histories[k].append(v)
        return dict.__setitem__(self, k, v)

    def history(self, k):
        try:
            return self._histories[k]
        except KeyError:
            if k in self.watched:
                return None
            else:
                raise KeyError, k

    def clearHistories(self, maxlen = None):
        h = dict()
        for k in self.watched:
            h[k] = collections.deque(maxlen = maxlen)
        self._histories = h

    def view(self, keys):
        out = dict()
        for k in keys:
            out[k] = self[k]
        return out

class withEnvironment:
    """
    >>> e = Monitor('x', 'y')
    >>> e['x'] = 1
    >>> e['y'] = 2
    >>> @withEnvironment(e)
    ... def foo(x, y):
    ...   return (x + y)
    >>> foo()
    3
    """

    def __init__(self, env, retain = None):
        self.env = env
        self.retain = retain

    def __call__(self, f):
        args, varargs, varkw, defaults = inspect.getargspec(f)
        if self.retain and not self.retain in args:
            raise KeyError, "Function does not depend on variable %s" % self.retain
        
        needed = set(args).difference(self.env.keys())
        if needed:
            raise ValueError, "Environment does not contain value for keys: %s" % list(needed)

        if self.retain:
            def wrapped(x):
                args = self.env.view(args)
                args[self.retain] = x
                return f(**args)
        else:        
            def wrapped():
                return f(**self.env.view(args))
        return wrapped

class Sampler:
    def __init__(self, env, key):
        self.frozen = False        
        self.env = env
        self.key = key
        
    def __iter__(self):
        return self
    
    def __call__(self, f):
        self.f = withEnvironment(self.env)(f)

    def updateEnv(self, v):
        self.env[self.key] = v
    
    def next(self):
        try:
            v = f()
            return v
        except AttributeError:
            raise Error, "Sampler not yet initialized with function."

    def freeze(self):
        self.frozen = True

class Gibbs(Sampler):
    pass

class Metropolis(Sampler):
    """
    Univariate normal metropolis sampler. Is passed a
    logLiklihood instead of a sampling function.
    """
    def __init__(self, env, key, sd = None, memory = 20):
        Sampler.__init__(self, env, key)
        self.memory = memory
        self.sd = sd or self._findSD()
        self.jumps = collections.deque(maxlen = self.memory)

    def __call__(self, f):
        self.f = withEnvironment(self.env, retain = self.key)(f)

    def x0(self):
        return self.env[self.key]

    def _findSD(self):
        """
        Find starting SD by computing variance at loglik mode and
        multiplying by 2.4.
        """
        eps1 = 0.0001
        eps2 = 0.00005
        mode = scipy.optimize.fmin(lambda x: -self.f(x), self.x0(), disp = False)
        curv = (self.f(mode + eps1 + eps2) -
                self.f(mode - eps1 + eps2) -
                self.f(mode + eps1 - eps2) +
                self.f(mode - eps1 - eps2))/(4*eps1*eps2)
        return -2.4/curv

    def accept(self):
        self.jumps.append(1)
        if np.random.rand() < 1.0/self.memory:
            self.updateScale()

    def reject(self):
        self.jumps.append(0)
        if np.random.rand() < 1.0/self.memory:
            self.updateScale()

    def updateScale(self):
        rate = float(sum(self.jumps))/self.memory
        self.sd *= 1+0.15*np.sqrt(self.memory)*(rate - 0.44)

    def metroCheck(self, log_p0, log_p1, log_q01 = 0, log_q10 = 0):
        """
        By default assumes a symmetric proposal distribution. Otherwise,
        qxy is the probability of x in (q given y).
        
        The metropolis algorithm accepts the proposed jump (x0 -> x1) iff
        
        a ~ U(0,1) and
        a < min(p1/p0 * q01/q10, 1)
        
        Or, in logarithmic form
        
        a < exp(min(log_p1 + log_q01 - log_p0 - log_q10, 0)).
                """
        a = np.random.uniform(low = 0, high = 1)
        paccept = exp(min(log_p1 + log_q01 - log_p0 - log_q10, 0))
        return a < paccept

    def __call__(self):
        x0 = self.x0()
        x1 = np.random.normal(x0, self.sd, 1)
        log_p0 = self.f(x0)
        log_p1 = self.f(x1)
        if self.metroCheck(log_p0, log_p1):
            if not self.frozen:
                self.accept()
            next = x1
        else:
            if not self.frozen:
                self.reject()
            next = x0
        self.updateEnv(next)
        return next

    next = __call__

class MCMC:

    def __init__(self, env):
        self.env = env
        self.samplers = []

    def inEnv(self, f, retain = None):
        return withEnvironment(self.env, retain = retain)(f)

    def gibbs(self, key):
        """
        Returns a stable Gibbs Sampler for the function that
        implements the iterator protocol.
        """
        s = Gibbs(self.env, key)
        # self.samplers.extend(s)
        return s

    def metropolis(self, key):
        s = Metropolis(self.env, key)
        self.samplers.extend(s)
        return s

    def updateAt(self, f, key):
        """
        Call the function f and update the environment at key with the
        output from the function.
        """
        def wrapped():
            val = f.next()
            self.env[key] = val
            return val
        return Sampler(wrapped)



if __name__ == "__main__":
    import doctest
    doctest.testmod()
