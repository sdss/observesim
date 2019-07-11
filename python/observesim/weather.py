import numpy as np
import numpy.fft as fft
import scipy.interpolate as interpolate

"""Weather module for simulations.

Generates simulated weather outcomes for a range of MJD.

Dependencies:

 numpy
 scipy

"""


class Weather(object):
    """Weather class

    Parameters:
    ----------

    mjd_start : float, np.float64
        Starting MJD to consider

    mjd_end : float, np.float64
        Ending MJD to consider

    dmjd_minutes : float, np.float64
        smallest resolution in time, in minutes (default 10)

    sigma : float, np.float64
        Gaussian cutoff of power spectrum in hours (default 2)

    alpha : float, np.float64
        Power law slope of power spectrum (default -0.75)

    fclear : float, np.float64
        fraction of clear time (default 0.5)

    seed : int
        random seed

    Methods:
    -------

    clear(mjd) : is it clear for this MJD, and how long

"""
    def __init__(self, mjd_start=None, mjd_end=None, dmjd_minutes=10.,
                 sigma=2., alpha=-0.75, fclear=0.5, seed=1):
        self.mjd_start = mjd_start
        self.mjd_end = mjd_end
        self.alpha = alpha
        self.sigma = sigma
        self.fclear = fclear
        self.dmjd = dmjd_minutes / 60. / 24.
        self.nmjd = np.int32(np.ceil((self.mjd_end - self.mjd_start) /
                                     self.dmjd))
        self.mjd = self.mjd_start + np.arange(self.nmjd) * self.dmjd
        self._initialize_conditions(seed=seed)

    def _initialize_conditions(self, seed=1):
        """Initialize the pattern of clear weather."""
        if seed is not None:
            print("!!using psuedo-random weather!!")
            np.random.seed(seed)
        nsigma = self.sigma / self.dmjd
        psigma = 1. / nsigma
        pk = np.zeros(self.nmjd // 2)
        pk[0] = 0.
        pk[1:len(pk)] = ((np.arange(len(pk) - 1) + 0.5)**self.alpha *
                         np.exp(- 0.5 * ((np.arange(len(pk) - 1) + 1.) /
                                         np.float64(len(pk)))**2 /
                                psigma**2))
        pk[1:len(pk)] = pk[1:len(pk)] - pk[len(pk) - 1]
        ufft = np.zeros(self.nmjd, dtype=np.complex64)
        ufft[0:len(pk):1] = (np.random.normal(size=len(pk)) +
                             1j * np.random.normal(size=len(pk))) * np.sqrt(pk)
        ufft[-1:-len(pk):-1] = np.conj(ufft[1:len(pk):1])
        udist = fft.ifft(ufft).real
        isort = np.argsort(udist)
        self._uvals = np.zeros(self.nmjd)
        self._uvals[isort] = ((np.arange(self.nmjd) + 0.5) /
                              np.float64(self.nmjd))
        self.clear_pattern = interpolate.interp1d(self.mjd, self._uvals,
                                                  fill_value='extrapolate')

    def clear(self, mjd=None, returnNext=True):
        """For a given MJD, return if it is clear and when next change is.

        Parameters:
        ----------

        mjd : float, np.float64
            MJD to check (days)

        returnNext : boolean
            Also return MJD of next change

        Returns:
        -------

        isclear : boolean
            Is the mjd clear?

        nextchange : float, np.float64
            MJD of next change of state

        Comments:
        --------
        Not high performance. Only takes a single MJD.
"""
        isclear = (self.clear_pattern(mjd) < self.fclear)
        if(returnNext is False):
            return(isclear)
        step = 5.
        ncheck = np.int32(np.ceil(step / (0.5 * self.dmjd)))
        dstep = ncheck * self.dmjd
        mjd_base = mjd
        while((mjd_base < self.mjd_end)):
            check_mjds = (mjd_base + self.dmjd * (np.arange(ncheck) + 1.))
            isclear_mjds = (self.clear_pattern(check_mjds) < self.fclear)
            different = np.where(isclear_mjds != isclear)[0]
            if(len(different) > 0):
                return(isclear, check_mjds[different[0]])
            mjd_base = mjd_base + dstep
        return(isclear, self.mjd_end)
