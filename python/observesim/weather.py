import os

import numpy as np
import numpy.fft as fft
import scipy.interpolate as interpolate


"""Weather module for simulations.

Generates simulated weather outcomes for a range of MJD.

Dependencies:

 numpy
 scipy

"""


class Weather3(object):
    """Weather class

    Parameters:
    ----------

    model_fname: str, pathlike
        path to saved model

    mjd_start : float, np.float64
        Starting MJD to consider

    mjd_end : float, np.float64
        Ending MJD to consider

    seed : int
        random seed

    Methods:
    -------

    clear() : is it clear for the current MJD, and how long
              until next change?
    """
    def __init__(self, model_fname=None,
                 mjd_start=None, mjd_end=None):
        self.mjd_start = mjd_start
        self.mjd_end = mjd_end

        self.model = np.genfromtxt(model_fname,
                                   dtype=None,
                                   delimiter=",",
                                   names=True,
                                   encoding="UTF8")

        # self.idx = np.argmin(np.abs(self.model["mjd"]-self.mjd_start))
        self.idx = 0
        self.mjd = self.model["mjd"][self.idx]
        self.state = self.model["state"][self.idx]

        assert self.mjd_end <= self.model["mjd"][-1]

    def _state_to_clear(self, state):
        """
        Converts a state integer to a clear / not clear boolean.

        Parameters
        ----------

        state : int
            An integer representing the fine-grained weather
            state

        Returns
        -------

        is_clear : boolean
            Whether the state corresponds to clear weather

        Comments
        --------
        The internal state is more fine-grained than a simple
        boolean clear / not clear rating. This function
        translates the internal state to a boolean rating.
        """
        return state < 2

    def _advance_time(self):
        """
        Advances time by one hour, and draws a new state.

        Comments
        --------
        The internal state in an integer, which represents
        the weather in a more fine-grained manner than
        clear / not clear. This state can then be translated
        to a clear / not clear rating.
        """
        self.idx += 1
        self.mjd = self.model["mjd"][self.idx]
        self.state = self.model["state"][self.idx]

    def clear(self, now=None, until=None):
        """
        Returns whether it is currently clear, and the MJD of
        the next change.

        Parameters
        ----------

        now : float
            MJD to start at, advances time as needed

        until : float
            MJD to stop at

        Returns
        -------

        current_clear : boolean
            Is the MJD clear?

        next_change : float
            MJD of the next change in clear status

        Comments
        --------
        Advances the internal state (both cloud cover and time)
        to the next change.
        """

        if now is not None:
            now = float(now)
            time = f"{now:.3f} {self.mjd:.3f} {self.mjd-now:.4f}"
            assert now + 1 / 24 >= self.mjd, "retreiving past weather not supported \n" + time
            while self.mjd < now or self.mjd - now < 1 / 24:
                self._advance_time()

        if until is None:
            # don't go a full day but enough for any long night
            until = self.mjd + 0.7
        current_clear = self._state_to_clear(self.state)
        current_state = int(self.state)
        while self.state == current_state and self.mjd <= until:
            self._advance_time()

        cloudy = self.state > 0

        return current_clear, self.mjd, cloudy


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

    def clear(self, now=None, until=None):
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

        cloudy: boolean
            Always false for this version

        Comments:
        --------
        Not high performance. Only takes a single MJD.
        """
        isclear = (self.clear_pattern(now) < self.fclear)
        # if(returnNext is False):
        #     return(isclear)
        step = 5.
        ncheck = np.int32(np.ceil(step / (0.5 * self.dmjd)))
        dstep = ncheck * self.dmjd
        mjd_base = now
        while((mjd_base < until)):
            check_mjds = (mjd_base + self.dmjd * (np.arange(ncheck) + 1.))
            isclear_mjds = (self.clear_pattern(check_mjds) < self.fclear)
            different = np.where(isclear_mjds != isclear)[0]
            if(len(different) > 0):
                return(isclear, check_mjds[different[0]], False)
            mjd_base = mjd_base + dstep
        # modify to match new fns so 3rd pos arg is "cloudy"
        return(isclear, until, False)
