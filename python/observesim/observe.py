import numpy as np

"""Scheduler module class.

Dependencies:

 numpy

"""


class Observe(object):
    """Observe class.

    Used to simulate observations of tiles.

    Methods:
    -------

    result() : return result of an observation

"""
    def __init__(self, defaultExp=None, cadencelist=None, cadences=None):
        """Return simulated result of an observation

        Parameters:
        ----------

        defaultExp : float
            default exposure length in days

        cadencelist : 
            the cadencelist object from scheduler

        cadences : 
            the cadences from scheduler (of length fields, matching that index)
        """

        self._result_dtype = [('sn2', np.float32),
                              ('mjd', np.float64),
                              ('duration', np.float64)]
        self._result0 = np.zeros(1, dtype=self._result_dtype)
        self.cadencelist = cadencelist
        self.cadences = cadences
        if defaultExp is None:
            self.defaultExp = np.float32(15. / 60. / 24.)
        else:
            self.defaultExp = defaultExp
        self.defaultSN2 = 3000.
        pass

    def result(self, fieldid=None, mjd=None, airmass=1, **kwargs):
        """Return simulated result of an observation

        Parameters:
        ----------

        fieldid : int, np.int32
            id of field to observe

        mjd : float, np.float64
            MJD start time of observation (days)

        airmass: float, np.float64
            The airmass of the observation

        **kwargs
            catch old duration kwarg quietly?

        Returns:
        -------

        result : dict
            has keys 'sn2', 'mjd', 'duration'
"""
        
        if self.cadences is None or self.cadencelist is None:
            duration = self.defaultExp
            sn2 = self.defaultSN2
        else:
            boss = self.cadencelist[self.cadences[fieldid]].requires_boss
            if boss:
                weight = 1
            else:
                weight = 0.05
            duration = self.defaultExp * airmass ** weight
            sn2 = self.defaultSN2  # + (np.random.randn() * 1000)

        fresult = self._result0
        fresult['sn2'] = sn2
        fresult['mjd'] = mjd
        fresult['duration'] = duration
        return(fresult)
