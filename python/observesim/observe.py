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
    def __init__(self):
        self._result_dtype = [('sn2', np.float32),
                              ('mjd', np.float64),
                              ('duration', np.float64)]
        self._result0 = np.zeros(1, dtype=self._result_dtype)
        pass

    def result(self, tileid=None, duration=None, mjd=None):
        """Return simulated result of an observation

        Parameters:
        ----------

        tileid : int, np.int32
            id of tile to observe

        duration : float, np.float64
            duration of observation (days)

        mjd : float, np.float64
            MJD start time of observation (days)

        Returns:
        -------

        result : dict
            has keys 'sn2', 'mjd', 'duration'
"""
        fresult = self._result0
        fresult['sn2'] = 3000.
        fresult['mjd'] = mjd
        fresult['duration'] = duration
        return(fresult)
