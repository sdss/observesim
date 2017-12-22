import numpy as np

"""Scheduler module class.

Dependencies:

 numpy

"""


class Observe(object):
    """Observe class.

    Methods:
    -------

"""
    def __init__(self):
        self._result_dtype = [('sn2', np.float32),
                              ('mjd', np.float64),
                              ('duration', np.float64)]
        self._result0 = np.zeros(1, dtype=self._result_dtype)
        pass

    def result(self, fieldid=None, duration=None, mjd=None):
        fresult = self._result0
        fresult['sn2'] = 3000.
        fresult['mjd'] = mjd
        fresult['duration'] = duration
        return(fresult)
