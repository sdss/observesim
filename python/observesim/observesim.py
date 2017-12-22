import numpy as np

"""Scheduler module class.

Dependencies:

 numpy

"""


class ObserveSim(object):
    """ObservSim class.

    Parameters:
    ----------

    Methods:
    -------

"""
    def __init__(self):
        self._result_dtype = [('sn2', np.float32)]
        self._result0 = np.zeros(1, dtype=self._result_dtype)
        pass

    def result(self, field=None):
        fresult = self._result0
        fresult['sn2'] = 10.
        return(fresult)
