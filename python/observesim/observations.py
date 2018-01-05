import numpy as np

"""Observations module class.

Dependencies:

 numpy
 scipy

"""


class Observations(object):
    """Observations class.

    Parameters:
    ----------

    observatory : string
       Observatory for fields ('apo' or 'lco'; default 'apo')

    Attributes:
    ----------

    nobservations : np.int32
        number of observations

    fieldid : ndarray of np.int32
        id of each field for observations

    mjd : ndarray of np.float64
        MJD of observation (days)

    duration : ndarray of np.float64
        duration of observation (days)

    sn2 : ndarray of np.float64
        duration of observation (days)

    Methods:
    -------

    add() : add an observation of a field 
    toarray() : Return ndarray of field properties

"""
    def __init__(self, observatory='apo'):
        self.nobservations = np.int32(0)
        self.observatory = observatory
        self.fieldid = np.zeros(0, dtype=np.int32)
        self.mjd = np.zeros(0, dtype=np.float64)
        self.duration = np.zeros(0, dtype=np.float64)
        self.sn2 = np.zeros(0, dtype=np.float64)
        return

    def add(self, fieldid=None, mjd=None, duration=None, sn2=None):
        self.fieldid = np.append(self.fieldid,
                                 np.array([np.float64(fieldid)]))
        self.mjd = np.append(self.mjd,
                             np.array([np.float64(mjd)]))
        self.duration = np.append(self.duration,
                                  np.array([np.float64(duration)]))
        self.sn2 = np.append(self.sn2,
                             np.array([np.float64(sn2)]))
        self.nobservations = len(self.fieldid)
        return

    def forfield(self, mjd=None, fieldid=None):
        indx = np.where((self.mjd <= mjd) &
                        (self.fieldid == fieldid))[0]
        return(self.toarray(indx=indx))

    def toarray(self, indx=None):
        obs0 = [('fieldid', np.int32),
                ('mjd', np.float64),
                ('duration', np.float64),
                ('sn2', np.float64)]
        if(indx is None):
            indx = np.arange(self.nobservations)
        nobs = len(indx)
        obs = np.zeros(nobs, dtype=obs0)
        if(nobs > 0):
            obs['fieldid'] = self.fieldid[indx]
            obs['mjd'] = self.mjd[indx]
            obs['duration'] = self.duration[indx]
            obs['sn2'] = self.sn2[indx]
        return(obs)
