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

    sn2 : ndarray of np.float32
        duration of observation (days)

    airmass : ndarray of np.float32
        airmass of observation

    lunation : ndarray of np.float32
        lunation of observation (Moon illumination fraction, or zero if
        below horizon)

    lst : ndarray of np.float32
        LST of observation

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
        self.sn2 = np.zeros(0, dtype=np.float32)
        self.airmass = np.zeros(0, dtype=np.float32)
        self.lunation = np.zeros(0, dtype=np.float32)
        self.lst = np.zeros(0, dtype=np.float32)
        return

    def add(self, fieldid=None, mjd=None, duration=None, sn2=None,
            lunation=None, airmass=None, lst=None):
        self.fieldid = np.append(self.fieldid,
                                 np.array([np.float64(fieldid)]))
        self.mjd = np.append(self.mjd,
                             np.array([np.float64(mjd)]))
        self.duration = np.append(self.duration,
                                  np.array([np.float64(duration)]))
        self.sn2 = np.append(self.sn2,
                             np.array([np.float32(sn2)]))
        self.lunation = np.append(self.lunation,
                                  np.array([np.float32(lunation)]))
        self.airmass = np.append(self.airmass,
                                 np.array([np.float32(airmass)]))
        self.lst = np.append(self.lst,
                             np.array([np.float32(lst)]))
        self.nobservations = len(self.fieldid)
        return(self.nobservations - 1)

    def forfield(self, mjd=None, fieldid=None):
        indx = np.where((self.mjd <= mjd) &
                        (self.fieldid == fieldid))[0]
        return(self.toarray(indx=indx))

    def toarray(self, indx=None):
        """Return observations as a record array

        Parameters:
        ----------

        indx : ndarray of np.int32
            indices of observations to return (default to all)

        Returns:
        -------

        observations : record array
            observation information
"""
        obs0 = [('fieldid', np.int32),
                ('mjd', np.float64),
                ('duration', np.float64),
                ('sn2', np.float32),
                ('airmass', np.float32),
                ('lunation', np.float32),
                ('lst', np.float32)]
        if(indx is None):
            indx = np.arange(self.nobservations)
        nobs = len(indx)
        obs = np.zeros(nobs, dtype=obs0)
        if(nobs > 0):
            obs['fieldid'] = self.fieldid[indx]
            obs['mjd'] = self.mjd[indx]
            obs['duration'] = self.duration[indx]
            obs['sn2'] = self.sn2[indx]
            obs['airmass'] = self.airmass[indx]
            obs['lunation'] = self.lunation[indx]
            obs['lst'] = self.lst[indx]
        return(obs)
