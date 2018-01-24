import numpy as np
import observesim.sloane as sloane

"""Fields module class.

Dependencies:

 numpy
 scipy
 observesim.sloane

"""


class Fields(object):
    """Fields class.

    Parameters:
    ----------

    observatory : string
       Observatory for fields ('apo' or 'lco'; default 'apo')

    nsloane : int, np.int32
       Number of Sloane fields to use as a base

    Attributes:
    ----------

    nfields : np.int32
        number of fields

    racen : ndarray of np.float64
        right ascension of each field

    deccen : ndarray of np.float64
        declination of each field

    fieldtype : list of strings
        type of each field

    fieldid : ndarray of np.int32
        id for each field

    priority : ndarray of np.int32
        priority number (smaller is higher priority)

    Methods:
    -------

    toarray() : Return ndarray of field properties

    Comments:
    --------

    In this implementation, the attributes that are lists or arrays
    are arranged so that attribute[fieldid] is the attribute for field
    fieldid.

    """
    def __init__(self, observatory='apo', nsloane=14522):
        self._limit = 9999
        self._default_duration = 15. / 60. / 24.  # in days
        self._set_fields(observatory=observatory, nsloane=nsloane)
        return

    def _set_fields(self, observatory='apo', nsloane=None):
        sl = sloane.Sloane(n=nsloane)
        (ra, dec) = (sl.ra, sl.dec)
        if(observatory == 'apo'):
            indx = np.where(dec >= -10.)[0]
            ra = ra[indx]
            dec = dec[indx]
        elif(observatory == 'lco'):
            indx = np.where(dec < -10.)[0]
            ra = ra[indx]
            dec = dec[indx]
        else:
            return
        self.nfields = len(ra)
        self.fieldid = np.arange(self.nfields, dtype=np.int32)
        self.racen = ra
        self.deccen = dec
        self.priority = np.zeros(self.nfields, dtype=np.int32) + self._limit
        self.duration = (np.zeros(self.nfields, dtype=np.float64) +
                         self._default_duration)
        self.fieldtype = list('standard' for indx in range(self.nfields))
        self.observatory = observatory
        return

    def toarray(self, indx=None):
        """Return fields as a record array

        Parameters:
        ----------

        indx : ndarray of np.int32
            fieldids to return (default to all)

        Returns:
        -------

        fields : record array
            field information
"""
        field0 = [('fieldid', np.int32),
                  ('racen', np.float64),
                  ('deccen', np.float64),
                  ('duration', np.float64)]
        if(indx is None):
            indx = np.arange(self.nfields)
        nfields = len(indx)
        fields = np.zeros(nfields, dtype=field0)
        fields['fieldid'] = self.fieldid[indx]
        fields['racen'] = self.racen[indx]
        fields['deccen'] = self.deccen[indx]
        fields['duration'] = self.duration[indx]
        return(fields)
