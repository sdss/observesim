import numpy as np
import fitsio
import observesim.cadence
import sys


# Class to define a singleton
class FieldsSingleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(FieldsSingleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Fields(object, metaclass=FieldsSingleton):
    """List of fields

    Parameters:
    ----------

    Attributes:
    ----------

    nfields : np.int32
        number of fields

    racen : ndarray of np.float64
        right ascension of each design (J2000 deg)

    deccen : ndarray of np.float64
        declination of each design (J2000 deg)

    fieldid : ndarray of np.int32
        id for each field

    cadence : string
        cadence for each field

    observations : list of ndarray of np.int32
        for each field, indices of observations
"""
    def __init__(self):
        self.nfields = 0
        self.racen = np.zeros(0, dtype=np.float64)
        self.deccen = np.zeros(0, dtype=np.float64)
        self.fieldid = np.zeros(0, dtype=np.int32)
        self.nextmjd = np.zeros(0, dtype=np.float64)
        self.cadence = []
        self.observations = []
        self.cadencelist = observesim.cadence.CadenceList()
        return

    def fromarray(self, fields_array=None):
        self.nfields = len(fields_array)
        self.racen = fields_array['racen']
        self.deccen = fields_array['deccen']
        self.fieldid = np.arange(self.nfields, dtype=np.int32)
        self.cadence = [c.decode().strip() for c in fields_array['cadence']]
        self.observations = [np.zeros(0, dtype=np.int32)] * self.nfields
        self.icadence = np.zeros(self.nfields, dtype=np.int32)
        self.nextmjd = np.zeros(self.nfields, dtype=np.float64)
        return

    def fromfits(self, filename=None):
        self.fields_fits = fitsio.read(filename)
        self.fromarray(self.fields_fits)
        return

    def add_observations(self, mjd=None, fieldid=None, iobs=None):
        self.observations[fieldid] = np.append(self.observations[fieldid],
                                               iobs)
        self.icadence[fieldid] = self.icadence[fieldid] + 1
        cadence = self.cadencelist.cadences[self.cadence[fieldid]]
        if(self.icadence[fieldid] < cadence.nexposures):
            self.nextmjd[fieldid] = (mjd +
                                     cadence.delta[self.icadence[fieldid]] -
                                     cadence.softness[self.icadence[fieldid]])
        else:
            self.nextmjd[fieldid] = 100000.
        return

    def toarray(self):
        """Return cadences as a record array

        Returns:
        -------

        fields : ndarray
            information on each field
"""
        maxn = np.array([len(x) for x in self.observations]).max()
        if(maxn == 1):
            maxn = 2
        fields0 = [('fieldid', np.int32),
                   ('racen', np.float64),
                   ('deccen', np.float64),
                   ('cadence', np.dtype('a20')),
                   ('nobservations', np.int32),
                   ('observations', np.int32, maxn)]
        fields = np.zeros(self.nfields, dtype=fields0)
        fields['fieldid'] = self.fieldid
        fields['racen'] = self.racen
        fields['deccen'] = self.deccen
        for indx in np.arange(self.nfields):
            fields['cadence'][indx] = self.cadence[indx]
            fields['nobservations'][indx] = len(self.observations[indx])
            if(fields['nobservations'][indx] > 0):
                fields['observations'][indx][0:fields['nobservations'][indx]] = self.observations[indx]
        return(fields)
