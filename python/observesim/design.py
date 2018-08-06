import os
import numpy as np
import fitsio
import observesim.sloane as sloane

"""Fields module class.

Dependencies:

 numpy
 scipy
 observesim.sloane

"""


class DesignBase(object):
    """Design base class.

    Parameters:
    ----------

    observatory : string
       Observatory for designs ('apo' or 'lco'; default 'apo')

    nsloane : int, np.int32
       Number of Sloane designs to use as a base

    Attributes:
    ----------

    ndesigns : np.int32
        number of designs

    racen : ndarray of np.float64
        right ascension of each design

    deccen : ndarray of np.float64
        declination of each design

    designtype : list of strings
        type of each design

    designid : ndarray of np.int32
        id for each design

    priority : ndarray of np.int32
        priority number (smaller is higher priority)

    Methods:
    -------

    toarray() : Return ndarray of design properties

    Comments:
    --------

    In this implementation, the attributes that are lists or arrays
    are arranged so that attribute[designid] is the attribute for design
    designid.

    Uses a uniform distribution of designs across the sky.

    """
    def __init__(self, observatory='apo', nsloane=8192):
        self._set_settings()
        self._set_designs(observatory=observatory, nsloane=nsloane)
        return

    def _set_settings(self):
        self._limit = 9999
        self._default_duration = 15. / 60. / 24.  # in days
        return

    def _set_designs(self, observatory='apo', nsloane=None):
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
        self.ndesigns = len(ra)
        self.designid = np.arange(self.ndesigns, dtype=np.int32)
        self.racen = ra
        self.deccen = dec
        self.priority = np.zeros(self.ndesigns, dtype=np.int32) + self._limit
        self.duration = (np.zeros(self.ndesigns, dtype=np.float64) +
                         self._default_duration)
        self.designtype = list('standard' for indx in range(self.ndesigns))
        self.observatory = observatory
        return

    def toarray(self, indx=None):
        """Return designs as a record array

        Parameters:
        ----------

        indx : ndarray of np.int32
            designids to return (default to all)

        Returns:
        -------

        designs : record array
            design information
"""
        design0 = [('designid', np.int32),
                   ('designtype', 'O'),
                   ('racen', np.float64),
                   ('deccen', np.float64),
                   ('duration', np.float64)]
        if(indx is None):
            indx = np.arange(self.ndesigns)
        ndesigns = len(indx)
        designs = np.zeros(ndesigns, dtype=design0)
        designs['designid'] = self.designid[indx]
        designs['designtype'] = self.designtype[indx]
        designs['racen'] = self.racen[indx]
        designs['deccen'] = self.deccen[indx]
        designs['duration'] = self.duration[indx]
        return(designs)


class DesignFile(DesignBase):
    """Design version reading from files

    Parameters:
    ----------

    observatory : string
       Observatory for designs ('apo' or 'lco'; default 'apo')

    filebase : int, np.int32
       Base name of files to read

    Attributes:
    ----------

    ndesigns : np.int32
        number of designs

    racen : ndarray of np.float64
        right ascension of each design

    deccen : ndarray of np.float64
        declination of each design

    designtype : list of strings
        type of each design

    designid : ndarray of np.int32
        id for each design

    lunation : ndarray of np.float32
        maximum illumination of Moon to observe at

    priority : ndarray of np.int32
        priority number (smaller is higher priority)

    Methods:
    -------

    toarray() : Return ndarray of design properties

    Comments:
    --------

    In this implementation, the attributes that are lists or arrays
    are arranged so that attribute[designid] is the attribute for design
    designid.

    Reads designs from files [filebase]-apo.fits and [filebase]-lco.fits

    """
    def __init__(self, observatory='apo', filebase=None):
        self._set_settings()
        self._set_designs(observatory=observatory, filebase=filebase)
        return

    def _set_designs(self, observatory='apo', filebase=None):
        filename = '{filebase}-{obs}.fits'.format(filebase=filebase,
                                                  obs=observatory)
        self._designs = fitsio.read(filename)
        self.ndesigns = len(self._designs)
        self.designid = np.arange(self.ndesigns, dtype=np.int32)
        self.designtype = self._designs['TILETYPE']
        self.lunation = self._designs['LUNATION']
        self.racen = self._designs['RA']
        self.deccen = self._designs['DEC']
        self.priority = np.zeros(self.ndesigns, dtype=np.int32) + self._limit
        self.duration = (np.zeros(self.ndesigns, dtype=np.float64) +
                         self._default_duration * self._designs['NEXP'])
        self.observatory = observatory
        return
