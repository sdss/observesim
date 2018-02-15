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


class TileBase(object):
    """Tile base class.

    Parameters:
    ----------

    observatory : string
       Observatory for tiles ('apo' or 'lco'; default 'apo')

    nsloane : int, np.int32
       Number of Sloane tiles to use as a base

    Attributes:
    ----------

    ntiles : np.int32
        number of tiles

    racen : ndarray of np.float64
        right ascension of each tile

    deccen : ndarray of np.float64
        declination of each tile

    tiletype : list of strings
        type of each tile

    tileid : ndarray of np.int32
        id for each tile

    priority : ndarray of np.int32
        priority number (smaller is higher priority)

    Methods:
    -------

    toarray() : Return ndarray of tile properties

    Comments:
    --------

    In this implementation, the attributes that are lists or arrays
    are arranged so that attribute[tileid] is the attribute for tile
    tileid.

    Uses a uniform distribution of tiles across the sky.

    """
    def __init__(self, observatory='apo', nsloane=8192):
        self._set_settings()
        self._set_tiles(observatory=observatory, nsloane=nsloane)
        return

    def _set_settings(self):
        self._limit = 9999
        self._default_duration = 15. / 60. / 24.  # in days
        return

    def _set_tiles(self, observatory='apo', nsloane=None):
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
        self.ntiles = len(ra)
        self.tileid = np.arange(self.ntiles, dtype=np.int32)
        self.racen = ra
        self.deccen = dec
        self.priority = np.zeros(self.ntiles, dtype=np.int32) + self._limit
        self.duration = (np.zeros(self.ntiles, dtype=np.float64) +
                         self._default_duration)
        self.tiletype = list('standard' for indx in range(self.ntiles))
        self.observatory = observatory
        return

    def toarray(self, indx=None):
        """Return tiles as a record array

        Parameters:
        ----------

        indx : ndarray of np.int32
            tileids to return (default to all)

        Returns:
        -------

        tiles : record array
            tile information
"""
        tile0 = [('tileid', np.int32),
                 ('tiletype', 'O'),
                 ('racen', np.float64),
                 ('deccen', np.float64),
                 ('duration', np.float64)]
        if(indx is None):
            indx = np.arange(self.ntiles)
        ntiles = len(indx)
        tiles = np.zeros(ntiles, dtype=tile0)
        tiles['tileid'] = self.tileid[indx]
        tiles['tiletype'] = self.tiletype[indx]
        tiles['racen'] = self.racen[indx]
        tiles['deccen'] = self.deccen[indx]
        tiles['duration'] = self.duration[indx]
        return(tiles)


class TileFile(TileBase):
    """Tile version reading from files

    Parameters:
    ----------

    observatory : string
       Observatory for tiles ('apo' or 'lco'; default 'apo')

    filebase : int, np.int32
       Base name of files to read

    Attributes:
    ----------

    ntiles : np.int32
        number of tiles

    racen : ndarray of np.float64
        right ascension of each tile

    deccen : ndarray of np.float64
        declination of each tile

    tiletype : list of strings
        type of each tile

    tileid : ndarray of np.int32
        id for each tile

    priority : ndarray of np.int32
        priority number (smaller is higher priority)

    Methods:
    -------

    toarray() : Return ndarray of tile properties

    Comments:
    --------

    In this implementation, the attributes that are lists or arrays
    are arranged so that attribute[tileid] is the attribute for tile
    tileid.

    Reads tiles from files [filebase]-apo.fits and [filebase]-lco.fits

    """
    def __init__(self, observatory='apo', filebase=None):
        self._set_settings()
        self._set_tiles(observatory=observatory, filebase=filebase)
        return

    def _set_tiles(self, observatory='apo', filebase=None):
        filename = '{filebase}-{obs}.fits'.format(filebase=filebase,
                                                  obs=observatory)
        self._tiles = fitsio.read(filename)
        self.ntiles = len(self._tiles)
        self.tileid = np.arange(self.ntiles, dtype=np.int32)
        self.tiletype = self._tiles['tiletype']
        self.racen = self._tiles['racen']
        self.deccen = self._tiles['deccen']
        self.priority = np.zeros(self.ntiles, dtype=np.int32) + self._limit
        self.duration = (np.zeros(self.ntiles, dtype=np.float64) +
                         self._default_duration)
        self.observatory = observatory
        return
