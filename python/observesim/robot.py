import os
import numpy as np
import astropy.io.ascii as ascii

"""Robot module class.

Dependencies:

 numpy
 astropy

"""


class Robot(object):
    """Robot class

    Parameters:
    ----------

    Attributes:
    ----------

    npositioner : np.int32
        number of positioners

    positionerid : ndarray of np.int32
        index number of positioners

    row : ndarray of np.int32
        row numbers

    pos : ndarray of np.int32
        position numbers

    xcen : ndarray of np.float64
        xfocal center position of fiber positioners (mm)

    ycen : ndarray of np.float64
        yfocal center position of fiber positioners (mm)

    assignment : list of strings
        assignment name

    optical : ndarray of np.boolean
        positioner has an optical fiber

    apogee : ndarray of np.boolean
        positioner has an APOGEE fiber

    inner_reach : np.float64
        inner annulus reach (mm)

    outer_reach : np.float64
        outer annulus reach (mm)

    Methods:
    -------

    positioners() : which positioners can reach a given x, y

    Notes:
    -----

    Some positions may be fiducials (neither optical nor apogee).
"""
    def __init__(self):
        self._read_config()
        self._set_parameters()
        return

    def _set_parameters(self):
        """Set basic paramaters"""
        self._ralpha = 7.4  # alpha arm radius in mm
        self._rbeta = 15.0  # beta arm radius in mm
        self._pitch = self._ralpha + self._rbeta
        self.inner_reach = self._rbeta - self._ralpha
        self.outer_reach = self._ralpha + self._rbeta
        self.exclusion = 11.

    def _read_config(self):
        """Read config file and set settings"""
        rfpfile = os.path.join(os.getenv("OBSERVESIM_DIR"),
                               "data", "fps_RTConfig.txt")
        rfp = ascii.read(rfpfile)
        self.positionerid = np.arange(len(rfp), dtype=np.int32)
        self.row = rfp['row']
        self.pos = rfp['pos']
        self.xcen = rfp['xcen']
        self.ycen = rfp['ycen']
        self.assignment = rfp['assignment']
        self.optical = ((self.assignment == "BOSS") |
                        (self.assignment == "BA"))
        self.apogee = (self.assignment == "BA")
        return

    def positioners(self, x=None, y=None):
        distances = np.sqrt((x - self.xcen)**2 + (y - self.ycen)**2)
        imatch = np.where((distances > self.inner_reach) &
                          (distances < self.outer_reach))[0]
        return self.positionerid[imatch]
