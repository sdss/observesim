import os
import numpy as np
import numpy.random as random
import astropy.io.ascii as ascii
import observesim.db.peewee.targetdb as targetdb

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
    def __init__(self, db=True, fps_layout='central_park'):
        if(db):
            self._read_db(fps_layout=fps_layout)
        else:
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
        self.positionerid = np.arange(len(rfp), dtype=np.int32) + 1
        self.indx = dict()
        for (indx, pid) in zip(np.arange(len(rfp)), self.positionerid):
            self.indx[pid] = indx
        self.row = rfp['row']
        self.pos = rfp['pos']
        self.xcen = rfp['xcen']
        self.ycen = rfp['ycen']
        self.assignment = rfp['assignment']
        self.optical = ((self.assignment == "BOSS") |
                        (self.assignment == "BA"))
        self.apogee = (self.assignment == "BA")
        return

    def _read_db(self, fps_layout='central_park'):
        """Read db and set settings (not functional yet)"""
        targetdb.database.connect_from_config('local')

        nactuators = targetdb.Fiber.select(targetdb.Actuator).count()
        actuators = (targetdb.Actuator.select()
                     .order_by(targetdb.Actuator.id)
                     .join(targetdb.FPSLayout,
                           on=(targetdb.Actuator.fps_layout_pk == targetdb.FPSLayout.pk))
                     .where(targetdb.FPSLayout.label == fps_layout)
                     .dicts())
        fibers = (targetdb.Fiber.select(targetdb.Fiber.fiberid,
                                        targetdb.Spectrograph.label.alias('spectrograph'),
                                        targetdb.Actuator.id,
                                        targetdb.FPSLayout.label.alias('fps_layout'))
                  .join(targetdb.Spectrograph,
                        on=(targetdb.Fiber.spectrograph_pk == targetdb.Spectrograph.pk))
                  .join(targetdb.Actuator,
                        on=(targetdb.Fiber.actuator_pk == targetdb.Actuator.pk))
                  .join(targetdb.FPSLayout,
                        on=(targetdb.Actuator.fps_layout_pk == targetdb.FPSLayout.pk))
                  .where(targetdb.FPSLayout.label == fps_layout)
                 ).dicts()

        self.positionerid = np.zeros(nactuators, dtype=np.int32)
        self.xcen = np.zeros(nactuators, dtype=np.float32)
        self.ycen = np.zeros(nactuators, dtype=np.float32)
        self.optical = np.zeros(nactuators, dtype=np.bool)
        self.apogee = np.zeros(nactuators, dtype=np.bool)
        self.indx = dict()

        indx = 0
        for indx, actuator in zip(np.arange(nactuators), actuators):
            self.positionerid[indx] = actuator['id']
            self.indx[self.positionerid[indx]] = indx
            self.xcen[indx] = actuator['xcen']
            self.ycen[indx] = actuator['ycen']

        for fiber in fibers:
            if(fiber['spectrograph'] == 'APOGEE'):
                self.apogee[self.indx[fiber['id']]] = True
            if(fiber['spectrograph'] == 'BOSS'):
                self.optical[self.indx[fiber['id']]] = True

        return

    def positioners(self, x=None, y=None):
        distances = np.sqrt((x - self.xcen)**2 + (y - self.ycen)**2)
        imatch = np.where((distances > self.inner_reach) &
                          (distances < self.outer_reach))[0]
        return self.positionerid[imatch]

    def targets(self, positionerid=None, x=None, y=None):
        xcen = self.xcen[self.indx[positionerid]]
        ycen = self.ycen[self.indx[positionerid]]
        distances = np.sqrt((x - xcen)**2 + (y - ycen)**2)
        imatch = np.where((distances > self.inner_reach) &
                          (distances < self.outer_reach))[0]
        return imatch

    def covered(self, x=None, y=None, type=None):
        if(type is 'optical'):
            iposid = np.where(self.optical > 0)[0]
        else:
            iposid = np.where(self.apogee > 0)[0]
        covered = np.zeros(len(x), dtype=np.int32)
        for (indx, positionerid) in zip(np.arange(len(iposid)),
                                        self.positionerid[iposid]):
            imatch = self.targets(positionerid=positionerid, x=x, y=y)
            covered[imatch] = 1
        return(covered)

    def assign(self, x=None, y=None, type=None):
        if(type is 'optical'):
            iposid = np.where(self.optical > 0)[0]
        else:
            iposid = np.where(self.apogee > 0)[0]
        positionerids = np.zeros(len(x), dtype=np.int32) - 1
        targets = np.zeros(len(self.positionerid), dtype=np.int32) - 1
        for (indx, positionerid) in zip(np.arange(len(iposid)),
                                        self.positionerid[iposid]):
            ifree = np.where(positionerids == -1)[0]
            imatch = self.targets(positionerid=positionerid,
                                  x=x[ifree], y=y[ifree])
            if(len(imatch)) > 0:
                random.shuffle(imatch)
                positionerids[ifree[imatch[0]]] = positionerid
                targets[iposid[indx]] = ifree[imatch[0]]
        return(positionerids, targets)
