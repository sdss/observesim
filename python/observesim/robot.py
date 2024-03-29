#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: robot.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import os
import warnings

import astropy.io.ascii as ascii
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.transforms
from matplotlib.patches import Wedge, Polygon
from matplotlib.collections import PatchCollection
import numpy as np
import numpy.random as random

try:
    import shapely.affinity
    import shapely.geometry
except ModuleNotFoundError:
    pass

import sdssdb.peewee.sdss5db.targetdb as targetdb
from observesim.utils import assign_targets_draining, xy2tp


__all__ = ['Robot', 'Configuration']

"""Robot module class.

Dependencies:

 numpy
 astropy
 matplotlib
 shapely

"""


# Class to define a singleton
class RobotSingleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(RobotSingleton,
                                        cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Robot(object, metaclass=RobotSingleton):
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

    boss : ndarray of np.boolean
        positioner has a BOSS fiber

    apogee : ndarray of np.boolean
        positioner has an APOGEE fiber

    inner_reach : np.float64
        inner annulus reach (mm)

    outer_reach : np.float64
        outer annulus reach (mm)

    Methods:
    -------

    reset() : reset the robot
    positioners() : which positioners can reach a given x, y
    targets() : which x, y positions are reachable by a positioner

    Notes:
    -----

    Some positions may be fiducials (neither boss nor apogee).
"""
    def __init__(self, **kwargs):
        self.reset(**kwargs)
        return

    def reset(self, db=True, fps_layout='filled_hex'):
        """Reset robot"""
        if(db):
            self._read_db(fps_layout=fps_layout)
        else:
            self._read_config()
        self._set_parameters()
        return

    def _set_parameters(self):
        """Set basic parameters"""
        self._ralpha = 7.4  # alpha arm radius in mm
        self._rbeta = 15.0  # beta arm radius in mm
        self._pitch = self._ralpha + self._rbeta
        self.inner_reach = self._rbeta - self._ralpha
        self.outer_reach = self._ralpha + self._rbeta
        self.exclusion = 11.
        self.phi_range = 180  # Either 180 or 360, the range of movement of the beta arm.

    def _read_config(self):
        """Read config file and set settings"""
        rfpfile = os.path.join(os.getenv("OBSERVESIM_DIR"),
                               "data", "fps_RTConfig.txt")
        rfp = ascii.read(rfpfile)
        self.npositioner = len(rfp)
        self.positionerid = np.arange(self.npositioner, dtype=np.int32) + 1
        self.indx = dict()
        for (indx, pid) in zip(np.arange(self.npositioner), self.positionerid):
            self.indx[pid] = indx
        self.row = rfp['row']
        self.pos = rfp['pos']
        self.xcen = rfp['xcen']
        self.ycen = rfp['ycen']
        self.assignment = rfp['assignment']
        self.boss = ((self.assignment == "BOSS") |
                     (self.assignment == "BA"))
        self.apogee = (self.assignment == "BA")
        return

    def _read_db(self, fps_layout='central_park'):
        """Read db and set settings (not functional yet)"""

        actuators = (targetdb.Actuator.select()
                     .order_by(targetdb.Actuator.id)
                     .join(targetdb.FPSLayout,
                           on=(targetdb.Actuator.fps_layout_pk == targetdb.FPSLayout.pk))
                     .where(targetdb.FPSLayout.label == fps_layout))

        nactuators = actuators.count()

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

        self.npositioner = nactuators
        self.positionerid = np.zeros(nactuators, dtype=np.int32)
        self.xcen = np.zeros(nactuators, dtype=np.float32)
        self.ycen = np.zeros(nactuators, dtype=np.float32)
        self.boss = np.zeros(nactuators, dtype=np.bool)
        self.apogee = np.zeros(nactuators, dtype=np.bool)
        self.fiducial = np.zeros(nactuators, dtype=np.bool)
        self.indx = dict()

        indx = 0
        for indx, actuator in enumerate(actuators):
            self.positionerid[indx] = actuator.id
            self.indx[self.positionerid[indx]] = indx
            self.xcen[indx] = actuator.xcen
            self.ycen[indx] = actuator.ycen
            self.fiducial[indx] = actuator.actuator_type.label == 'Fiducial'

        for fiber in fibers:
            if(fiber['spectrograph'] == 'APOGEE'):
                self.apogee[self.indx[fiber['id']]] = True
            if(fiber['spectrograph'] == 'BOSS'):
                self.boss[self.indx[fiber['id']]] = True

        return

    def positioner_overlaps(self, positionerid):
        """For a ``positionerid`` returns a list of overlapping positionerids.

        A positioner is consider to overlap with ``positionerid`` if their
        beta arms can tocuh in any way, that is, if their distance is less
        than twice the added lenghts of the alpha and beta arms.

        """

        xcen = self.xcen
        ycen = self.ycen

        pos = np.array([xcen, ycen]).T
        pos[self.fiducial] = np.nan

        idx = self.indx[positionerid]
        this_pos = pos[idx]

        pos_to_index = np.array(list(self.indx.items()), dtype=np.int)

        distance = np.sqrt(np.sum((pos - this_pos)**2, axis=1))

        collides = distance <= (2 * (self._ralpha + self._rbeta))
        collides[idx] = False  # Don't collide with ourselves

        return np.array([pos_to_index[np.where(pos_to_index[:, 1] == ii)][0][0]
                         for ii in np.where(collides)[0]])

    def corners(self):
        rmax = np.sqrt(self.xcen**2 + self.ycen**2).max() + self.outer_reach
        hsqrt3 = np.sqrt(3.) * 0.5
        xcorners = np.array([- rmax, - 0.5 * rmax, 0.5 * rmax, rmax,
                             0.5 * rmax, - 0.5 * rmax, - rmax],
                            dtype=np.float32)
        ycorners = np.array([0., - hsqrt3 * rmax, - hsqrt3 * rmax, 0.,
                             hsqrt3 * rmax, hsqrt3 * rmax, 0.],
                            dtype=np.float32)
        return(xcorners, ycorners)

    def within_corners(self, x=None, y=None):
        xc, yc = self.corners()
        within = np.ones(len(x), dtype=np.int32)
        for indx in np.arange(len(xc) - 1):
            xe = 0.5 * (xc[indx] + xc[indx + 1])
            ye = 0.5 * (yc[indx] + yc[indx + 1])
            d = (xe * x + ye * y) / (xe * xe + ye * ye)
            within[d > 1.] = 0
        return(within)

    def positioners(self, x=None, y=None):
        distances = np.sqrt((x - self.xcen)**2 + (y - self.ycen)**2)
        imatch = np.where((distances > self.inner_reach) &
                          (distances < self.outer_reach))[0]
        return self.positionerid[imatch]

    def targets(self, positionerid=None, x=None, y=None,
                requires_apogee=None, requires_boss=None):
        xcen = self.xcen[self.indx[positionerid]]
        ycen = self.ycen[self.indx[positionerid]]
        distances = np.sqrt((x - xcen)**2 + (y - ycen)**2)
        within = ((distances > self.inner_reach) &
                  (distances < self.outer_reach))
        istype = np.ones(len(x), dtype=np.bool)
        if(requires_apogee is not None):
            if(self.apogee[self.indx[positionerid]] == 0):
                iapogee = np.where(requires_apogee)[0]
                istype[iapogee] = 0
        if(requires_boss is not None):
            if(self.boss[self.indx[positionerid]] == 0):
                iboss = np.where(requires_boss)[0]
                istype[iboss] = 0
        imatch = np.where(within & istype)[0]
        return imatch

    def covered(self, x=None, y=None, type=None):
        if(type is 'boss'):
            iposid = np.where(self.boss > 0)[0]
        else:
            iposid = np.where(self.apogee > 0)[0]
        covered = np.zeros(len(x), dtype=np.int32)
        for (indx, positionerid) in zip(np.arange(len(iposid)),
                                        self.positionerid[iposid]):
            imatch = self.targets(positionerid=positionerid, x=x, y=y)
            covered[imatch] = 1
        return(covered)

    def assign(self, x=None, y=None, type=None):
        if(type is 'boss'):
            iposid = np.where(self.boss > 0)[0]
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

    def plot(self, xoff=0., yoff=0.):
        fig, ax = plt.subplots()

        width = (self.outer_reach - self.inner_reach)

        ax.scatter(self.xcen, self.ycen, s=4)

        indx = np.where(self.apogee)[0]
        patches_apogee = []
        for xcen, ycen in zip(self.xcen[indx], self.ycen[indx]):
            patches_apogee.append(Wedge((xcen, ycen), self.outer_reach, 0, 360,
                                        width=width))
        p = PatchCollection(patches_apogee, alpha=0.3, color='red')
        ax.add_collection(p)

        indx = np.where(self.boss & ~self.apogee)[0]
        patches_boss = []
        for xcen, ycen in zip(self.xcen[indx], self.ycen[indx]):
            patches_boss.append(Wedge((xcen, ycen), self.outer_reach, 0, 360,
                                        width=width))
        p = PatchCollection(patches_boss, alpha=0.3, color='blue')
        ax.add_collection(p)

        plt.xlabel('xfocal (mm)')
        plt.ylabel('xfocal (mm)')

        return


class Configuration(object):
    """Represents a configuration of the robot around a set of targets.

    Parameters:
        robot (`.Robot` object):
            The `.Robot` object describing the FPS and its layout.
        targets (`numpy.ndarray`):
            A ``Nx2`` array of target positions on the focal plane, in mm.
        reassign (bool):
            If ``True``, uses `~observesim.utils.assign_targets_draining`
            to find an optimal allocation of targets and positioners.

    Attributes:
        target_to_positioner (`numpy.ndarray`):
            A 1-d array with the same length as the number of targets. In the
            same order as ``targets``, contains the index of the positioner
            allocated to it.
        theta_phi (`numpy.ndarray`):
            A 2-d array with the ``(theta, phi)`` values for each positioner.
            Values for indices corresponding to fiducials are set to NaN.

    """

    def __init__(self, robot, targets, reassign=False):

        self.robot = robot
        self.targets = targets
        self._polygons = None

        assert isinstance(self.targets, np.ndarray) and self.targets.ndim == 2

        n_real_positioners = np.sum(~self.robot.fiducial)

        if n_real_positioners < self.targets.shape[0]:
            warnings.warn('more targets than positioners. Only the first '
                          f'{n_real_positioners} targets will be allocated', UserWarning)
            self.targets = targets[0:n_real_positioners]

        self.target_to_positioner = np.zeros(self.targets.shape[0], dtype=np.int)

        if reassign:
            positioner_to_targets = assign_targets_draining(robot, targets)
            for positioner in positioner_to_targets:
                pos_targets = positioner_to_targets[positioner]
                assert len(pos_targets) <= 1
                self.target_to_positioner[pos_targets[0]] = positioner
        else:
            n_target = 0
            for positioner in np.where(~self.robot.fiducial)[0]:
                self.target_to_positioner[n_target] = int(positioner)
                n_target += 1

        self.reset_positioners()

    def reset_positioners(self):
        """Resets positioners to folded."""

        n_positioners = len(self.robot.xcen)

        self.theta_phi = np.zeros((n_positioners, 2), dtype=np.float32)
        self.theta_phi[:, 1] = 180

        # Marks fiducials with NaNs
        self.theta_phi[self.robot.fiducial] = np.nan

    def compute(self):
        """Computes a new set of ``(theta, phi)`` for each positioner."""

        xcen = self.robot.xcen[self.target_to_positioner]
        ycen = self.robot.ycen[self.target_to_positioner]

        theta_phi_init = xy2tp(self.targets[:, 0] - xcen, self.targets[:, 1] - ycen)
        theta_phi_init = theta_phi_init[0, :, :]  # Gets the first configuration

        # Resets positioners
        self.reset_positioners()

        for ii in range(self.targets.shape[0]):
            self.theta_phi[self.target_to_positioner[ii], :] = theta_phi_init[ii, :]

    def get_polygons(self, arm_width=4, only_collision=False, fraction=2/3.):
        """Returns a list of `~shapely.MultiPolygon`

        Parameters:
            arm_width (float):
                The width, in mm, of each of the arms.
            only_collision (bool):
                If ``False``, returns a `~shapely.MultiPolygon` with both
                arms (one `~shapely.Polygon` for each arm). If ``True``,
                returns a single `~shapely.Polygon` with the region of the
                beta arm that can collide with other arms.
            fraction (float):
                The fraction of the beta arm that can collide with other
                beta arms.

        """

        if self._polygons is not None and only_collision is False:
            return self._polygons

        # Arm templates. We'll rotate and translate them below.

        alpha_arm = shapely.geometry.Polygon([(0, arm_width / 2.),
                                              (0, -arm_width / 2.),
                                              (self.robot._ralpha, -arm_width / 2.),
                                              (self.robot._ralpha, +arm_width / 2.),
                                              (0, arm_width / 2.)])

        beta_arm = shapely.geometry.Polygon([(0, arm_width / 2.),
                                             (0, -arm_width / 2.),
                                             (self.robot._rbeta, -arm_width / 2.),
                                             (self.robot._rbeta, +arm_width / 2.),
                                             (0, arm_width / 2.)])

        # The part of the beta arm that we consider that can collide with other
        # arms. Defined as the final "fraction" of the arm.
        beta_arm_collision = shapely.geometry.Polygon(
            [(self.robot._rbeta * (1 - fraction), arm_width / 2.),
             (self.robot._rbeta * (1 - fraction), -arm_width / 2.),
             (self.robot._rbeta, -arm_width / 2.),
             (self.robot._rbeta, +arm_width / 2.),
             (self.robot._rbeta * (1 - fraction), arm_width / 2.)])

        # A fiducial is just a shapely Point. This is mostly a placeholder
        # right now since we consider fiducials cannot collide with
        # positioners.
        actuator = shapely.geometry.Point(0, 0)

        polygons = []

        n_positioner = 0
        for theta, phi in self.theta_phi:

            xcen = self.robot.xcen[n_positioner]
            ycen = self.robot.ycen[n_positioner]

            if self.robot.fiducial[n_positioner]:
                polygons.append(shapely.affinity.translate(actuator, xoff=xcen, yoff=ycen))
                n_positioner += 1
                continue

            alpha_arm_end = (xcen + self.robot._ralpha * np.cos(np.radians(theta)),
                             ycen + self.robot._ralpha * np.sin(np.radians(theta)))

            if only_collision is False:

                alpha_arm_t = shapely.affinity.rotate(alpha_arm, theta, origin=(0, 0))
                alpha_arm_t = shapely.affinity.translate(alpha_arm_t, xoff=xcen, yoff=ycen)

                beta_arm_t = shapely.affinity.rotate(beta_arm, theta + phi, origin=(0, 0))
                beta_arm_t = shapely.affinity.translate(
                    beta_arm_t, xoff=alpha_arm_end[0], yoff=alpha_arm_end[1])

                polygons.append(shapely.geometry.MultiPolygon([alpha_arm_t, beta_arm_t]))

            else:

                beta_arm_collision_t = shapely.affinity.rotate(beta_arm_collision,
                                                               theta + phi, origin=(0, 0))
                beta_arm_collision_t = shapely.affinity.translate(
                    beta_arm_collision_t, xoff=alpha_arm_end[0], yoff=alpha_arm_end[1])

                polygons.append(beta_arm_collision_t)

            n_positioner += 1

        # We don't overwrite self._polygons for only_collision
        if only_collision:
            return polygons

        self._polygons = polygons

        return polygons

    def get_collisions(self):
        """Returns a list of positioners that are collided."""

        collisioned = np.zeros(self.theta_phi.shape[0], dtype=np.bool)
        collision_polygons = self.get_polygons(only_collision=True)

        for ii in range(self.theta_phi.shape[0]):
            if self.robot.fiducial[ii]:  # Skips fiducials
                continue
            beta_arm_ii = collision_polygons[ii]
            for jj in range(ii + 1, self.theta_phi.shape[0]):
                if self.robot.fiducial[ii]:
                    continue
                beta_arm_jj = collision_polygons[jj]
                if beta_arm_jj.intersects(beta_arm_ii):
                    collisioned[ii] = True
                    collisioned[jj] = True

        return collisioned

    def plot(self):
        """Plots the configuration.

        Note that the actuator shapes are slightly modified for stylistic
        purposes and some collisions my not seem evident from the plot. The
        collisions are calculated accurately using
        `~Configuration.get_collisions`.

        """

        with plt.style.context(['seaborn-deep']):

            __, ax = plt.subplots()

            arm_width = 4.

            xcen = self.robot.xcen[~self.robot.fiducial]
            ycen = self.robot.ycen[~self.robot.fiducial]

            # Removes fiducials
            theta_phi = self.theta_phi[~self.robot.fiducial]
            collisions = self.get_collisions()[~self.robot.fiducial]

            # Calculates the positions of the articulation between alpha and
            # beta arm.
            theta_rad = np.radians(theta_phi[:, 0])
            joints = np.array([xcen + (self.robot._ralpha + 0.25) * np.cos(theta_rad),
                               ycen + (self.robot._ralpha + 0.25) * np.sin(theta_rad)]).T

            # Plots the positioner axes, joints, and targets.
            ax.scatter(xcen, ycen, s=1, color='k')
            ax.scatter(joints[:, 0], joints[:, 1], s=0.1, color='k')
            ax.scatter(self.targets[:, 0], self.targets[:, 1], s=1, color='b')

            # Now we create patches for each arm on each positioner.
            n_positioner = 0
            for theta, phi in theta_phi:

                is_collisioned = collisions[n_positioner]

                theta_rad = np.radians(theta)

                alpha_arm_end = (xcen[n_positioner] + self.robot._ralpha * np.cos(theta_rad),
                                 ycen[n_positioner] + self.robot._ralpha * np.sin(theta_rad))

                ec = 'k' if not is_collisioned else 'r'

                # We use FancyBboxPatch to make the rectangles rounded.
                # The dimensions of the boxes are slightly modified to look
                # good on the plot.
                alpha_arm = matplotlib.patches.FancyBboxPatch(
                    (xcen[n_positioner] - 1, ycen[n_positioner] - arm_width / 2.),
                    self.robot._ralpha + 2, arm_width, boxstyle='round,pad=0.3,rounding_size=2',
                    fc='none', ec=ec, lw=0.2)

                # We need to add the patch before modify it so that ax.transData works.
                ax.add_patch(alpha_arm)

                alpha_transform = matplotlib.transforms.Affine2D().rotate_deg_around(
                    xcen[n_positioner], ycen[n_positioner], theta)
                alpha_arm.set_transform(alpha_transform + ax.transData)

                beta_arm = matplotlib.patches.FancyBboxPatch(
                    (alpha_arm_end[0] - 1, alpha_arm_end[1] - arm_width / 2.),
                    self.robot._rbeta + 2, arm_width, boxstyle='round,pad=0.3,rounding_size=2',
                    fc='none', ec=ec, lw=0.2)

                ax.add_patch(beta_arm)

                beta_transform = matplotlib.transforms.Affine2D().rotate_deg_around(
                    alpha_arm_end[0], alpha_arm_end[1], theta + phi)
                beta_arm.set_transform(beta_transform + ax.transData)

                n_positioner += 1

        ax.set_xlabel('xFocal (mm)')
        ax.set_xlabel('yFocal (mm)')

        return ax
