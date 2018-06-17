#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: field.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import os
import warnings

import numpy as np
import fitsio
import matplotlib.pyplot as plt

import observesim.cadence as cadence
import observesim.robot as robot

__all__ = ['Field']

"""Field module class.

Dependencies:

 numpy
 astropy
 matplotlib

"""


class Field(object):
    """Field class

    Parameters:
    ----------

    racen : np.float64
        boresight RA, J2000 deg

    deccen : np.float64
        boresight Dec, J2000 deg

    Attributes:
    ----------

    racen : np.float64
        boresight RA, J2000 deg

    deccen : np.float64
        boresight Dec, J2000 deg

    field_cadence : int, np.int32
        index of field cadence in cadencelist

    robot : Robot class
        instance of Robot

    target_x : ndarray of np.float64
        x positions of targets

    target_y : ndarray of np.float64
        y positions of targets

    target_cadence : ndarray of np.int32
        cadences of targets

    target_type : ndarray of strings
        target types ('boss' or 'apogee')

    Methods:
    -------

    read_cadences() : read cadences from a file
    read_targets() : read targets from a file
    assign() : assign targets to robots for cadence

    Notes:
    -----

    This class is definitely going to need to be refactored (please
    email me if you are reading this comment in 2025 ...).

"""
    def __init__(self, racen=None, deccen=None,
                 db=True, fps_layout='filled_hex'):
        self.robot = robot.Robot(db=db, fps_layout=fps_layout)
        self.racen = racen
        self.deccen = deccen
        self.cadencelist = cadence.CadenceList()
        self.field_cadence = 0
        return

    def _arrayify(self, quantity=None, dtype=np.float64):
        """Cast quantity as ndarray of numpy.float64"""
        try:
            length = len(quantity)
        except TypeError:
            length = 1
        return np.zeros(length, dtype=dtype) + quantity

    def radec2xy(self, ra=None, dec=None):
        # Yikes!
        scale = 218.
        x = (ra - self.racen) * np.cos(self.deccen * np.pi / 180.) * scale
        y = (dec - self.deccen) * scale
        return(x, y)

    def read_targets(self, filename=None):
        self.target_fits = fitsio.read(filename)
        self.ntarget = len(self.target_fits)
        self.target_ra = self.target_fits['ra']
        self.target_dec = self.target_fits['dec']
        self.target_x, self.target_y = self.radec2xy(self.target_ra,
                                                     self.target_dec)
        self.target_cadence = self.target_fits['cadence']
        self.target_type = np.array([t.decode().strip()
                                     for t in self.target_fits['type']])
        return

    def plot(self, epochs=None):
        if(epochs is None):
            epochs = np.arange(self.assignments.shape[1])
        else:
            epochs = self._arrayify(epochs, dtype=np.int32)
        colors = ['black', 'red', 'green', 'blue', 'cyan', 'purple']
        plt.scatter(self.robot.xcen, self.robot.ycen, s=3, color='black')
        plt.scatter(self.target_x, self.target_y, s=3, color='red')
        for irobot in np.arange(self.assignments.shape[0]):
            for iepoch in np.array(epochs):
                icolor = iepoch % len(colors)
                itarget = self.assignments[irobot, iepoch]
                if(itarget >= 0):
                    xst = self.robot.xcen[irobot]
                    yst = self.robot.ycen[irobot]
                    xnd = self.target_x[itarget]
                    ynd = self.target_y[itarget]
                    plt.plot([xst, xnd], [yst, ynd], color=colors[icolor])

    def assign(self):
        nexposures = self.cadencelist.cadences[self.field_cadence].nexposures
        self.assignments = (np.zeros((self.robot.npositioner, nexposures),
                                     dtype=np.int32) - 1)
        got_target = np.zeros(self.ntarget, dtype=np.int32)
        ok_cadence = np.zeros(self.target_cadence.max() + 1, np.bool)
        for icadence in np.unique(self.target_cadence):
            ok_cadence[icadence] = self.cadencelist.cadence_consistency(icadence, self.field_cadence)
        iok = np.where(ok_cadence[self.target_cadence])[0]
        if(len(iok) == 0):
            return
        for indx in np.arange(self.robot.npositioner):
            positionerid = self.robot.positionerid[indx]
            ileft = np.where(got_target[iok] == 0)[0]
            if(len(ileft) > 0):
                it = self.robot.targets(positionerid=positionerid,
                                        x=self.target_x[iok[ileft]],
                                        y=self.target_y[iok[ileft]],
                                        type=self.target_type[iok[ileft]])
                if(len(it) > 0):
                    itarget = self.cadencelist.pack_targets(self.target_cadence[iok[ileft[it]]], self.field_cadence)
                    got_target[iok[ileft[it[itarget]]]] = 1
                    self.assignments[indx, :] = iok[ileft[it[itarget]]]
        return
