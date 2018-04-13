#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2018-04-10
# @Filename: utils.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# @Copyright: José Sánchez-Gallego


import numpy as np


__all__ = ['xy2tp', 'tp2xy']


def xy2tp(x, y, r_alpha=7.4, r_beta=15.0):
    """Converts ``(x, y)`` positions on the focal plane to ``(theta, phi)``.

    Because each position in the patrol area of the actuator can be reached
    with two different configurations this routine returns both solutions for
    each set of ``(x, y)`` positions. Note that when ``phi=0`` both solutions
    are identical.

    Parameters:
        x,y (`~numpy.array` or float):
            The x and y position of the target with respect to the axis of
            the positioner. Either a float or an array of floats.
        r_alpha,r_beta (float):
            The lengths of the actuator arms in millimetres.

    Returns:
        result:
            A 3-d array in which the first two dimensions are the
            ``(theta, phi)`` pairs for each input ``(x, y)``. The third
            dimension contains the two possible solutions for each ``(x, y)``
            pair. Unreachable positions return NaN values.

    """

    x = np.atleast_1d(x)
    y = np.atleast_1d(y)

    # We use the law of the cosines on the triangle formed by the two
    # positioner arms and the vector from the centre of the positioner to the
    # target (x, y) position.

    rr = np.hypot(x, y)
    rr[rr == 0] = np.nan

    cos_A = (r_alpha**2 + r_beta**2 - rr**2) / (2. * r_alpha * r_beta)
    cos_B = (rr**2 + r_alpha**2 - r_beta**2) / (2. * r_alpha * rr)
    cos_D = x / rr

    # Sets invalid angles to NaN
    cos_A[np.abs(cos_A) > 1] = np.nan
    cos_B[np.abs(cos_B) > 1] = np.nan
    cos_D[np.abs(cos_D) > 1] = np.nan

    D = np.arccos(cos_D)
    D[y < 0] = 2. * np.pi - D[y < 0]

    theta_1 = np.rad2deg(D - np.arccos(cos_B))
    phi_1 = np.rad2deg(np.pi - np.arccos(cos_A))

    theta_2 = np.rad2deg(D + np.arccos(cos_B))
    phi_2 = np.rad2deg(np.pi + np.arccos(cos_A))

    theta_phi = np.array([np.array([theta_1, phi_1]).T,
                          np.array([theta_2, phi_2]).T])

    valid = ~np.isnan(theta_phi)
    theta_phi[valid][theta_phi[valid] < 0] += 360.
    theta_phi[valid] = theta_phi[valid] % 360.

    return theta_phi


def tp2xy(theta, phi, r_alpha=7.4, r_beta=15.0):
    """Converts ``(theta, pi)`` robot positions to ``(x, y)``.

    Assumes that the axis of the robot positioner is at ``(0, 0)``.

    Parameters:
        theta,phi (`~numpy.array` or float):
            The theta and phi angles of the positioner arms.
            Either a float or an array of floats.
        r_alpha,r_beta (float):
            The lengths of the actuator arms in millimetres.

    Returns:
        result:
            A 2-d array of ``(x, y)`` pairs for each input ``(theta, phi)``.

    """

    x = r_alpha * np.cos(np.radians(theta)) + r_beta * np.cos(np.radians(phi + theta))
    y = r_alpha * np.sin(np.radians(theta)) + r_beta * np.sin(np.radians(phi + theta))

    return np.array([x, y]).T
