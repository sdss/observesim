#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2018-04-10
# @Filename: utils.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# @Copyright: José Sánchez-Gallego


import numpy as np


__all__ = ['xy2tp', 'tp2xy', 'generate_mock_targets',
           'assign_targets_draining']


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


def generate_mock_targets(robot, min_distance=3, one_per_positioner=True):
    """An function to generate fully allocatable sets of targets.

    Generates a list of ``(x, y)`` positions on the focal plane ensuring that
    every actuator has a target within reach and that the targets are separated
    a ``min_distance``.

    Parameters:
        robot (`~observesim.robot.Robot` object):
            The `~observesim.robot.Robot` layout and geometry.
        min_distance (float):
            The minimum distance, in millimetres, that two targets can be
            separated.
        one_per_positioner (bool):
            If ``True``, ensures that each positioner has an allocated target
            (although a target may be in reach of multiple positioners).
            Otherwise, targets are distributed randomly but always at a
            position on the focal plane where at least a positioner can reach
            them.

    Returns:
        targets (`~numpy.ndarray`):
            A Numpy array with the ``(x, y)`` position on the focal plane. The
            ordering is such that it matches the order of the actuators in
            ``robot``.

    """

    positioners = np.array([robot.xcen[~robot.fiducial],
                            robot.ycen[~robot.fiducial]]).T
    n_positioners = positioners.shape[0]

    targets = np.zeros(positioners.shape)

    for nn in range(n_positioners):

        valid = False

        while not valid:

            # Determines whether the new target will be placed around a
            # positioner in sequential or random order.
            if one_per_positioner:
                positioner = positioners[nn, :]
            else:
                positioner = positioners[np.random.randint(n_positioners)]

            new_target_tp = np.random.uniform(0, 360, size=2)
            new_target_xy = tp2xy(new_target_tp[0], new_target_tp[1],
                                  r_alpha=robot._ralpha, r_beta=robot._rbeta) + positioner

            current_targets = targets[:nn]

            if len(current_targets) > 0:
                distances = np.sqrt(np.sum((current_targets - new_target_xy)**2, axis=1))
                if np.any(distances <= min_distance):
                    continue

            valid = True
            targets[nn, :] = new_target_xy

    return targets


def assign_targets_draining(robot, targets, return_target_to_positioners=False):
    """Implements the target-fibre allocation method by Morales et al. (2012)

    Uses a draining algorithm to optimally assign targets to
    fibres/positioners. The method implements the following steps:

    i) For each target, creates a list of available positioners.
    ii) Creates a list of positioners. Assigns each target to the valid
        positioner that has a shorter list of targets.
    iii) Moves targets between positioners (as permitted by the lists defined
         in ``(i)`` until all list are as short as possible.

    Parameters:
        robot (`~observesim.robot.Robot` object):
            The `~observesim.robot.Robot` layout and geometry.
        targets (`numpy.ndarray`):
            An ``Nx2`` array with the ``(x, y)`` positions of the ``N`` targets
            on the focal plane (in mm).
        return_target_to_positioners (bool):
            If ``True``, also returns a dictionary of indices from with valid
            positioners for each target in ``targets``.

    Returns:
        indices (dict):
            A dictionary in which the keys are the indices of the actuators in
            ``robot``. For each element, the value is a list on the ``target``
            indices associated to that actuator.

    Example:
        ::

            >> from observesim.robot import Robot
            >> from observesim.utils import assign_targets_draining, generate_mock_targets
            >> robot = Robot()
            >> targets = generate_mock_targets(robot)
            >> assignment = assign_targets_draining(robot, targets)
            >> print(assignment)
            {0: [72], 1: [15], ...}
            >> actuator_0_xcen = robot.xcen[0]
            >> targets_for_actuator_0 = [targets[ii, :] for ii in assignment[0]]

    """

    targets = np.atleast_2d(targets)
    assert targets.shape[1] == 2, 'invalid dimensions for target array.'

    pos_idx = np.where(~robot.fiducial)[0]
    pos_xy = np.array([robot.xcen[pos_idx], robot.ycen[pos_idx]]).T

    target_to_positioners = {}
    positioner_to_targets = {ii: [] for ii in pos_idx}
    for target_idx in range(targets.shape[0]):

        # Determines which actuators have a (theta, phi) configuration that
        # would allow them to reach the target
        target_xy = targets[target_idx, :]
        target_distance = target_xy - pos_xy
        theta_phi = xy2tp(target_distance[:, 0], target_distance[:, 1],
                          r_alpha=robot._ralpha, r_beta=robot._rbeta)

        # Valid positioners are those that have (theta, phi) not NaN
        valid_positioners = np.unique(np.where(~np.isnan(theta_phi))[1])

        # valid_positioners are indices of pos_xy. We need to convert them
        # to pos_idx values.
        valid_positioner_idx = [pos_idx[ii] for ii in valid_positioners]

        # Assigns the list of valid positioners for this target
        target_to_positioners[target_idx] = valid_positioner_idx

        # Selects the valid positioner with fewer targets and assigns this
        # target to it.
        n_target_in_valid_positioners = [len(positioner_to_targets[ii])
                                         for ii in valid_positioner_idx]

        min_positioner = valid_positioner_idx[np.argmin(n_target_in_valid_positioners)]
        positioner_to_targets[min_positioner].append(target_idx)

    # Implements step (iii), reassigning targets to positioners until all
    # lists are optimally short.

    changed = True
    while changed is True:

        changed = False
        new_positioner_to_targets = positioner_to_targets.copy()

        for positioner in positioner_to_targets:
            for target in positioner_to_targets[positioner]:
                valid_positioners = target_to_positioners[target]
                for new_positioner in valid_positioners:
                    origin_length = len(new_positioner_to_targets[positioner])
                    dest_length = len(new_positioner_to_targets[new_positioner])
                    if dest_length + 1 < origin_length:
                        new_positioner_to_targets[new_positioner] += [target]
                        new_positioner_to_targets[positioner].remove(target)
                        changed = True
                        break

    if return_target_to_positioners is False:
        return positioner_to_targets
    else:
        return positioner_to_targets, target_to_positioners
