#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2018-04-11
# @Filename: test_utils.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# @Copyright: José Sánchez-Gallego


import numpy as np
import pytest

from ..utils import tp2xy, xy2tp


def test_xy2tp(robot):
    """Tests the utils.xy2tp function."""

    n_actuators = np.sum(robot.fiducial)

    thetas = np.random.sample(n_actuators) * 360.
    phis = np.random.sample(n_actuators) * 360.

    xy = tp2xy(thetas, phis)

    theta_phi_recovered = xy2tp(xy[:, 0], xy[:, 1])

    for ii in range(n_actuators):

        # Checks that one of the recovered solutions matches the original (theta, phi)
        solutions = theta_phi_recovered[:, ii, :]

        solution1_comp = (pytest.approx(solutions[0, 0]) == thetas[ii] and
                          pytest.approx(solutions[0, 1]) == phis[ii])
        solution2_comp = (pytest.approx(solutions[1, 0]) == thetas[ii] and
                          pytest.approx(solutions[1, 1]) == phis[ii])

        assert solution1_comp or solution2_comp

        # Checks that both solutions recover the (x, y) pair.
        xy_from_solution1 = tp2xy(solutions[0, 0], solutions[0, 1])
        xy_from_solution2 = tp2xy(solutions[1, 0], solutions[1, 1])

        assert pytest.approx(xy_from_solution1[0]) == xy[ii][0]
        assert pytest.approx(xy_from_solution1[1]) == xy[ii][1]

        assert pytest.approx(xy_from_solution2[0]) == xy[ii][0]
        assert pytest.approx(xy_from_solution2[1]) == xy[ii][1]
        print(xy_from_solution1, xy[ii])
        print(xy_from_solution2, xy[ii])
