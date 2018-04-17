#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2018-04-11
# @Filename: test_utils.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# @Copyright: José Sánchez-Gallego


import numpy as np
import pytest
from scipy.spatial.distance import cdist

from .. import utils


class TestThetaPhi(object):

    def test_xy2tp(self, robot):
        """Tests the utils.xy2tp function."""

        n_actuators = np.sum(robot.fiducial)

        thetas = np.random.sample(n_actuators) * 360.
        phis = np.random.sample(n_actuators) * 360.

        xy = utils.tp2xy(thetas, phis)

        theta_phi_recovered = utils.xy2tp(xy[:, 0], xy[:, 1])

        for ii in range(n_actuators):

            # Checks that one of the recovered solutions matches the original (theta, phi)
            solutions = theta_phi_recovered[:, ii, :]

            solution1_comp = (pytest.approx(solutions[0, 0]) == thetas[ii] and
                              pytest.approx(solutions[0, 1]) == phis[ii])
            solution2_comp = (pytest.approx(solutions[1, 0]) == thetas[ii] and
                              pytest.approx(solutions[1, 1]) == phis[ii])

            assert solution1_comp or solution2_comp

            # Checks that both solutions recover the (x, y) pair.
            xy_from_solution1 = utils.tp2xy(solutions[0, 0], solutions[0, 1])
            xy_from_solution2 = utils.tp2xy(solutions[1, 0], solutions[1, 1])

            assert pytest.approx(xy_from_solution1[0]) == xy[ii][0]
            assert pytest.approx(xy_from_solution1[1]) == xy[ii][1]

            assert pytest.approx(xy_from_solution2[0]) == xy[ii][0]
            assert pytest.approx(xy_from_solution2[1]) == xy[ii][1]

    def test_xy2tp_invalid(self):
        """Tests the utils.xy2tp function with invalid points."""

        # A point that is inside the inner ring
        assert np.all(np.isnan(utils.xy2tp(0, 0)))

        # A point that is not reachable
        assert np.all(np.isnan(utils.xy2tp(100, 100)))

    def test_xy2tp_phi_zero(self):
        """Tests the utils.xy2tp function when phi=0."""

        results = utils.xy2tp(22.4, 0., r_alpha=7.4, r_beta=15.0)

        assert pytest.approx(results[0, 0, 0]) == 0.
        assert pytest.approx(results[1, 0, 0]) == 0.

        assert pytest.approx(results[0, 0, 1], abs=1e-5) == 0.
        assert pytest.approx(results[1, 0, 1], abs=1e-5) == 360.


class TestTargetAllocation(object):

    def test_generate_mock_targets(self, robot):

        targets = utils.generate_mock_targets(robot, min_distance=3)

        assert len(targets) == (len(robot.xcen) - np.sum(robot.fiducial))

        distances = cdist(targets, targets, 'euclidean')
        distances_triu = distances[np.triu_indices(distances.shape[0], 1)]
        assert np.all(distances_triu > 3)

    def test_assign_targets_draining(self, robot):

        # Creates a long list of targets
        targets = np.vstack([utils.generate_mock_targets(robot, one_per_positioner=False)
                             for __ in range(5)])

        positioner_to_targets, target_to_positioners = utils.assign_targets_draining(
            robot, targets, return_target_to_positioners=True)

        assert len(positioner_to_targets) == len(robot.xcen) - np.sum(robot.fiducial)

        for positioners in target_to_positioners.values():
            assert len(positioners) > 0

        # Checks that each positioner has the shortest possible list of targets
        for positioner in positioner_to_targets:
            assigned_targets = positioner_to_targets[positioner]
            for target in assigned_targets:
                valid_positioners = target_to_positioners[target]
                assert positioner in valid_positioners
                for valid_positioner in valid_positioners:
                    n_target_valid_positioner = len(positioner_to_targets[valid_positioner])
                    assert len(assigned_targets) <= (n_target_valid_positioner + 1)
