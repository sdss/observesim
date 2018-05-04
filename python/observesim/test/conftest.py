#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2018-04-12
# @Filename: conftest.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# @Copyright: José Sánchez-Gallego


import pathlib

import pytest

from ..robot import Robot


def pytest_addoption(parser):
    parser.addoption('--plot', action='store_true', default=False, help='outputs test plots')


def pytest_configure(config):
    """Runs during configuration of conftest."""

    do_plot = config.getoption('--plot')

    if do_plot:
        plots_path = (pathlib.Path(__file__).parent / 'plots')
        if plots_path.exists():
            for fn in plots_path.glob('*'):
                fn.unlink()


@pytest.fixture
def plot(request):
    return request.config.getoption('--plot')


@pytest.fixture
def robot():
    yield Robot()
