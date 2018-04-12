#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2018-04-12
# @Filename: __init__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# @Copyright: José Sánchez-Gallego


import pathlib


def save_test_plot(fig, filename):
    """Saves a plot to test/plots."""

    plot_dir = pathlib.Path(__file__).parent / 'plots'

    if not plot_dir.exists():
        plot_dir.mkdir()

    fig.savefig(plot_dir / filename)
