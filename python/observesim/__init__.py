#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

"""
observesim
"""

import os
import yaml

# Inits the logging system. Only shell logging, and exception and warning catching.
# File logging can be started by calling log.start_file_logger(name).
# from .misc import log
# import observesim.scheduler
# import observesim.tile
# import observesim.observations
# import observesim.weather
# import observesim.observe


def merge(user, default):
    """Merges a user configuration with the default one."""

    if isinstance(user, dict) and isinstance(default, dict):
        for kk, vv in default.items():
            if kk not in user:
                user[kk] = vv
            else:
                user[kk] = merge(user[kk], vv)

    return user


NAME = 'observesim'

# # Loads config
# config_file = os.path.dirname(__file__) + '/etc/{0}.cfg'.format(NAME)
# config = yaml.load(open(config_file), Loader=yaml.FullLoader)

# # If there is a custom configuration file, updates the defaults using it.
# custom_config_fn = os.path.expanduser('~/.{0}/{0}.cfg'.format(NAME))
# if os.path.exists(custom_config_fn):
#     config = merge(yaml.load(open(custom_config_fn), Loader=yaml.FullLoader),
#                    config)


__version__ = '0.1.0'
