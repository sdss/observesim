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
from .misc import log
import observesim.scheduler
import observesim.fields
import observesim.observations
import observesim.weather
import observesim.observe

NAME = 'observesim'

# Loads config
config = yaml.load(open(os.path.dirname(__file__) + '/etc/{0}.cfg'.format(NAME)))


__version__ = '0.1.0'
