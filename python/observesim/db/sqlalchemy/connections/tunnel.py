#!/usr/bin/env python
# encoding: utf-8
#
# @Author: José Sánchez-Gallego
# @Date: Feb 8, 2018
# @Filename: tunnel.py
# @License: BSD 3-Clause
# @Copyright: José Sánchez-Gallego


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from observesim import log, config

from . import DatabaseConnection


profile = config['database']['tunnel']

user = profile['user']
dbname = profile['dbname']
host = profile['host']
port = profile['port']

database_connection_string = f'postgresql+psycopg2://{user}@{host}:{port}/{dbname}'

# Intial database connection creation and instances to be exported.

db = DatabaseConnection(database_connection_string=database_connection_string)
engine = db.engine
metadata = db.metadata
Session = db.Session
Base = db.Base


log.info(f'connected to database {dbname!r} using profile \'tunnel\'')
