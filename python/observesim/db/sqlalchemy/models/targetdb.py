#!/usr/bin/env python
# encoding: utf-8
#
# @Author: José Sánchez-Gallego
# @Date: Feb 8, 2018
# @Filename: targetdb.py
# @License: BSD 3-Clause
# @Copyright: José Sánchez-Gallego


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from sqlalchemy.orm import configure_mappers, relation

from ..connections import DatabaseConnection


def __repr__(self):
    """A custom repr for targetdb models.

    By default it always prints pk, name, and label, if found. Models can
    define they own ``__print_fields__`` as a list of field to be output in the
    repr.

    """

    fields = ['pk={0!r}'.format(self.pk)]

    for ff in self.__print_fields__:
        if hasattr(self, ff):
            fields.append('{0}={1!r}'.format(ff, getattr(self, ff)))

    return '<{0}: {1}>'.format(self.__class__.__name__, ', '.join(fields))


db = DatabaseConnection()
Base = db.Base

Base.__print_fields__ = ['label', 'name']
Base.__repr__ = __repr__


class Target(Base):

    __tablename__ = 'target'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class File(Base):

    __tablename__ = 'file'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Field(Base):

    __tablename__ = 'field'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class TargetType(Base):

    __tablename__ = 'target_type'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class TargetCompletion(Base):

    __tablename__ = 'target_completion'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Magnitude(Base):

    __tablename__ = 'magnitude'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class StellarParams(Base):

    __tablename__ = 'stellar_params'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Program(Base):

    __tablename__ = 'program'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Spectrograph(Base):

    __tablename__ = 'spectrograph'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class TargetCadence(Base):

    __tablename__ = 'target_cadence'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Lunation(Base):

    __tablename__ = 'lunation'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Survey(Base):

    __tablename__ = 'survey'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Tile(Base):

    __tablename__ = 'tile'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class TargetToTile(Base):

    __tablename__ = 'target_to_tile'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Simulation(Base):

    __tablename__ = 'simulation'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Exposure(Base):

    __tablename__ = 'exposure'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Weather(Base):

    __tablename__ = 'weather'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Spectrum(Base):

    __tablename__ = 'spectrum'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class FiberConfiguration(Base):

    __tablename__ = 'fiber_configuration'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Fiber(Base):

    __tablename__ = 'fiber'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class FiberStatus(Base):

    __tablename__ = 'fiber_status'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class Actuator(Base):

    __tablename__ = 'actuator'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class ActuatorStatus(Base):

    __tablename__ = 'actuator_status'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


class ActuatorType(Base):

    __tablename__ = 'actuator_type'
    __table_args__ = {'autoload': True, 'schema': 'targetdb'}


Target.file = relation(File, backref='targets')
Target.field = relation(Field, backref='targets')
Target.target_type = relation(TargetType, backref='targets')
Target.target_completion = relation(TargetCompletion, backref='targets')
Target.magnitude = relation(Magnitude, backref='target')
Target.stellar_params = relation(StellarParams, backref='target')
Target.program = relation(Program, backref='targets')
Target.spectrograph = relation(Spectrograph, backref='targets')
Target.target_cadence = relation(TargetCadence, backref='targets')
Target.lunation = relation(Lunation, backref='targets')

Program.survey = relation(Survey, backref='programs')

Tile.targets = relation(Target, secondary=TargetToTile.__table__, backref='tiles')

Tile.simulation = relation(Simulation, backref='tiles')

Exposure.tile = relation(Tile, backref='exposures')
Exposure.weather = relation(Weather, backref='exposures')
Exposure.simulation = relation(Simulation, backref='exposures')

Spectrum.exposure = relation(Exposure, backref='spectra')
Spectrum.fiber_configuration = relation(FiberConfiguration, backref='spectrum')

FiberConfiguration.fiber = relation(Fiber, backref='fiber_configurations')
FiberConfiguration.tile = relation(Tile, backref='fiber_configurations')
FiberConfiguration.target = relation(Target, backref='fiber_configurations')

Fiber.fiber_status = relation(FiberStatus, backref='fibers')
Fiber.spectrograph = relation(Spectrograph, backref='fibers')
Fiber.actuator = relation(Actuator, backref='fibers')

Actuator.actuator_status = relation(ActuatorStatus, backref='actuators')
Actuator.actuator_type = relation(ActuatorType, backref='actuators')


try:
    configure_mappers()
except RuntimeError as error:
    raise RuntimeError('An error occurred when verifying the relationships '
                       'between the database tables. Most likely this is an error '
                       'in the definition of the SQLAlchemy relationships '
                       'see the error message below for details.\n\n' + str(error))
