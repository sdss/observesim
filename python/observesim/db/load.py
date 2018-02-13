#!/usr/bin/env python
# encoding: utf-8
#
# @Author: José Sánchez-Gallego
# @Date: Feb 6, 2018
# @Filename: load.py
# @License: BSD 3-Clause
# @Copyright: José Sánchez-Gallego


import logging
import pathlib

import numpy as np

from astropy import table

from observesim import log
from .peewee import targetdb


__all__ = ['load_targets', 'load_fibres']


def _create_program_record(program, survey):
    """Creates or retrieves a ``program`` for a given ``survey``.

    Return the ``pk`` of the program.

    """

    program_dbo, created = targetdb.Program.get_or_create(label=program.strip())

    if created:
        log.debug(f'created record for program {program!r}')
    else:
        log.debug(f'found record pk={program_dbo.pk} for program {program!r}')

    if program_dbo.survey is not None:
        assert program_dbo.survey.label == survey, \
            'the program survey in the DB and the input value do not match.'
    else:
        with targetdb.database.atomic():
            program_dbo.survey_pk = targetdb.Survey.get(label=survey.strip()).pk
            program_dbo.save()
            log.debug(f'associated program {program!r} to survey {survey!r}')

    return program_dbo.pk


def _create_cadence_records(catalogue, row):
    """Creates cadence records and returns a list of pks to associate to ``catalogue``."""

    target_cadence_pks = np.zeros(len(catalogue))

    n_epochs = catalogue[row['n_epochs']]

    if row['n_exps'] == '-':
        catalogue.add_column(table.Column([1] * len(catalogue), 'n_exps'))
        row['n_exps'] = 'n_exps'

    n_exps = catalogue[row['n_exps']]

    # We assume that either cadence or cadence_code exist
    if row['cadence'] != '-':

        cadences = catalogue[row['cadence']]

        target_cadence_data = np.array([cadences, n_epochs, n_exps]).T

        # Iterates on each unique combination of cadence parameters.
        for cadence, n_epoch, n_exp in np.unique(target_cadence_data, axis=0):
            target_cadence_dbo, __ = targetdb.TargetCadence.get_or_create(cadence=cadence,
                                                                          n_epochs=n_epoch,
                                                                          n_exp_per_epoch=n_exp)

            # Determines the indices in the catalogue that correspond to this combination
            # of cadence parameters.
            target_where = np.where((catalogue[row['cadence']] == cadence) &
                                    (catalogue[row['n_exps']] == n_exp) &
                                    (catalogue[row['n_epochs']] == n_epoch))

            target_cadence_pks[target_where] = target_cadence_dbo.pk

    elif row['cadence_code'] != '':

        cadence_codes = catalogue[row['cadence_code']]

        target_cadence_data = np.array([cadence_codes, n_epochs, n_exps]).T
        for cadence_code, n_epoch, n_exp in np.unique(target_cadence_data, axis=0):
            target_cadence_dbo, __ = targetdb.TargetCadence.get_or_create(
                cadence_code=cadence_code,
                n_epochs=n_epoch,
                n_exp_per_epoch=n_exp)

            target_where = np.where((catalogue[row['cadence_code']] == cadence_code) &
                                    (catalogue[row['n_exps']] == n_exp) &
                                    (catalogue[row['n_epochs']] == n_epoch))

            target_cadence_pks[target_where] = target_cadence_dbo.pk

    else:

        raise ValueError('cannot find cadence or cadence_code column')

    return target_cadence_pks


def _load_table_from_cols(columns, catalogue, row, db_model):
    """Loads a set of columns into a model."""

    catalogue_columns = [row[col] for col in columns]
    db_columns = [getattr(db_model, col) for col in columns]

    data = zip(*[catalogue[col].tolist() for col in catalogue_columns])

    pks = db_model.insert_many(data, fields=db_columns).execute()

    return list(zip(*list(pks)))[0]


def _load_magnitude(catalogue, row):
    """Loads the magnitude table and returns a list with the inserted pks."""

    valid_columns = [col for col in row.colnames if '_mag' in col and row[col] != '-']
    if len(valid_columns) == 0:
        return None

    return _load_table_from_cols(valid_columns, catalogue, row, targetdb.Magnitude)


def _load_stellar_params(catalogue, row):
    """Loads the stellar params table and returns a list of pks."""

    possible_columns = ['distance', 'teff', 'logg', 'mass', 'spectral_type', 'age']

    valid_columns = [col for col in possible_columns if row[col] != '-']
    if len(valid_columns) == 0:
        return None

    return _load_table_from_cols(valid_columns, catalogue, row, targetdb.StellarParams)


def load_targets(filename, verbose=False, remove=False):
    """Populates ``targetdb.target`` and associated tables.

    Parameters:
        filename (str):
            The path to the table describing the list of files and columns to
            insert.
        verbose (bool):
            Determines the level of verbosity.
        remove (bool):
            If ``True``, the target, magnitude, field, program, target cadence,
            file, and stellar paramater tables will be truncated before
            inserting new records.

    """

    if verbose:
        log.sh.setLevel(logging.DEBUG)
    else:
        log.sh.setLevel(logging.INFO)

    filename = pathlib.Path(filename)
    assert filename.exists()

    target_files = table.Table.read(filename, format='ascii')

    if remove:
        models = [targetdb.Target, targetdb.TargetCadence,
                  targetdb.Field, targetdb.File, targetdb.Program,
                  targetdb.Magnitude, targetdb.StellarParams]

        for model in models:
            log.info('Deleting all records in table {} ...'.format(model._meta.table_name))
            model.delete().execute()

    for row in target_files:

        fn = pathlib.Path(row['filename'])
        log.info(f'processing file {fn}')
        assert fn.exists(), f'cannot find catalogue {fn!s}'

        # Reads the catalogue
        catalogue = table.Table.read(fn)

        # Determines the spectrograph
        spectrograph = row['instrument']
        log.debug(f'determining spectrograph pk for {spectrograph}')
        spectrograph_dbo = targetdb.Spectrograph.get_or_none(label=spectrograph)
        assert spectrograph_dbo is not None, 'cannot find spectrograph record.'
        spectrograph_pk = spectrograph_dbo.pk
        log.debug(f'spectrograph pk={spectrograph_pk}')

        # Inserts or recovers the program pk
        program_pk = _create_program_record(row['program'], row['survey'])

        # Inserts or recovers the file pk
        file_dbo, __ = targetdb.File.get_or_create(filename=str(fn.name))
        file_pk = file_dbo.pk

        # Creates lists for each one of the columns to insert
        ra = catalogue[row['ra']].tolist()
        dec = catalogue[row['dec']].tolist()
        spectrograph_pks = [spectrograph_pk] * len(catalogue)
        program_pks = [program_pk] * len(catalogue)
        file_pks = [file_pk] * len(catalogue)
        file_index = list(range(len(catalogue)))
        target_completion_pks = [0] * len(catalogue)  # Sets it to automatic

        # If the field column is set, determines the corresponding field
        # for each target.
        if row['field'] == '-':
            field_pks = [None] * len(catalogue)
        else:
            field_pks = np.zeros(len(catalogue))

            for field in np.unique(catalogue[row['field']]):
                field_dbo, __ = targetdb.Field.get_or_create(label=field.strip())
                field_pk = field_dbo.pk
                field_pks[np.where(catalogue[row['field']] == field)] = field_pk

        # Creates unique cadence rows and returns a list of pks to assign to target_cadence_pk
        log.debug('adding cadence data')
        target_cadence_pks = _create_cadence_records(catalogue, row)

        # Loads magnitude data
        magnitude_pks = _load_magnitude(catalogue, row)
        if magnitude_pks is None:
            magnitude_pks = [None] * len(catalogue)
        else:
            log.debug('magnitude data added')

        # Loads stellar parameters
        stellar_params_pks = _load_stellar_params(catalogue, row)
        if stellar_params_pks is None:
            stellar_params_pks = [None] * len(catalogue)
        else:
            log.debug('stellar parameters data added')

        # Gets the pk for science target type and assigns it to all targets
        science_target_pk = targetdb.TargetType.get(label='Science').pk
        target_type_pks = [science_target_pk] * len(catalogue)

        # Bulk loads all the data
        rows = zip(ra, dec, spectrograph_pks, program_pks, file_pks,
                   file_index, target_completion_pks, field_pks, target_cadence_pks,
                   magnitude_pks, stellar_params_pks, target_type_pks)

        targetdb.Target.insert_many(rows, fields=[targetdb.Target.ra,
                                                  targetdb.Target.dec,
                                                  targetdb.Target.spectrograph_pk,
                                                  targetdb.Target.program_pk,
                                                  targetdb.Target.file_pk,
                                                  targetdb.Target.file_index,
                                                  targetdb.Target.target_completion_pk,
                                                  targetdb.Target.field_pk,
                                                  targetdb.Target.target_cadence_pk,
                                                  targetdb.Target.magnitude_pk,
                                                  targetdb.Target.stellar_params_pk,
                                                  targetdb.Target.target_type_pk]).execute()


def load_fibres(filename, verbose=False, remove=False):
    """Loads fibres and actuators into ``targetdb``.

    Parameters:
        filename (str):
            Path to the file containing the positions of the actuators and
            the fibres associated to each one.
        verbose (bool):
            Determines the level of verbosity.
        remove (bool):
            If ``True``, the fiber, actuator, and fiber configuation tables
            will be truncated before inserting new records.

    """

    if verbose:
        log.sh.setLevel(logging.DEBUG)
    else:
        log.sh.setLevel(logging.INFO)

    filename = pathlib.Path(filename)
    assert filename.exists()

    positions = table.Table.read(filename, format='ascii.commented_header',
                                 names=['row', 'pos', 'x', 'y', 'assignment'])

    if remove:
        log.info('Deleting all actuator and fibre records ...')
        targetdb.Actuator.delete().execute()
        targetdb.Fiber.delete().execute()
        targetdb.FiberConfiguration.delete().execute()

    boss_spec_pk = targetdb.Spectrograph.get(label='BOSS').pk
    apogee_spec_pk = targetdb.Spectrograph.get(label='APOGEE').pk

    fibre_status_ok = targetdb.FiberStatus.get(label='OK').pk
    actuator_status_ok = targetdb.ActuatorStatus.get(label='OK').pk

    log.info('loading actuator and fibre data ... ')

    fiberid = 1

    with targetdb.database.atomic():

        for ii, position in enumerate(positions):

            if position['assignment'] == 'Fiducial':
                actuator_type = 'Fiducial'
                specs = []
            else:
                actuator_type = 'Robot'
                if position['assignment'] == 'BA':
                    specs = ['BOSS', 'APOGEE']
                elif position['assignment'] == 'BOSS':
                    specs = ['BOSS']
                else:
                    raise ValueError('invalid position of type {}'.format(position['assignment']))

            xcen = position['x']
            ycen = position['y']

            actuator = targetdb.Actuator.create(
                id=(ii + 1), xcen=xcen, ycen=ycen,
                actuator_status_pk=actuator_status_ok,
                actuator_type_pk=targetdb.ActuatorType.get(label=actuator_type).pk)

            if verbose:
                log.debug(f'created actuator id={actuator.id} of type '
                          f'{actuator_type!r} at ({xcen}, {ycen})')

            for spec in specs:
                fibre = targetdb.Fiber.create(
                    fiberid=fiberid, throughput=1,
                    spectrograph_pk=(boss_spec_pk if spec == 'BOSS' else apogee_spec_pk),
                    fiber_status_pk=fibre_status_ok, actuator_pk=actuator.pk)

                if verbose:
                    log.debug(f'created fibre id={fibre.fiberid} associated with '
                              f'spectrograph {spec} and actuator id={actuator.id}')
