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

    n_epochs = catalogue[row['n_epochs_column']]

    if row['n_exps_column'] == '-':
        catalogue.add_column(table.Column([1] * len(catalogue), 'n_exps'))
        row['n_exps_column'] = 'n_exps'

    n_exps = catalogue[row['n_exps_column']]

    # We assume that either cadence or cadence_code exist
    if row['cadence_column'] != '-':

        cadences = catalogue[row['cadence_column']]

        target_cadence_data = np.array([cadences, n_epochs, n_exps]).T

        # Iterates on each unique combination of cadence parameters.
        for cadence, n_epoch, n_exp in np.unique(target_cadence_data, axis=0):
            target_cadence_dbo, __ = targetdb.TargetCadence.get_or_create(cadence=cadence,
                                                                          n_epochs=n_epoch,
                                                                          n_exp_per_epoch=n_exp)

            # Determines the indices in the catalogue that correspond to this combination
            # of cadence parameters.
            target_where = np.where((catalogue[row['cadence_column']] == cadence) &
                                    (catalogue[row['n_exps_column']] == n_exp) &
                                    (catalogue[row['n_epochs_column']] == n_epoch))

            target_cadence_pks[target_where] = target_cadence_dbo.pk

    elif row['cadence_code_column'] != '':

        cadence_codes = catalogue[row['cadence_code_column']]

        target_cadence_data = np.array([cadence_codes, n_epochs, n_exps]).T
        for cadence_code, n_epoch, n_exp in np.unique(target_cadence_data, axis=0):
            target_cadence_dbo, __ = targetdb.TargetCadence.get_or_create(
                cadence_code=cadence_code,
                n_epochs=n_epoch,
                n_exp_per_epoch=n_exp)

            target_where = np.where((catalogue[row['cadence_code_column']] == cadence_code) &
                                    (catalogue[row['n_exps_column']] == n_exp) &
                                    (catalogue[row['n_epochs_column']] == n_epoch))

            target_cadence_pks[target_where] = target_cadence_dbo.pk

    else:

        raise ValueError('cannot find cadence or cadence_code column')

    return target_cadence_pks


def load_targetdb(filename, verbose=False):
    """Populates targetdb.

    Parameters:
        filename (str):
            The path to the table describing the list of files and columns to
            insert.
        verbose (bool):
            Determines the level of verbosity.

    """

    if verbose:
        log.sh.setLevel(logging.DEBUG)

    filename = pathlib.Path(filename)
    assert filename.exists()

    target_files = table.Table.read(filename, format='ascii')

    log.info('Deleting all target records ...')
    targetdb.Target.delete().execute()

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
        ra = catalogue[row['ra_column']].tolist()
        dec = catalogue[row['dec_column']].tolist()
        spectrograph_pks = [spectrograph_pk] * len(ra)
        program_pks = [program_pk] * len(ra)
        file_pks = [file_pk] * len(ra)
        file_index = list(range(len(ra)))
        target_completion_pks = [0] * len(ra)  # Sets it to automatic

        # If the field column is set, determines the corresponding field
        # for each target.
        if row['field_column'] == '-':
            field_pks = [None] * len(ra)
        else:
            field_pks = np.zeros(len(ra))

            for field in np.unique(catalogue[row['field_column']]):
                field_dbo, __ = targetdb.Field.get_or_create(label=field.strip())
                field_pk = field_dbo.pk
                field_pks[np.where(catalogue[row['field_column']] == field)] = field_pk

        # Creates unique cadence rows and returns a list of pks to assign to target_cadence_pk
        log.debug('adding cadence data')
        target_cadence_pks = _create_cadence_records(catalogue, row)

        # Bulk loads all the data
        rows = zip(ra, dec, spectrograph_pks, program_pks, file_pks,
                   file_index, target_completion_pks, field_pks, target_cadence_pks)

        targetdb.Target.insert_many(rows, fields=[targetdb.Target.ra,
                                                  targetdb.Target.dec,
                                                  targetdb.Target.spectrograph_pk,
                                                  targetdb.Target.program_pk,
                                                  targetdb.Target.file_pk,
                                                  targetdb.Target.file_index,
                                                  targetdb.Target.target_completion_pk,
                                                  targetdb.Target.field_pk,
                                                  targetdb.Target.target_cadence_pk]).execute()
