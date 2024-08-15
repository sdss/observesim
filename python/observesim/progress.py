import os

import numpy as np
import pandas as pd
from peewee import JOIN

from sdssdb.peewee.sdss5db import database
database.set_profile('operations')
from sdssdb.peewee.sdss5db import targetdb


def fieldsFromDB(obs="APO", plan="eta-9"):
    os.environ["OBSERVATORY"] = obs.upper()
    from sdssdb.peewee.sdss5db import opsdb
    opsdb.database.connect()

    Field = targetdb.Field
    Cadence = targetdb.Cadence
    Version = targetdb.Version
    Observatory = targetdb.Observatory
    f2p = opsdb.FieldToPriority

    query = Field.select(Field.pk.alias("field_pk"), Field.field_id.alias("field_id"),
                    Field.racen, Field.deccen,
                    Field.slots_exposures,
                    Cadence.label.alias("cadence"))\
            .join(Cadence)\
            .switch(Field).join(Version)\
            .switch(Field).join(Observatory)\
            .switch(Field).join(f2p, JOIN.LEFT_OUTER)\
            .where(Version.plan == plan,
                    Observatory.label == obs.upper(),
                    f2p.field_priority_pk == None)\
            .order_by(Field.pk).dicts()

    fieldTable = pd.DataFrame(query)

    fields_model = [('pk', np.int32),
                    ('field_id', np.int32),
                    ('racen', np.float64),
                    ('deccen', np.float64),
                    ('nfilled', np.int32),
                    ('flag', np.int32),
                    ('slots_exposures', np.int32, (24, 2)),
                    ('cadence', np.dtype('a25'))]

    fields_sum = np.zeros(len(fieldTable["field_id"]), dtype=fields_model)

    fields_sum["pk"] = fieldTable["field_pk"].values
    fields_sum["field_id"] = fieldTable["field_id"].values
    fields_sum["racen"] = fieldTable["racen"].values
    fields_sum["deccen"] = fieldTable["deccen"].values
    fields_sum["nfilled"] = [0 for i in fields_sum["pk"]]
    fields_sum["slots_exposures"] = [i for i in fieldTable["slots_exposures"].values]
    fields_sum["cadence"] = fieldTable["cadence"].values

    Design = targetdb.Design
    d2f = targetdb.DesignToField

    dquery = Design.select(Design.design_id, d2f.exposure,
                           d2f.field_pk)\
                   .join(d2f).join(Field).join(Version)\
                   .switch(Field).join(Observatory)\
                   .where(Version.plan == plan,
                          Observatory.label == obs.upper()).dicts()
    
    fieldTable = pd.DataFrame(dquery)

    arrayThatIndexes = fieldTable.to_numpy()

    design_ids = arrayThatIndexes[:, 0]
    exposures = arrayThatIndexes[:, 1]
    field_pks = arrayThatIndexes[:, 2]

    all_designs = list()

    for i, f in enumerate(fields_sum["pk"]):
        f_designs = np.where(field_pks == f)
        exps = exposures[f_designs]
        args_in_order = f_designs[0][np.argsort(exps)]
        designs = design_ids[args_in_order]
        fields_sum["nfilled"][i] = len(designs)
        all_designs.append(list(designs))

    return fields_sum, all_designs


def doneForObs(obs="APO", plan="zeta-3"):
    os.environ["OBSERVATORY"] = obs.upper()
    from sdssdb.peewee.sdss5db import opsdb
    opsdb.database.connect()
    
    Field = targetdb.Field
    d2f = targetdb.DesignToField
    Cadence = targetdb.Cadence
    Version = targetdb.Version
    d2s = opsdb.DesignToStatus
    Status = opsdb.CompletionStatus
    f2p = opsdb.FieldToPriority

    complete = Status.get(label="done")

    query = d2s.select(Field.pk.alias("field_pk"),
                       Field.field_id,
                       Field.racen, Field.deccen,
                       d2s.mjd, Cadence.label, d2s.design_id)\
               .join(d2f, on=(d2f.design_id == d2s.design_id))\
               .join(Field).join(Cadence)\
               .switch(d2s)\
               .join(Status)\
               .switch(Field).join(Version)\
               .switch(Field).join(f2p, JOIN.LEFT_OUTER)\
               .where(Status.pk == complete.pk,
                      Version.plan == plan,
                      f2p.field_priority_pk == None)\
                .order_by(d2s.mjd).dicts()

    # field_mjds = pd.DataFrame(query)

    return query
