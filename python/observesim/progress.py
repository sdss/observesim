import os

import numpy as np
import pandas as pd
from astropy.time import Time

from sdssdb.peewee.sdss5db import database
database.set_profile('operations')
from sdssdb.peewee.sdss5db import targetdb


def fieldsFromDB(obs="APO", plan="eta-5"):
    Field = targetdb.Field
    Cadence = targetdb.Cadence
    Version = targetdb.Version
    Observatory = targetdb.Observatory

    query = Field.select(Field.pk, Field.field_id,
                    Field.racen, Field.deccen,
                    Field.slots_exposures,
                    Cadence.label.alias("cadence"))\
            .join(Cadence)\
            .switch(Field).join(Version)\
            .switch(Field).join(Observatory)\
            .where(Version.plan == plan,
                   Observatory.label == obs.upper())\
            .order_by(Field.pk).dicts()

    fieldTable = pd.DataFrame(query)

    fields_model = [('pk', np.int32),
                    ('field_id', np.int32),
                    ('racen', np.float64),
                    ('deccen', np.float64),
                    ('nfilled', np.int32),
                    ('flag', np.int32),
                    ('slots_exposures', np.int32, (24, 2)),
                    ('cadence', np.dtype('a20'))]

    fields_sum = np.zeros(len(fieldTable["field_id"]), dtype=fields_model)

    fields_sum["pk"] = fieldTable["pk"].values
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

    for f in fields_sum["pk"]:
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
    
    Field = targetdb.Field
    d2f = targetdb.DesignToField
    Cadence = targetdb.Cadence
    Version = targetdb.Version
    d2s = opsdb.DesignToStatus
    Status = opsdb.CompletionStatus

    complete = Status.get(label="done")

    query = d2s.select(Field.pk, Field.field_id,
                       Field.racen, Field.deccen,
                       d2s.mjd, Cadence.label, d2s.design_id)\
               .join(d2f, on=(d2f.design_id == d2s.design_id))\
               .join(Field).join(Cadence)\
               .switch(d2s)\
               .join(Status)\
               .switch(Field).join(Version)\
               .where(Status.pk == complete.pk,
                      Version.plan == plan)\
                .order_by(d2s.mjd).dicts()

    field_mjds = pd.DataFrame(query)

    return field_mjds




    # exp = opsdb.Exposure
    # cfg = opsdb.Configuration
    # cf = opsdb.CameraFrame
    # camera = opsdb.Camera
    # flavor = opsdb.ExposureFlavor

    # science = flavor.get(label="Science")
    # if obs.upper() == "APO":
    #     blue = camera.get(label="b1")
    # else:
    #     blue = camera.get(label="b2")

    # # for obs hist
    # query = exp.select(cfg.design_id, exp.start_time,
    #                    Field.pk.alias("field_pk"))\
    #            .join(cfg)\
    #            .join(d2s, on=(d2s.design_id == cfg.design_id))\
    #            .join(Status)\
    #            .switch(exp).join(cf)\
    #            .switch(d2s)\
    #            .join(d2f, on=(d2f.design_id == d2s.design_id))\
    #            .join(Field).join(Version)\
    #            .where(Status.pk == complete.pk, 
    #                   exp.exposure_flavor_pk == science.pk,
    #                   cf.camera_pk == blue.pk,
    #                   Version.plan == plan).dicts()

    # design_ids = list()
    # start_times = list()
    # for q in query:
    #     design_ids.append(q["design"])
    #     start_times.append(q["start_time"])

    # dt_times = Time(start_times)
    # dt_times.format = "mjd"
    # mjds = [t.value for t in dt_times]

    # return mjds, field_mjds

PK
FIELDID
RACEN
DECCEN
CADENCE
NOBSERVATIONS
OBSERVATIONS


# OBS
FIELD_PK
MJD
CADENCE
NFILLED
RACEN
DECCEN

NEXP_CUMUL
DESIGN_ID

AIRMASS
SKYBRIGHTNESS
LST
HA

DURATION
APGSN2
RSN2
BSN2