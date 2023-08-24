import os

import numpy as np
import pandas as pd
from astropy.time import Time

from sdssdb.peewee.sdss5db import database
database.set_profile('operations')
from sdssdb.peewee.sdss5db import targetdb


def doneForObs(obs="APO"):
    os.environ["OBSERVATORY"] = obs.upper()
    from sdssdb.peewee.sdss5db import opsdb
    
    Field = targetdb.Field
    d2f = targetdb.DesignToField
    Cadence = targetdb.Cadence
    d2s = opsdb.DesignToStatus
    Status = opsdb.CompletionStatus

    complete = Status.get(label="done")

    query = d2s.select(Field.pk, Field.field_id,
                    Field.racen, Field.deccen,
                    d2s.mjd, Cadence.label)\
            .join(d2f, on=(d2f.design_id == d2s.design_id))\
            .join(Field).join(Cadence)\
            .switch(d2s)\
            .join(Status)\
            .where(Status.pk == complete.pk).dicts()

    field_mjds = pd.DataFrame(query)

    print(len(field_mjds), len(np.unique(field_mjds["pk"])))

    exp = opsdb.Exposure
    cfg = opsdb.Configuration

    query = exp.select(cfg.design_id, exp.start_time)\
            .join(cfg)\
            .join(d2s, on=(d2s.design_id == cfg.design_id))\
            .join(Status)\
            .where(Status.pk == complete.pk).dicts()

    design_ids = list()
    start_times = list()
    for q in query:
        design_ids.append(q["design"])
        start_times.append(q["start_time"])

    dt_times = Time(start_times)
    dt_times.format = "mjd"
    mjds = [t.value for t in dt_times]
