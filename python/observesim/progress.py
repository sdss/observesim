import os

import numpy as np
import pandas as pd
from astropy.time import Time

from sdssdb.peewee.sdss5db import database
database.set_profile('operations')
from sdssdb.peewee.sdss5db import targetdb


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
                      Version.plan == plan).dicts()

    field_mjds = pd.DataFrame(query)

    exp = opsdb.Exposure
    cfg = opsdb.Configuration
    cf = opsdb.CameraFrame
    camera = opsdb.Camera
    flavor = opsdb.ExposureFlavor

    science = flavor.get(label="Science")
    if obs.upper() == "APO":
        blue = camera.get(label="b1")
    else:
        blue = camera.get(label="b2")

    # for obs hist
    query = exp.select(cfg.design_id, exp.start_time)\
               .join(cfg)\
               .join(d2s, on=(d2s.design_id == cfg.design_id))\
               .join(Status)\
               .switch(exp).join(cf)\
               .where(Status.pk == complete.pk, 
                      exp.exposure_flavor_pk == science.pk,
                      cf.camera_pk == blue.pk).dicts()

    design_ids = list()
    start_times = list()
    for q in query:
        design_ids.append(q["design"])
        start_times.append(q["start_time"])

    dt_times = Time(start_times)
    dt_times.format = "mjd"
    mjds = [t.value for t in dt_times]

    return mjds, field_mjds
