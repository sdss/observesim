import os
import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import fitsio
from PyAstronomy.pyasl.asl.astroTimeLegacy import daycnv
import astropy.coordinates as coord
import astropy.units as u

__all__ = ["cartonCoverageByYearMaps", "cartonCoverageCumulative",
           "runAllCumulativeEpochs", "cumulativeWebpage"]

def sortEpochByCarton(carton, obs_data, rsTarg, cad2exp):
    obs_carton = obs_data[np.where(obs_data["carton"] == carton)]
    targs = rsTarg[np.where(rsTarg["carton"] == carton)]

    print(f"working on {carton} with {len(targs)} targs and {len(obs_carton)} obs")

    cids = np.unique(obs_carton["catalogid"])
    carton_epoch_mjds = list()
    done_ra = list()
    done_dec = list()
    field_pk = list()
    cid = list()
    for ci in cids:
        # loop through ids, there is an id per target-exposure in obs_data
        # assumes cadences are respected and grabs the mjd each epoch is completed
        targ = targs[np.where(targs["catalogid"] == ci)]
        obs = obs_carton[np.where(obs_carton["catalogid"] == ci)]
        nexp = cad2exp[str(targ["cadence"][0])]

        # this magic bit of indexing grabs every nexp_th element
        # syntax is [starting_index::step_size], don't forget it's 0-indexed
        epoch_mjds = obs["obs_mjd"][nexp-1::nexp]
        ra = obs["ra"][nexp-1::nexp]
        dec = obs["dec"][nexp-1::nexp]
        pks = obs["field_pk"][nexp-1::nexp]
        catalog_id = obs["catalogid"][nexp-1::nexp]
        carton_epoch_mjds.extend(epoch_mjds)
        done_ra.extend(ra)
        done_dec.extend(dec)
        field_pk.extend(pks)
        cid.extend(catalog_id)

    return carton_epoch_mjds, done_ra, done_dec, field_pk, cid


def cartonCoverageByYearMaps(base, plan, rs_base, version=None):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    res_base = v_base + "byYearCartons"

    try:
        os.makedirs(res_base)
    except:
        pass

    obs_data_apo = fitsio.read(v_base + f"obsTargets-{plan}-apo-0.fits",
                               columns=["obs_mjd", "catalogid", "ra", "dec",
                                        "carton", "field_pk"])
    obs_data_lco = fitsio.read(v_base + f"obsTargets-{plan}-lco-0.fits",
                               columns=["obs_mjd", "catalogid", "ra", "dec",
                                        "carton", "field_pk"])

    rsTargApo = fitsio.read(rs_base +
                            f"/{plan}/final/rsTargetsFinal-{plan}-apo.fits",
                            columns=["catalogid", "cadence", "carton"])

    rsCadenceApo = fitsio.read(rs_base +
                            f"/{plan}/final/rsCadencesFinal-{plan}-apo.fits")

    rsTargLco = fitsio.read(rs_base +
                            f"/{plan}/final/rsTargetsFinal-{plan}-lco.fits",
                            columns=["catalogid", "cadence", "carton"])

    rsCadenceLco = fitsio.read(rs_base +
                            f"/{plan}/final/rsCadencesFinal-{plan}-lco.fits")

    cartons = np.union1d(np.unique(obs_data_apo["carton"]),
                         np.unique(obs_data_lco["carton"]))
    
    cad2exp_apo = {c["CADENCE"].split("_v")[0]: c["NEXP"][0] for c in rsCadenceApo}
    cad2exp_lco = {c["CADENCE"].split("_v")[0]: c["NEXP"][0] for c in rsCadenceLco}

    mk = "H"
    mks = 15
    alpha = 0.4
    color = "dodgerblue"

    dtype = [('mjd', np.float64),
             ('ra', np.float64),
             ('dec', np.float64),
             ('field_pk', np.int32),
             ('catalog_id', np.int64)]

    for c in cartons:
        if "ops" in c:
            continue

        apo_file_name = os.path.join(res_base, f"{c}_mjds_apo.fits")
        lco_file_name = os.path.join(res_base, f"{c}_mjds_lco.fits")

        if os.path.isfile(apo_file_name):
            try:
                apo = fitsio.read(apo_file_name)
            except OSError:
                apo = np.zeros(0, dtype=dtype)
        else:
            mjd, ra, dec, fpk, cid = sortEpochByCarton(c, obs_data_apo, 
                                                       rsTargApo,
                                                       cad2exp_apo)
            apo = np.zeros(len(mjd), dtype=dtype)
            apo["mjd"] = np.array(mjd)
            apo["ra"] = np.array(ra)
            apo["dec"] = np.array(dec)
            apo["field_pk"] = np.array(fpk)
            apo["catalog_id"] = np.array(cid)
            fitsio.write(apo_file_name, apo)
        if os.path.isfile(lco_file_name):
            try:
                lco = fitsio.read(lco_file_name)
            except OSError:
                lco = np.zeros(0, dtype=dtype)
        else:
            mjd, ra, dec, fpk, cid = sortEpochByCarton(c, obs_data_lco, 
                                                       rsTargLco,
                                                       cad2exp_lco)
            lco = np.zeros(len(mjd), dtype=dtype)
            lco["mjd"] = np.array(mjd)
            lco["ra"] = np.array(ra)
            lco["dec"] = np.array(dec)
            lco["field_pk"] = np.array(fpk)
            lco["catalog_id"] = np.array(cid)
            fitsio.write(lco_file_name, lco)

        jan_2022 = 59580
        cuts = jan_2022 + np.arange(365, 365.25*6, 365.25)

        apo_ra_raw = coord.Angle(-(np.array(apo["ra"]) + 90) * u.degree)
        apo_ra = apo_ra_raw.wrap_at(180 * u.degree)
        apo_dec = coord.Angle(np.array(apo["dec"]) * u.degree)

        lco_ra_raw = coord.Angle(-(np.array(lco["ra"]) + 90) * u.degree)
        lco_ra = lco_ra_raw.wrap_at(180 * u.degree)
        lco_dec = coord.Angle(np.array(lco["dec"]) * u.degree)

        for y, cut in enumerate(cuts):
            w, h = 16, 9
            f = plt.figure(figsize=(w, h))
            ax1 = plt.subplot(111, projection='mollweide')

            w_apo = np.where(apo["mjd"] <= cut)
            w_lco = np.where(lco["mjd"] <= cut)

            ax1.scatter(lco_ra[w_lco].radian, lco_dec[w_lco].radian, s=mks,
                        alpha=alpha, c=color, marker=mk)
            ax1.scatter(apo_ra[w_apo].radian, apo_dec[w_apo].radian, s=mks,
                        alpha=alpha, c=color, marker=mk)

            f.suptitle(f"{c}_year{y+1}", fontsize=15)

            plt.savefig(f"{res_base}/{c}_year{y+1}.png")

            plt.close()

    return None

def cartonCoverageCumulative(base, plan, version=None, cartons=[], name=None):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    res_base = v_base + "byYearCartons"

    carton_data = dict()

    start = 99999
    stop = 0

    for c in cartons:
        apo_file_name = os.path.join(res_base, f"{c}_mjds_apo.fits")
        lco_file_name = os.path.join(res_base, f"{c}_mjds_lco.fits")

        apo = fitsio.read(apo_file_name)
        lco = fitsio.read(lco_file_name)

        mjds = np.union1d(apo["mjd"], lco["mjd"])

        min_mjd = np.min(mjds)
        max_mjd = np.max(mjds)

        if min_mjd < start:
            start = int(min_mjd)
        if max_mjd > stop:
            stop = int(max_mjd)

        carton_data[f"{c}_apo"] = apo
        carton_data[f"{c}_lco"] = lco

    assert start < stop, "mjd acquisition failed"

    mjds = np.arange(start, stop, 10)

    done_cartons_apo = dict()
    done_cartons_lco = dict()
    for c in cartons:
        done_carton_apo = list()
        done_carton_lco = list()
        for m in mjds:
            done_apo = np.where(carton_data[f"{c}_apo"]["mjd"] <= m)
            done_lco = np.where(carton_data[f"{c}_lco"]["mjd"] <= m)
            done_carton_apo.append(len(done_apo[0]))
            done_carton_lco.append(len(done_lco[0]))
        done_cartons_apo[c] = done_carton_apo
        done_cartons_lco[c] = done_carton_lco
    

    plt.figure(figsize=(10,8))

    ax = plt.subplot(111)

    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    years_fmt = mdates.DateFormatter('%Y')

    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(years_fmt)
    ax.xaxis.set_minor_locator(months)

    cal_days = daycnv(mjds+2400000.5, mode="dt")

    for c in cartons:
        ax.plot(cal_days, done_cartons_apo[c], label=f"{c}_apo")
        ax.plot(cal_days, done_cartons_lco[c], label=f"{c}_lco")
    
    ax.legend(loc="upper left")
    ax.set_xlabel("Date", fontsize=16)
    ax.set_ylabel("N epochs", fontsize=16)

    if name is None:
        name = "_".join([c.split("_")[1] for c in cartons])


    plt.savefig(f"{res_base}/{name}_cumulative.png")

    plt.close()
    

def cartonEpochCumulative(base, plan, version=None, carton=None, loc="apo"):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    res_base = v_base + "byYearCartons"

    loc = loc.lower()
    file_name = os.path.join(res_base, f"{carton}_mjds_{loc}.fits")
    try:
        obs = fitsio.read(file_name)
    except OSError:
        return

    print("epoch cumulative", carton, loc)

    start = np.min(obs["mjd"])
    stop = np.max(obs["mjd"])
    mjds = np.arange(start, stop, 10)

    done_one = list()
    done_two = list()
    done_4plus = list()
    for m in mjds:
        done = obs[np.where(obs["mjd"] <= m)]
        unique, counts = np.unique(done["catalog_id"], return_counts=True)
        one = np.where(counts >= 1)
        two = np.where(counts >= 2)
        plus = np.where(counts > 2)
        done_one.append(len(one[0]))
        done_two.append(len(two[0]))
        done_4plus.append(len(plus[0]))

    plt.figure(figsize=(10,8))

    ax = plt.subplot(111)

    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    years_fmt = mdates.DateFormatter('%Y')

    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(years_fmt)
    ax.xaxis.set_minor_locator(months)

    cal_days = daycnv(mjds+2400000.5, mode="dt")

    ax.plot(cal_days, done_one, label="1 epoch")
    ax.plot(cal_days, done_two, label="2 epochs")
    ax.plot(cal_days, done_4plus, label="4+ epochs")
    ax.legend()
    ax.set_title(f"{carton} {loc.upper()}")
    ax.set_xlabel("Date")
    ax.set_ylabel("N targets with M epochs done")
    plt.savefig(f"{res_base}/{carton}_{loc}_N_epochs.pdf")
    plt.savefig(f"{res_base}/{carton}_{loc}_N_epochs.png")


def runAllCumulativeEpochs(base, plan, version=None):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    res_base = v_base + "byYearCartons"

    cachedFiles = glob.glob(f"{res_base}/*mjds*.fits")

    for f in cachedFiles:
        relative = f.split("/")[-1]
        carton = relative.split("_mjds")[0]
        loc = relative[-8:-5]
        cartonEpochCumulative(base, plan, version=version, carton=carton, loc=loc)
    

table_heads = """<html><head><meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
<h2>Cumulative Exposures</h2>
<p> The plots below show the cumulative epochs carton over time.</p>

<table><tbody>
<tr><td><h4>APO</h4></td> <td><h4>LCO </h4></td></tr>
"""

table_row ="""<tr><td><a href="{apo_pdf}"><img src="{apo_png}" width="450/"></a></td>
<td><a href="{lco_pdf}"><img src="{lco_png}" width="450/"></a></td></tr>"""

tail = """</tbody></table>
</body></html>"""

def cumulativeWebpage(base, plan, version=None):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    res_base = v_base + "byYearCartons"

    pngs = glob.glob(f"{res_base}/*N_epochs.png")

    html = table_heads

    for p in pngs:
        if "_apo_" in p:
            fname = p
        else:
            continue

        aname = fname.split("/")[-1]
        
        lname = aname.replace("apo", "lco")
        pics = {"apo_png": aname, "apo_pdf": aname.replace("png", "pdf"),
                "lco_png": lname, "lco_pdf": lname.replace("png", "pdf")}
        html += table_row.format(**pics)
    
    with open(res_base + "/cumalativeSummary.html", "w") as html_file:
        print(html + tail, file=html_file)
