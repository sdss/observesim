import os
import numpy as np
import fitsio
import matplotlib.pyplot as plt


psuedo_cads = ["dark_10x4_4yr",
               "dark_174x8",
               "dark_100x8",
               "dark_2x1",
               "dark_2x2",
               "dark_2x4",
               "bright_x1",
               "bright_x2",
               "bright_x4"]

def countCartons(v_base, plan, rs_base, idx=0):
    obs_data_apo = fitsio.read(v_base + f"obsTargets-{plan}-apo-{idx}.fits",
                           columns=["carton"])
    obs_data_lco = fitsio.read(v_base + f"obsTargets-{plan}-lco-{idx}.fits",
                           columns=["carton"])
    comp_data = fitsio.read(rs_base + 
                            f"/{plan}/final/rsCompletenessFinal-{plan}-both.fits",
                            columns=["program", "carton",
                                     "satisfied_apo", "satisfied_lco"])

    cartons = np.unique(comp_data["carton"])

    # header = "carton, RS satisfied obs, observed\n"
    
    carton_counts = list()

    for c in cartons:
        if "ops" in c:
            continue
        print("summarizing: ", c)
        obs_carton_apo = obs_data_apo[np.where(obs_data_apo["carton"] == c)]
        obs_carton_lco = obs_data_lco[np.where(obs_data_lco["carton"] == c)]
        comp_carton = comp_data[np.where(comp_data["carton"] == c)]
        comp_apo = comp_carton[np.where(comp_carton["satisfied_apo"])]
        comp_lco = comp_carton[np.where(comp_carton["satisfied_lco"])]
        carton_counts.append(
            {"carton": c,
            "rs_apo": len(comp_apo),
            "done_apo": len(obs_carton_apo),
            "rs_lco": len(comp_lco),
            "done_lco":len(obs_carton_lco),}
        )

    return carton_counts


def matchCadence(cadence):
    base = cadence.split("_v")[0]
    for cad in psuedo_cads:
        if base == cad:
            return cad
    if "x1" in base:
        return "bright_x1"
    elif "x2" in base:
        return "bright_x2"
    elif "x4" in base:
        return "bright_x4"
    else:
        # print(f"missing cadence! {base}")
        return None


def fieldCounts(v_base, plan, rs_base, loc="apo"):
    fields = fitsio.read(v_base + f"{plan}-{loc}-fields-0.fits")

    tabulated = {c : [0, 0] for c in psuedo_cads}

    for f in fields:
        cad = matchCadence(f["cadence"])

        if cad is None:
            continue

        tabulated[cad][0] += f["nobservations"]
    
    allocation = fitsio.read(rs_base + 
                 f"/{plan}/final/rsAllocationFinal-{plan}-{loc}.fits")
    
    if "nfilled" in allocation.dtype.names:
        nfilled = allocation["nfilled"]
    elif "nallocated_full" in allocation.dtype.names:
        nfilled = allocation["nallocated_full"]
    else:
        print("WARN: strange rsAllocation format, estimating nfilled")
        nfilled = np.sum(np.sum(fits_dat["slots_exposures"], axis=1),
                         axis=1, dtype=int)

    for i, f in enumerate(allocation):
        cad = matchCadence(f["cadence"])

        if cad is None:
            continue
           
        tabulated[cad][1] += nfilled[i]
    
    return tabulated


def cumulativeDesigns(v_base, plan, rs_base, loc="apo", idx=0, hist_mjd=None):

    time_file = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + f"/etc/time_avail_{loc}.csv"
    time_array = np.genfromtxt(time_file, names=True, delimiter=",", dtype=None, encoding="UTF-8")

    if loc == "apo":
        weather = 0.5
        dark_design = 23 / 60
        bright_design = 18 / 60

        bright_factor = 1.1
        dark_factor = 1.2
    else:
        weather = 0.7
        dark_design = 23 / 60
        bright_design = 21 / 60

        bright_factor = 1.1
        dark_factor = 1.4

    max_bright = time_array["cum_bright"] / bright_design / bright_factor * weather
    max_dark = time_array["cum_dark"] / dark_design / dark_factor * weather

    sim_data = fitsio.read(v_base + f"{plan}-{loc}-fields-{idx}.fits")

    obs_data = fitsio.read(v_base + f"{plan}-{loc}-observations-{idx}.fits")

    dark = list()
    bright = list()
    for f in sim_data:
        obs_idx = f["observations"][:int(f["nobservations"])]
        mjds = [obs_data[i]["mjd"] for i in obs_idx]

        cad = f["cadence"]
        if "bright" in cad:
            bright.extend(mjds)
        else:
            dark.extend(mjds)

    mjd_start = np.min(bright + dark)
    mjd_end = np.max(bright + dark) + 2

    mjds = np.arange(mjd_start, mjd_end, 1)

    bright = np.array(bright)
    dark = np.array(dark)

    hist_offset_bright = 0
    hist_offset_dark = 0

    cum_bright = list()
    cum_dark = list()
    for m in mjds:
        w_bright = np.where(bright < m)
        bright_now = len(w_bright[0])
        cum_bright.append(bright_now)
        w_dark = np.where(dark < m)
        dark_now = len(w_dark[0])
        cum_dark.append(dark_now)
        w_now = np.where(time_array["mjd"] == int(m))
        if m <= hist_mjd and len(w_now[0]):
            offset_bright_now = int(max_bright[w_now]) - bright_now
            if offset_bright_now > hist_offset_bright:
                hist_offset_bright = offset_bright_now
            offset_dark_now = int(max_dark[w_now]) - dark_now
            if offset_dark_now > hist_offset_dark:
                hist_offset_dark = offset_dark_now
            max_bright[w_now] = bright_now
            max_dark[w_now] = dark_now
        else:
            max_bright[w_now] -= hist_offset_bright
            max_dark[w_now] -= hist_offset_dark

    cum_bright = np.array(cum_bright)
    cum_dark = np.array(cum_dark)

    plt.figure(figsize=(8,6))
    plt.plot(time_array["mjd"], max_bright, c="r", linestyle="--", label="theoretical bright")
    plt.plot(time_array["mjd"], max_dark, c="b", linestyle="--", label="theoretical dark")
    plt.plot(time_array["mjd"], max_dark+max_bright, c="k", linestyle="--", label="theoretical total")

    plt.plot(mjds, cum_bright, c="r", linestyle="-", label="sim bright")
    plt.plot(mjds, cum_dark, c="b", linestyle="-", label="sim dark")
    plt.plot(mjds, cum_bright+cum_dark, c="k", linestyle="-", label="sim total")

    plt.xlabel("MJD")
    plt.ylabel("N Designs")
    plt.title(f"{loc}")
    plt.legend(loc="best")
    plt.savefig(f"{v_base}/{plan}-{loc}-sim-v-theory-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-sim-v-theory-{idx}.pdf")
    plt.close()


def lstSummary(v_base, plan, rs_base, loc="apo", idx=0):
    fields = fitsio.read(v_base + f"{plan}-{loc}-fields-{idx}.fits")

    obs = fitsio.read(v_base + f"{plan}-{loc}-observations-{idx}.fits")
    lst = fitsio.read(v_base + f"{plan}-{loc}-lst-{idx}.fits")
    alloc = fitsio.read(rs_base + f"/{plan}/final/rsAllocationFinal-{plan}-{loc}.fits")

    bins = np.arange(0, 25, 1)

    missed_bright = []
    missed_dark = []

    sum_missed = 0

    for f in fields:
        rs = alloc[f["pk"]]
        if rs["nallocated"] <= int(f["nobservations"]) or rs["nallocated"] < 1:
            continue
        sum_missed += (rs["nallocated"] - int(f["nobservations"]))
        assert rs['racen']-f['racen'] < 0.01 and rs['deccen']-f['deccen'] < 0.01

        obs_idx = f["observations"][:int(f["nobservations"])]
        dark_lsts = [obs[i]["lst"] / 15 for i in obs_idx if obs[i]["skybrightness"] <= 0.35]
        bright_lsts = [obs[i]["lst"] / 15 for i in obs_idx if obs[i]["skybrightness"] > 0.35]
        dark_hist, bin_edges = np.histogram(dark_lsts, bins=bins)
        bright_hist, bin_edges = np.histogram(bright_lsts, bins=bins)
        dark_slots = rs["slots_exposures"][:,0]
        bright_slots = rs["slots_exposures"][:,1]

        dark_diff = dark_slots - dark_hist
        w_missed = np.where(dark_diff > 0)[0]
        for w, l in zip(w_missed, dark_diff[w_missed]):
            missed_dark.extend([w for i in range(int(l))])

        bright_diff = bright_slots - bright_hist
        w_missed = np.where(bright_diff > 0)[0]
        for w, l in zip(w_missed, bright_diff[w_missed]):
            missed_bright.extend([w for i in range(int(l))])

    print(f"TOTAL MISSED DESIGNS {loc}", sum_missed)

    unused = lst[np.where(np.logical_and(lst["field_pk"] < 0, ~lst["weather"]))]
    bins = np.arange(0, 25, 1)
    dark = unused[np.where(unused["bright"] <= 0.35)]
    bright = unused[np.where(unused["bright"] > 0.35)]
    plt.figure()
    plt.hist(dark["lst"] / 15, bins=bins, label="unused dark time", alpha=0.6)
    plt.hist(bright["lst"] / 15, bins=bins, label="unused bright time", alpha=0.6)

    plt.hist(missed_dark, bins=bins, label="RS missed dark", alpha=0.6)
    plt.hist(missed_bright, bins=bins, label="RS missed bright", alpha=0.6)
    plt.legend(loc="best")
    plt.title("Sim VS RS LST")
    plt.savefig(f"{v_base}/{plan}-{loc}-sim_vs_rs_lst-{idx}.png")


def quickSummary(base, plan, rs_base, version=None, idx=0, hist_mjd=None):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"
    print(f"!! {version}, {v_base}")
    header = """<html><head><meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
    <style>
        table, th, tr, td {{border: 2px solid black}}
    </style>
    </head><body>
    """

    not_used_currently = """
    <table><tbody>

    <tr> <th>carton</th> <th>RS planned APO targs</th> <th>Done APO targs</th>
    <th>RS planned LCO targs</th> <th>Done LCO targs</th></tr>
    """

    carton_table_row = """<tr><td>{carton}</td> <td>{rs_apo}</td> <td>{done_apo}</td>
    <td>{rs_lco}</td> <td>{done_lco}</td> </tr>"""

    # next_tab = "</tbody></table>"

    next_tab = """<table><tbody>

    <tr> <th>cadence</th> <th>RS planned APO field_exps</th> <th>Done APO field_exps</th>
    <th>RS planned LCO field_exps</th> <th>Done LCO fiefield_expslds</th></tr>
    """

    field_table_row = """<tr><td>{cadence}</td> <td>{rs_apo}</td> <td>{done_apo}</td>
    <td>{rs_lco}</td> <td>{done_lco}</td> </tr>"""

    tail = f"""</tbody></table>

    <a href="{plan}-apo-sim-v-theory.pdf"><img src="{plan}-apo-sim-v-theory-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-sim-v-theory.pdf"><img src="{plan}-lco-sim-v-theory-{idx}.png" width="600px/"> </a>
    </body></html>"""

    # carton_counts = countCartons(v_base, plan, rs_base)
    tabulated_apo = fieldCounts(v_base, plan, rs_base, loc="apo")
    tabulated_lco = fieldCounts(v_base, plan, rs_base, loc="lco")
    cumulativeDesigns(v_base, plan, rs_base, loc="apo", idx=idx, hist_mjd=hist_mjd)
    cumulativeDesigns(v_base, plan, rs_base, loc="lco", idx=idx, hist_mjd=hist_mjd)
    lstSummary(v_base, plan, rs_base, loc="apo", idx=idx)
    lstSummary(v_base, plan, rs_base, loc="lco", idx=idx)

    # for c in carton_counts:
    #     header += carton_table_row.format(**c)
    
    header += next_tab

    total_apo = [0, 0]
    total_lco = [0, 0]
    for k, v in tabulated_apo.items():
        total_apo[0] += v[0]
        total_apo[1] += v[1]
    for k, v in tabulated_lco.items():
        total_lco[0] += v[0]
        total_lco[1] += v[1]

    tabulated_apo["total"] = total_apo
    tabulated_lco["total"] = total_lco

    for k in tabulated_apo.keys():
        row_dict = {
            "cadence": k,
            "rs_apo": tabulated_apo[k][1],
            "done_apo": tabulated_apo[k][0],
            "rs_lco": tabulated_lco[k][1],
            "done_lco": tabulated_lco[k][0],
        }
        header += field_table_row.format(**row_dict)
    
    with open(f"{v_base}/obsSimStats-{plan}-{idx}.html", "w") as webPage:
        print(header + tail, file=webPage)


def multiQuickSummary(base, plan, rs_base, version=None, hist_mjd=None):
    for i in range(4):
        quickSummary(base, plan, rs_base, version=version, idx=i, hist_mjd=hist_mjd)
