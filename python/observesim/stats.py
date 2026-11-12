import os

import astropy.coordinates as coord
import astropy.units as u
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


def fieldCounts(v_base, plan, rs_base, loc="apo", idx=0):
    fields = fitsio.read(v_base + f"{plan}-{loc}-fields-{idx}.fits")

    tabulated = {c : [0, 0] for c in psuedo_cads}

    for f in fields:
        cad = matchCadence(f["cadence"])

        if cad is None:
            continue

        tabulated[cad][0] += f["nobservations"]
    
    alloc_file = rs_base + f"/{plan}/final/rsAllocationFinal-{plan}-{loc}.fits"
    if not os.path.isfile(alloc_file):
        altbase = "/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/sandbox/observesim/alternateInputs/"
        alloc_file = altbase + f"/{plan}/final/rsAllocationFinal-{plan}-{loc}.fits"
    allocation = fitsio.read(alloc_file)
    
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


def translate(name):
    if "bright" in name:
        return "bright"
    else:
        return name.replace("_no_apogee_skies", "").replace("_v2", "")


def lstSummary(v_base, plan, rs_base, loc="apo", idx=0):
    fields = fitsio.read(v_base + f"{plan}-{loc}-fields-{idx}.fits")

    obs = fitsio.read(v_base + f"{plan}-{loc}-observations-{idx}.fits")
    lst = fitsio.read(v_base + f"{plan}-{loc}-lst-{idx}.fits")
    alloc_file = rs_base + f"/{plan}/final/rsAllocationFinal-{plan}-{loc}.fits"
    if not os.path.isfile(alloc_file):
        altbase = "/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/sandbox/observesim/alternateInputs/"
        alloc_file = altbase + f"/{plan}/final/rsAllocationFinal-{plan}-{loc}.fits"
    alloc = fitsio.read(alloc_file)

    bins = np.arange(0, 25, 1)

    extra = []

    missed_bright = []
    missed_dark = []

    missed_bright_ra = []
    missed_dark_ra = []

    sum_missed = 0
    sum_missed2 = 0

    for f in fields:
        if f["base_priority"] < 0:
            extra.append(f["pk"])
            continue
        rs = alloc[f["pk"]]
        if rs["nallocated"] <= int(f["nobservations"]) or rs["nallocated"] < 1:
            continue
        sum_missed += (rs["nallocated"] - int(f["nobservations"]))
        debug = f"{f['pk']} {f['fieldid']} {rs['fieldid']} {f['racen']} {f['deccen']} {rs['racen']} {rs['deccen']}"
        assert rs['racen']-f['racen'] < 0.01 and rs['deccen']-f['deccen'] < 0.01, debug

        obs_idx = f["observations"][:int(f["nobservations"])]
        dark_lsts = [obs[i]["lst"] / 15 for i in obs_idx if obs[i]["skybrightness"] <= 0.35]
        bright_lsts = [obs[i]["lst"] / 15 for i in obs_idx if obs[i]["skybrightness"] > 0.35]
        dark_hist, bin_edges = np.histogram(dark_lsts, bins=bins)
        bright_hist, bin_edges = np.histogram(bright_lsts, bins=bins)
        dark_slots = rs["slots_exposures"][:,0]
        bright_slots = rs["slots_exposures"][:,1]

        unfinished = rs["nallocated"] - int(f["nobservations"])
        if unfinished > 1000 or rs["nallocated"] > 99:
            continue

        dark_diff = dark_slots - dark_hist
        w_missed = np.where(dark_diff > 0)[0]
        for w, l in zip(w_missed, dark_diff[w_missed]):
            # print(f["fieldid"], f["cadence"], w, l)
            missed_dark.extend([w for i in range(int(l))])

        bright_diff = bright_slots - bright_hist
        w_missed = np.where(bright_diff > 0)[0]
        for w, l in zip(w_missed, bright_diff[w_missed]):
            missed_bright.extend([w for i in range(int(l))])
        sum_missed2 += unfinished
        if "bright" in f["cadence"]:
            missed_bright_ra.extend([f["racen"] / 15 for i in range(unfinished)])
        else:
            # print(f["fieldid"], f["cadence"], f["racen"], unfinished)
            missed_dark_ra.extend([f["racen"] / 15 for i in range(unfinished)])

    print(f"TOTAL MISSED DESIGNS {loc}", sum_missed, sum_missed2)
    print(len(missed_dark_ra), "dark missed, ", len(missed_bright_ra), "bright missed")

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
    plt.title(f"Sim VS RoboStrategy LST ({sum_missed} missed, {len(unused)} unused)")
    plt.savefig(f"{v_base}/{plan}-{loc}-sim_vs_rs_lst-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-sim_vs_rs_lst-{idx}.pdf")

    plt.close()

    plt.figure()
    # plt.hist(dark["lst"] / 15, bins=bins, label="unused dark time", alpha=0.6)
    plt.hist(bright["lst"] / 15, bins=bins, label="unused bright time", alpha=0.6)

    # plt.hist(missed_dark_ra, bins=bins, label="missed dark RA", alpha=0.6)
    plt.hist(missed_bright_ra, bins=bins, label="missed bright RA", alpha=0.6)
    plt.legend(loc="best")
    plt.title(f"Sim VS RoboStrategy LST ({sum_missed} missed, {len(unused)} unused)")
    plt.savefig(f"{v_base}/{plan}-{loc}-sim_vs_lst_ra-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-sim_vs_lst_ra-{idx}.pdf")

    plt.close()

    unique_cadences = np.unique([translate(c) for c in lst["cadence"] if "x" in c])

    unique_cadences = list(unique_cadences)
    unique_cadences.append("dark_all")

    if loc.lower() == "apo":
        sjd_offset = 0.5
    else:
        sjd_offset = 0.4
    mjd = np.floor(lst["mjd"] - sjd_offset)

    unique_mjds = np.sort(np.unique(mjd))

    bright = list()
    dark = list()
    weather = list()
    skipped = list()
    bright_in_dark = list()
    weather_hours = list()
    onsky_hours = list()
    skipped_hours = list()
    bright_hours = list()
    dark_hours = list()
    over_alloc = list()
    overplan = list()

    cadences = {c: list() for c in unique_cadences}
    cadences_hours = {c: list() for c in unique_cadences}

    for m in unique_mjds:
        night = lst[np.where(mjd == m)]
        w_weather = np.where(night["weather"])
        weather.append(len(w_weather[0]))
        weather_hours.append(np.sum(night[w_weather]["duration"]) * 24)
        onsky = night[np.where(~night["weather"])]
        w_skipped = np.where(onsky["field_pk"] < 0)
        skipped.append(len(w_skipped[0]))
        skipped_hours.append(np.sum(night[w_skipped]["duration"]) * 24)
        ontarg = onsky[np.where(onsky["field_pk"] > 0)]
        onsky_hours.append(np.sum(ontarg["duration"]) * 24)
        b_mask = ["bright" in m for m in ontarg["mode"]]
        bright.append(len(np.where(b_mask)[0]))
        if len(b_mask) > 0:
            this_dark = ontarg[np.invert(b_mask)]
            # dark.append(len(np.where(np.invert(b_mask))[0]))
            dark.append(len(this_dark))
            dark_hours.append(np.sum(this_dark["duration"] * 24))
        else:
            dark.append(0)
            dark_hours.append(0)
        this_bright = ontarg[b_mask]
        bright_hours.append(np.sum(this_bright["duration"] * 24))
        b_in_d = this_bright[this_bright["bright"] < 0.35]
        bright_in_dark.append(len(b_in_d))
        if len(b_in_d):
            o = len(np.where(b_in_d["remaining_dark_lst"] < 0)[0])
            over_alloc.append(o)
        else:
            over_alloc.append(0)

        for cad in unique_cadences:
            cadences[cad].append(0)
            cadences_hours[cad].append(0)

        for o in ontarg:
            cad = translate(o["cadence"])
            cadences[cad][-1] += 1
            cadences_hours[cad][-1] += o["duration"] * 24
            if "dark" in cad:
                cadences["dark_all"][-1] += 1
                cadences_hours["dark_all"][-1] += o["duration"] * 24

        w_over = np.isin(ontarg["field_pk"], extra)
        overplan.append(len(np.where(w_over)[0]))

    skipped_fin = np.cumsum(skipped)[-1]
    weather_fin = np.cumsum(weather)[-1]
    bright_fin = np.cumsum(bright)[-1]
    dark_fin = np.cumsum(dark)[-1]
    bright_in_dark_fin = np.cumsum(bright_in_dark)[-1]
    either_fin = (np.cumsum(dark)+np.cumsum(bright))[-1]
    over_alloc_fin = np.cumsum(over_alloc)[-1]
    overplan_fin = np.cumsum(overplan)[-1]

    output_base = os.environ.get("OBSERVESIM_OUTPUT_BASE")
    time_file = os.path.join(output_base, f'time_avail_{loc}.csv')
    time_array = np.genfromtxt(time_file, names=True, delimiter=",", dtype=None, encoding="UTF-8")
    time_array.dtype

    if loc.lower() == "apo":
        dark_design = 23 / 60
        bright_design = 21 / 60
        dark_factor = 1.2
        bright_factor = 1.1
        weather_factor = 0.5
    else:
        dark_design = 24 / 60
        bright_design = 21 / 60
        dark_factor = 1.1
        bright_factor = 1.1
        weather_factor = 0.7

    start_date = np.min(unique_mjds)

    subset = time_array[time_array["mjd"] >= start_date]

    pred_bright = (subset["bright"] + subset["twilight"]) / bright_design * weather_factor
    adjusted_bright = np.cumsum(pred_bright) / bright_factor
    pred_dark = subset["dark"] / dark_design * weather_factor
    adjusted_dark = np.cumsum(pred_dark) / dark_factor

    realistic_total = adjusted_bright + adjusted_dark

    inter_mjd, i_avail, i_lst = np.intersect1d(subset['mjd'], unique_mjds, return_indices=True)

    pred_skip = realistic_total[i_avail] - (np.cumsum(dark)+np.cumsum(bright))[i_lst]

    pred_skip_fin = int(pred_skip[-1])
    adjusted_bright_fin = int(adjusted_bright[-1])
    adjusted_dark_fin = int(adjusted_dark[-1])
    realistic_total_fin = int(realistic_total[-1])

    plt.plot(unique_mjds, np.cumsum(bright_in_dark), c="tab:brown", label=f"bright in dark ({bright_in_dark_fin})")
    plt.plot(unique_mjds, np.cumsum(skipped), c="tab:red", label=f"nothing ({skipped_fin})")
    plt.plot(unique_mjds, np.cumsum(weather), c="tab:gray", label=f"weather ({weather_fin})")
    plt.plot(unique_mjds, np.cumsum(bright), c="tab:olive", label=f"bright ({bright_fin})")
    plt.plot(unique_mjds, np.cumsum(dark), c="tab:blue", label=f"dark ({dark_fin})")
    plt.plot(unique_mjds, np.cumsum(over_alloc), c="tab:green", label=f"encroaching ({over_alloc_fin})")
    plt.plot(unique_mjds, np.cumsum(overplan), c="tab:pink", label=f"overplan ({overplan_fin})")
    plt.plot(unique_mjds, np.cumsum(dark)+np.cumsum(bright), c="tab:purple", label=f"bright+dark ({either_fin})")
    plt.axhline(0, c="k", linewidth=0.5)
    plt.legend(fontsize=10)
    plt.title("Cumulative time usage (exps)")
    plt.xlabel("MJD")
    plt.ylabel("N")
    plt.savefig(f"{v_base}/{plan}-{loc}-cumaltive_time-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-cumaltive_time-{idx}.pdf")

    plt.plot(inter_mjd, pred_skip, c="tab:red", linestyle="--", label=f"nothing pred ({pred_skip_fin})")
    plt.plot(subset["mjd"], adjusted_bright, linestyle="--", c="tab:olive", label=f"bright avail ({adjusted_bright_fin})")
    plt.plot(subset["mjd"], adjusted_dark, linestyle="--", c="tab:blue", label=f"dark avail ({adjusted_dark_fin})")
    plt.plot(subset["mjd"], realistic_total, linestyle="--", c="tab:purple", label=f"bright+dark avail ({realistic_total_fin})")
    
    plt.savefig(f"{v_base}/{plan}-{loc}-cumaltive_time_w_model-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-cumaltive_time_w_model-{idx}.pdf")
    plt.close()

    sim_all = np.cumsum(weather_hours) + np.cumsum(onsky_hours) + np.cumsum(skipped_hours)

    avail_all = np.cumsum(subset["bright"] + subset["dark"] + subset["twilight"])
    pred_weather = avail_all * (1-weather_factor)

    onsky_days_fin = np.sum(onsky_hours)
    skipped_fin = np.sum(skipped_hours)
    weather_fin = np.sum(weather_hours)
    dark_fin = np.sum(dark_hours) / dark_fin * 60
    bright_fin = np.sum(bright_hours) / bright_fin * 60

    plt.plot(unique_mjds, np.cumsum(onsky_hours), c="tab:purple", label=f"on sky ({int(onsky_days_fin)})")
    plt.plot(unique_mjds, np.cumsum(skipped_hours), c="tab:red", label=f"nothing ({int(skipped_fin)})")
    plt.plot(unique_mjds, np.cumsum(weather_hours), c="tab:gray", label=f"weather ({int(weather_fin)})")
    plt.plot(subset["mjd"], pred_weather, c="tab:gray", linestyle="--")
    plt.plot(unique_mjds, np.cumsum(dark_hours), c="tab:blue", label=f"dark per des {dark_fin:.1f}")
    plt.plot(unique_mjds, np.cumsum(bright_hours), c="tab:olive", label=f"bright per des {bright_fin:.1f}")
    plt.plot(unique_mjds, sim_all, c="tab:purple", label=f"sim all ({int(sim_all[-1])})")
    plt.plot(subset["mjd"], avail_all, c="tab:purple", linestyle="--", label=f"available ({int(avail_all[-1])})")
    # plt.axhline(0, c="k", linewidth=0.5)
    plt.legend(fontsize=10)
    plt.title("Cumulative time usage in hours")
    plt.xlabel("MJD")
    plt.ylabel("hours")

    plt.savefig(f"{v_base}/{plan}-{loc}-cumaltive_days-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-cumaltive_days-{idx}.pdf")
    plt.close()

    for name, hist in cadences.items():
        plt.plot(unique_mjds, np.cumsum(hist), label=f"{name} ({np.cumsum(hist)[-1]})")
    plt.legend(fontsize=10)
    plt.title("Cumulative designs by cadence")
    plt.xlabel("MJD")
    plt.ylabel("N")

    plt.savefig(f"{v_base}/{plan}-{loc}-cumaltive_cadence-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-cumaltive_cadence-{idx}.pdf")
    plt.close()

    for name, hist in cadences_hours.items():
        plt.plot(unique_mjds, np.cumsum(hist), label=f"{name} ({int(np.cumsum(hist)[-1])})")
    plt.legend(fontsize=10)
    plt.title("Cumulative hours by cadence")
    plt.xlabel("MJD")
    plt.ylabel("hours")

    plt.savefig(f"{v_base}/{plan}-{loc}-cumaltive_cadence_hours-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-cumaltive_cadence_hours-{idx}.pdf")
    plt.close()


def fieldCompletion(v_base, plan, loc="apo", idx=0):
    sim_data = fitsio.read(v_base + f"{plan}-{loc}-fields-{idx}.fits")

    missed = sim_data["nfilled"] - sim_data["nobservations"]
    w_xtra = np.where(missed < -1)
    # for f in sim_data[w_xtra]:
    #     print(f["nfilled"], f["nobservations"], f["cadence"])

    # w_prob = np.where(np.logical_and(missed > 0, sim_data["nfilled"] < 200))
    w_prob = np.where((missed > 0) & (sim_data["nfilled"] < 200) & (sim_data["base_priority"] >= 0))

    plt.figure()
    bins = np.arange(1, 20, 1)
    plt.hist(missed[w_prob], bins=bins)
    plt.title("Fields with N missing observations")
    plt.savefig(f"{v_base}/{plan}-{loc}-missing_obs_count-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-missing_obs_count-{idx}.pdf")
    plt.close()

    problem = sim_data[w_prob]
    # for p in problem[np.argsort(problem["cadence"])]:
    #     print(p["fieldid"], p["cadence"], p["nfilled"], p["nobservations"])
    plt.figure()
    fig, ax = plt.subplots()
    im = ax.scatter(problem["nfilled"], missed[w_prob],
               vmax=10, c=problem["nobservations"])
    # ax.set_xlim(0,20)
    # ax.set_ylim(0,20)
    ax.set_xlabel("N designs planned")
    ax.set_ylabel("N designs missed")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("N observations")
    plt.savefig(f"{v_base}/{plan}-{loc}-missing_obs_v_planned-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-missing_obs_v_planned-{idx}.pdf")
    plt.close()

    ra = coord.Angle(-(problem['racen']+90)*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(problem['deccen']*u.degree)

    # idx_sort = np.argsort(problem['racen'])
    # for indx in idx_sort:
    #     if loc != "apo" or ra[indx].degree < 60:
    #         continue
    #     print(int(ra[indx].degree), dec[indx], problem["cadence"][indx], missed[w_prob][indx])

    f = plt.figure()
    
    ax = plt.subplot(111, projection='mollweide')
    im = ax.scatter(ra.radian, dec.radian, s=3, vmax=10, c=missed[w_prob])
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("N missed")
    plt.savefig(f"{v_base}/{plan}-{loc}-missing_obs_map-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-missing_obs_map-{idx}.pdf")
    plt.close()
    
    ra = coord.Angle(-(sim_data['racen']+90)*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(sim_data['deccen']*u.degree)

    f = plt.figure()
    
    ax = plt.subplot(111, projection='mollweide')
    im = ax.scatter(ra.radian, dec.radian, s=3, c=sim_data["base_priority"])
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Priority")
    plt.savefig(f"{v_base}/{plan}-{loc}-priority_map-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-priority_map-{idx}.pdf")
    plt.close()

    f = plt.figure()
    ax = plt.subplot(111, projection='mollweide')
    im = ax.scatter(ra.radian, dec.radian, s=3, vmax=3, c=sim_data["nobservations"])
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("N done")
    plt.savefig(f"{v_base}/{plan}-{loc}-done_map-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-done_map-{idx}.pdf")
    plt.close()
    f = plt.figure()

    
    ra = coord.Angle((sim_data['racen'])*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(sim_data['deccen']*u.degree)

    ax = plt.subplot(111, projection='mollweide')
    im = ax.scatter(ra.radian, dec.radian, s=3, vmax=3, c=sim_data["nobservations"])
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("N done")
    plt.savefig(f"{v_base}/{plan}-{loc}-done_map-realra-{idx}.png")
    plt.savefig(f"{v_base}/{plan}-{loc}-done_map-realra-{idx}.pdf")
    plt.close()



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

    missing = """
    <table><tbody>
    <tr><td>
    <a href="{plan}-apo-cumaltive_cadence-{idx}.pdf"><img src="{plan}-apo-cumaltive_cadence-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-cumaltive_cadence-{idx}.pdf"><img src="{plan}-lco-cumaltive_cadence-{idx}.png" width="600px/"> </a>
    </td></tr>
    <tr><td>
    <a href="{plan}-apo-cumaltive_cadence_hours-{idx}.pdf"><img src="{plan}-apo-cumaltive_cadence_hours-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-cumaltive_cadence_hours-{idx}.pdf"><img src="{plan}-lco-cumaltive_cadence_hours-{idx}.png" width="600px/"> </a>
    </td></tr>
    <tr><td>
    <a href="{plan}-apo-cumaltive_time-{idx}.pdf"><img src="{plan}-apo-cumaltive_time-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-cumaltive_time-{idx}.pdf"><img src="{plan}-lco-cumaltive_time-{idx}.png" width="600px/"> </a>
    </td></tr>
    <tr><td>
    <a href="{plan}-apo-cumaltive_days-{idx}.pdf"><img src="{plan}-apo-cumaltive_days-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-cumaltive_days-{idx}.pdf"><img src="{plan}-lco-cumaltive_days-{idx}.png" width="600px/"> </a>
    </td></tr>
    <tr><td>
    <a href="{plan}-apo-sim_vs_rs_lst-{idx}.pdf"><img src="{plan}-apo-sim_vs_rs_lst-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-sim_vs_rs_lst-{idx}.pdf"><img src="{plan}-lco-sim_vs_rs_lst-{idx}.png" width="600px/"> </a>
    </td></tr>
    <tr><td>
    <a href="{plan}-apo-missing_obs_count-{idx}.pdf"><img src="{plan}-apo-missing_obs_count-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-missing_obs_count-{idx}.pdf"><img src="{plan}-lco-missing_obs_count-{idx}.png" width="600px/"> </a>
    </td></tr>
    <tr><td>
    <a href="{plan}-apo-missing_obs_v_planned-{idx}.pdf"><img src="{plan}-apo-missing_obs_v_planned-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-missing_obs_v_planned-{idx}.pdf"><img src="{plan}-lco-missing_obs_v_planned-{idx}.png" width="600px/"> </a>
    </td></tr>
    <tr><td>
    <a href="{plan}-apo-missing_obs_map-{idx}.pdf"><img src="{plan}-apo-missing_obs_map-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-missing_obs_map-{idx}.pdf"><img src="{plan}-lco-missing_obs_map-{idx}.png" width="600px/"> </a>
    </td></tr>
    <tr><td>
    <a href="{plan}-apo-priority_map-{idx}.pdf"><img src="{plan}-apo-priority_map-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-priority_map-{idx}.pdf"><img src="{plan}-lco-priority_map-{idx}.png" width="600px/"> </a>
    </td></tr>
    </tbody></table>
    """.format(plan=plan, idx=idx)

    missing_html = "<html><body>" + missing + "</body></html>"

    with open(f"{v_base}/incomplete-{plan}-{idx}.html", "w") as webPage:
        print(missing_html, file=webPage)

    next_tab = """<table><tbody>

    <tr> <th>cadence</th> <th>RS planned APO field_exps</th> <th>Done APO field_exps</th>
    <th>RS planned LCO field_exps</th> <th>Done LCO fiefield_expslds</th></tr>
    """

    field_table_row = """<tr><td>{cadence}</td> <td>{rs_apo}</td> <td>{done_apo}</td>
    <td>{rs_lco}</td> <td>{done_lco}</td> </tr>"""

    tail = f"""</tbody></table>

    <a href="{plan}-apo-sim-v-theory-{idx}.pdf"><img src="{plan}-apo-sim-v-theory-{idx}.png" width="600px/"> </a>
    <a href="{plan}-lco-sim-v-theory-{idx}.pdf"><img src="{plan}-lco-sim-v-theory-{idx}.png" width="600px/"> </a>
    """
    tail += missing + "\n \n </body></html>"

    # carton_counts = countCartons(v_base, plan, rs_base)
    tabulated_apo = fieldCounts(v_base, plan, rs_base, loc="apo", idx=idx)
    tabulated_lco = fieldCounts(v_base, plan, rs_base, loc="lco", idx=idx)
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

    fieldCompletion(v_base, plan, loc="apo", idx=idx)
    fieldCompletion(v_base, plan, loc="lco", idx=idx)

    with open(f"{v_base}/obsSimStats-{plan}-{idx}.html", "w") as webPage:
        print(header + tail, file=webPage)


def multiQuickSummary(base, plan, rs_base, version=None, hist_mjd=None, n=1, loc=None):
    if loc is not None:
        if version is not None:
            v_base = os.path.join(base, version)
            v_base += "/"
        else:
            v_base = os.path.join(base, plan)
            v_base += "/"
        for i in range(n):
            fieldCounts(v_base, plan, rs_base, loc=loc)
            cumulativeDesigns(v_base, plan, rs_base, loc=loc, idx=i, hist_mjd=hist_mjd)
            lstSummary(v_base, plan, rs_base, loc=loc, idx=i)
            fieldCompletion(v_base, plan, loc=loc, idx=i)
        return
    for i in range(n):
        quickSummary(base, plan, rs_base, version=version, idx=i, hist_mjd=hist_mjd)
