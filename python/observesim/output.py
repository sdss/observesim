import time
import os
import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import fitsio
import yaml
import healpy as hp
from PyAstronomy.pyasl.asl.astroTimeLegacy import daycnv
import astropy.coordinates as coord
import astropy.units as u

from roboscheduler.fields import epochs_completed

__all__ = ["cumulativePlot", "doHist", "plotTargMetric", "writeWebPage",
           "countFields", "spiders_area_for_program_time",
           "cartonCoverageByYearMaps", "cartonCoverageCumulative",
           "runAllCumulativeEpochs"]


def read_field(fname, alloc, exp_to_mjd):
    # fetch info, match mjd to exp

    w_names = fitsio.read(fname, ext=1)

    w_idx = fitsio.read(fname, ext=2)

    catalog_ids = list()
    cadences = list()
    cartons = list()
    programs = list()
    target_pks = list()
    carton_pks = list()
    c2t_pks = list()
    categories = list()
    ras = list()
    decs = list()
    mjds = list()

    for m, i in zip(exp_to_mjd, range(alloc["iexpst"], alloc["iexpnd"]+1)):
        # zip will end when exp_to_mjd ends if it is shorter than
        # nplanned (i.e. the range)
        if len(w_idx["equivRobotID"].shape) == 1:
            assert len(exp_to_mjd) == 1, "obs len != plan len"
            w_assigned = np.where(w_idx["equivRobotID"] != -1)
        else:
            w_assigned = np.where(w_idx["equivRobotID"][:, i] != -1)

        # assert len(w_assigned[0]) <= 500, "more ids than robots!"

        catalog_ids.extend(list(w_names["catalogid"][w_assigned]))
        cadences.extend(list(w_names["cadence"][w_assigned]))
        cartons.extend(list(w_names["carton"][w_assigned]))
        programs.extend(list(w_names["program"][w_assigned]))
        target_pks.extend(list(w_names["target_pk"][w_assigned]))
        carton_pks.extend(list(w_names["carton_pk"][w_assigned]))
        c2t_pks.extend(list(w_names["carton_to_target_pk"][w_assigned]))
        categories.extend(list(w_names["category"][w_assigned]))
        ras.extend(list(w_names["ra"][w_assigned]))
        decs.extend(list(w_names["dec"][w_assigned]))
        mjds.extend([m for i in w_assigned[0]])

    return (catalog_ids, cadences, cartons, programs, target_pks,
            carton_pks, c2t_pks, categories, ras, decs, mjds)


def countFields(res_base, rs_base, plan, version=None, loc="apo", N=0, save=True):
    """Create obsTargets summary file, combining robostrategy planned visits
       with observesim visits.
    """
    if version is not None:
        v_base = os.path.join(res_base, version)
        v_base += "/"
    else:
        v_base = os.path.join(res_base, plan)
        v_base += "/"

    allocation = fitsio.read(rs_base + "{plan}/final/rsAllocationFinal-{plan}-{loc}.fits".format(plan=plan, loc=loc))

    sim_data = fitsio.read(v_base + f"{plan}-{loc}-fields-{N}.fits")
    obs_data = fitsio.read(v_base + f"{plan}-{loc}-observations-{N}.fits")

    # prep out struct
    all_targs = list()
    all_cads = list()
    all_mjds = list()
    all_field_ids = list()
    all_field_pks = list()
    all_programs = list()
    all_target_pks = list()
    all_cartons = list()
    all_carton_pks = list()
    all_c2t_pks = list()
    all_categories = list()
    all_ras = list()
    all_decs = list()

    start_time = time.time()

    print("matching fields for {} {}".format(plan, loc))
    for f in sim_data:
        obs_idx = f["observations"][:int(f["nobservations"])]
        mjds = [obs_data[i]["mjd"] for i in obs_idx]
        cad = f["cadence"]

        # field_pk is an index into rsAllocationFinal for each observatory
        # so definitely keep doing that in roboscheduler
        alloc = allocation[f["pk"]]

        fname = "{plan}/final/rsFieldAssignmentsFinal-{plan}-{loc}-{fieldid}.fits"
        fname = rs_base + fname.format(plan=plan, loc=loc, fieldid=f["fieldid"])

        cat_ids, cadences, cartons, programs, target_pks, carton_pks, c2t_pks,\
            categories, ras, decs, targ_mjds = read_field(fname, alloc, mjds)
        all_targs.extend(cat_ids)
        all_cads.extend(cadences)
        all_mjds.extend(targ_mjds)
        all_field_ids.extend([f["fieldid"] for i in cat_ids])
        all_field_pks.extend([f["pk"] for i in cat_ids])
        all_programs.extend(programs)
        all_target_pks.extend(target_pks)
        all_cartons.extend(cartons)
        all_carton_pks.extend(carton_pks)
        all_c2t_pks.extend(c2t_pks)
        all_categories.extend(categories)
        all_ras.extend(ras)
        all_decs.extend(decs)

        # print(f["pk"], f["fieldid"], "GOT", len(ids), len(all_targs))

    assert len(all_targs) == len(all_cads), "targ != cad!!"
    assert len(all_targs) == len(all_mjds), "targ != mjd!!"
    assert len(all_targs) == len(all_field_ids), "targ != field!!"


    # completeness = fitsio.read(rs_base + "{plan}/rsCompleteness-{plan}-{loc}.fits".format(
    #                                       loc=loc, plan=plan),
    #                            columns=["catalogid", "carton", "program", "assigned", "cadence", "ra", "dec"])
    # completeness = completeness[np.where(completeness["catalogid"] != -1)]

    # comp_dict = dict()
    # for c in completeness:
    #     comp_dict[c["catalogid"]] = (c["catalogid"], c["carton"], c["program"], c["assigned"], c["cadence"], c["ra"], c["dec"])

    # for k, c in zip(all_targs, all_cads):
    #     targetid, carton, program, got, cadence, ra, dec = comp_dict[k]
    #     targ_ids.append(targetid)
    #     programs.append(program)
    #     cartons.append(carton)
    #     # gots.append(got)
    #     ras.append(ra)
    #     decs.append(dec)
    #     # if cadence != c:
    #     #     print("field cad: {}, targ cad: {}".format(c, cadence))
    #     #     print("pk: {}, targ id: {}".format(k, targetid))

    dtype = [('catalogid', np.int64),
             ('target_pk', np.int64),
             ('carton_to_target_pk', np.int64),
             ('carton_pk', np.int64),
             ('cadence', np.dtype('a40')),
             ('program', np.dtype('a40')),
             ('carton', np.dtype('a40')),
             ('category', np.dtype('a40')),
             ('field_id', np.int32),
             ('field_pk', np.int32),
             # ('assigned', np.int32),
             ('obs_mjd', np.float64),
             ('ra', np.float64),
             ('dec', np.float64)]
    obs_targs = np.zeros(len(all_targs), dtype=dtype)

    obs_targs["catalogid"] = all_targs
    obs_targs["cadence"] = all_cads
    obs_targs["field_id"] = all_field_ids
    obs_targs["field_pk"] = all_field_pks
    obs_targs["obs_mjd"] = all_mjds
    # obs_targs["assigned"] = gots
    obs_targs["program"] = all_programs
    obs_targs["carton"] = all_cartons
    obs_targs["category"] = all_categories
    obs_targs["target_pk"] = all_target_pks
    obs_targs["carton_to_target_pk"] = all_c2t_pks
    obs_targs["carton_pk"] = all_carton_pks
    obs_targs["ra"] = all_ras
    obs_targs["dec"] = all_decs

    # catch those -1 obs mjds!!! Oops
    # UPDATE: with rs*Final updates, this should be unnecessary, but harmless
    where_not_observed = np.where(obs_targs["obs_mjd"] > 5e4)
    cleaned_targs = obs_targs[where_not_observed]

    # fitsio.write(v_base + "allObsTargets-{plan}-{loc}.fits".format(plan=plan, loc=loc),
    #              obs_targs, clobber=True)

    if save:
        fitsio.write(v_base + f"obsTargets-{plan}-{loc}-{N}.fits",
                     cleaned_targs, clobber=True)

    return obs_targs


header = """<html><head><meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
<style>
    table, th, tr, td {{border: 1px solid black}}
    .bg-green  {{background-color:green; color:black;}}
    .bg-orange {{background-color:orange; color:black;}}
    .bg-red    {{background-color:red; color:black;}}
</style>
</head><body>

<h1>Observesim: {plan}</h1>
<p>This page summarizes the observesim results for the {plan} robostrategy runs <a href="https://data.sdss.org/sas/sdss5/sandbox/robostrategy/allocations/{plan}">found here</a></p>

<p>A video animation of this simulation is available! <a href="moviePngsAllSky/sdss5_sim.mp4" download>Download</a></p>

<h2>Summary</h2>
<p> The results for each observatory are summarized below. The average completion % for each field cadence class is shown, as well as the number of exposures falling into each class (note the log scale). </p>"""

bar_plots = """<a href="{plan}-apo-cadence_bar.pdf"><img src="{plan}-apo-cadence_bar.png" width="600px/"> </a>
<a href="{plan}-lco-cadence_bar.pdf"><img src="{plan}-lco-cadence_bar.png" width="600px/"> </a>"""

fits_out = """<h2> Output files: </h2>
<p> Each input field from the rsAllocation file is recorded in a fields file for  <a href="{plan}-apo-fields-0.fits">APO</a> and <a href="{plan}-lco-fields-0.fits">LCO</a>.
The fields files contain indices into an observations to record info for each observation. <a href="{plan}-apo-observations-0.fits"> APO observations</a>, and <a href="{plan}-lco-observations-0.fits"> LCO observations</a>. </p>
<p> A record of the priority decisions made by roboscheduler is now <a href="priorityLogs/"> available </a></p>"""

target_table = """<h2>Summary of Observed Targets </h2>
<p>For each observation of a field, the targets planned to be observed are specified in an
rsAssignments file. The targets can therefor be matched to field observations. The results of
this matching are here for <a href="obsTargets-{plan}-apo-0.fits">APO</a>
and <a href="obsTargets-{plan}-apo-0.fits">LCO</a>. They can then be sorted by cadence to determine how many
targets where observed total. The plot below shows a summary of total targets per cadence, note the log
scale.</p>

<a href="{plan}-target_summary.pdf"><img src="{plan}-target_summary.png" width="900px/"> </a>

<p>Below the total number of targets observed, as well as the number observed at each observatory,
is shown. A csv file containing this information is also <a href="{plan}-target_summary.txt">available</a>.</p>

<table><tbody>
<tr> <th colspan=2></th> <th colspan=4>robostrategy</th>
<th colspan=3>observesim</th> </tr>

<tr> <th>carton</th> <th>required</th> <th>input</th>
<th>assigned</th> <th>assign_apo</th> <th>assign_lco</th>
<th>total</th> <th>apo</th> <th>lco</th> </tr>
"""

targ_table_row = """<tr><td>{prog}</td> <td>{req}</td> <td>{input}</td>
<td>{assign}</td> <td>{assign_apo}</td> <td>{assign_lco}</td>
<td>{total}</td> <td class="{apo_flag}">{apo}</td> <td class="{lco_flag}">{lco}</td> </tr>"""


agn_metrics = """</tbody></table>
<p>Additional coverage metrics are computed for bhm_spiders_agn </p>

<table><tbody>
<tr> <th> loc </th> <th> area planned (deg<sup>2</sup>) </th> <th> area obs (deg<sup>2</sup>) </th>
<tr> <td> total </td> <td> {total_plan:.2f} </td> <td> {total_obs:.2f} </td> </tr>
<tr> <td> apo </td> <td> {apo_plan:.2f} </td> <td> {apo_obs:.2f} </td> </tr>
<tr> <td> lco </td> <td> {lco_plan:.2f} </td> <td> {lco_obs:.2f} </td> </tr>
</tbody></table>

<a href="{plan}-apo-spiders_v_time.pdf"><img src="{plan}-apo-spiders_v_time.png" width="600px/"> </a>
<a href="{plan}-lco-spiders_v_time.pdf"><img src="{plan}-lco-spiders_v_time.png" width="600px/"> </a>
"""

table_heads = """<h2>Cumulative Exposures</h2>
<p> The plots below show the cumulative exposures for each field cadence class over time. A pdf showing all the cumulative plots is available for <a href="{plan}-apo-cumulative.pdf">APO</a> and <a href="{plan}-lco-cumulative.pdf">LCO</a></p>

<table><tbody>
<tr><td><h4>APO</h4></td> <td><h4>LCO </h4></td></tr>
"""

table_row ="""<tr><td><a href="{apo_png}"><img src="{apo_png}" width="450/"></a></td>
<td><a href="{lco_png}"><img src="{lco_png}" width="450/"></a></td></tr>"""

tail = """</tbody></table>
</body></html>"""


def comp_to_css(assign, obs):
    if assign == 0:
        return "bg-green"
    else:
        frac = obs/assign
    if frac >= 0.9:
        return "bg-green"
    elif frac < 0.9 and frac >= 0.7:
        return "bg-orange"
    else:
        return "bg-red"


def writeWebPage(base, rs_base, plan, version=None):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    html = header + "\n" + bar_plots + "\n" + fits_out + "\n" + target_table + "\n"
    html = html.format(plan=plan)

    targ_sum_file = v_base + plan + "-target_summary.txt"

    targ_sum = np.genfromtxt(targ_sum_file, names=True, delimiter=",",
                             dtype=None, encoding=None)

    for t in targ_sum:
        html += targ_table_row.format(prog=t["carton"], req=t["required"],
                                      total=t["total"], apo=t["apo"], lco=t["lco"],
                                      assign=t["assigned"], assign_apo=t["assign_apo"],
                                      assign_lco=t["assign_lco"], input=t["input"],
                                      apo_flag=comp_to_css(t["assign_apo"], t["apo"]),
                                      lco_flag=comp_to_css(t["assign_lco"], t["lco"]))

    agn_args = dict()
    agn_args["apo_plan"], agn_args["apo_obs"] = spiders_area_for_program(base, rs_base, plan,
                                                                     version=version, loc="apo")
    agn_args["lco_plan"], agn_args["lco_obs"] = spiders_area_for_program(base, rs_base, plan,
                                                                     version=version, loc="lco")

    agn_args["total_plan"] = agn_args["apo_plan"] + agn_args["lco_plan"]
    agn_args["total_obs"] = agn_args["apo_obs"] + agn_args["lco_obs"]
    agn_args["plan"] = plan

    html += agn_metrics.format(**agn_args)

    html += "\n" + table_heads + "\n"

    files = os.listdir(v_base)
    cum_pngs = [f for f in files if "cumulative.png" in f]

    headchars = len(plan) + 5  # 5 for "-lco-" or "-apo"
    tailschars = len("-cumulative.png")

    progs_plotted = list()
    for c in cum_pngs:
        prog = c[headchars:-tailschars]
        if prog not in progs_plotted:
            progs_plotted.append(prog)

    progs_plotted.sort()

    for p in progs_plotted:
        a_file = "{}-apo-{}-cumulative.png".format(plan, p)

        l_file = "{}-lco-{}-cumulative.png".format(plan, p)

        html += table_row.format(apo_png=a_file, lco_png=l_file) + "\n"

    with open(v_base + "/summary.html", "w") as html_file:
        print(html + tail, file=html_file)


def getCounts(res_base, rs_base, plan, version=None, loc="apo"):
    if version is not None:
        v_base = os.path.join(res_base, version)
        v_base += "/"
    else:
        v_base = os.path.join(res_base, plan)
        v_base += "/"
    files = os.listdir(v_base)
    res_files = [f for f in files if "{plan}-{loc}-fields".format(plan=plan, loc=loc) in f]
    counts = list()
    print("reading {} result files for {}".format(len(res_files), v_base))
    for f in res_files:
        sim_data = fitsio.read(v_base + f)
        counts.append(sim_data['nobservations'])

    fields = fitsio.read(rs_base+"{plan}/final/rsAllocationFinal-{plan}-{loc}.fits".format(plan=plan, loc=loc),
                         columns=["fieldid", "needed", "cadence"])

    planned = fields['needed']

    counts = np.mean(np.array(counts), axis=0)

    assert len(sim_data) == len(fields), "THIS IS NOT APPLES AND APPLES!!"

    for s, f in zip(sim_data, fields):
        assert s["fieldid"] == f["fieldid"], "fields mixed up"

    return counts, planned, [c for c in fields["cadence"]]


def countEpochs(res_base, rs_base, plan, version=None, loc="apo"):
    if version is not None:
        v_base = os.path.join(res_base, version)
        v_base += "/"
    else:
        v_base = os.path.join(res_base, plan)
        v_base += "/"
    files = os.listdir(v_base)

    N = 0

    sim_data = fitsio.read(v_base + "{plan}-{loc}-fields-{n}.fits".format(plan=plan, loc=loc, n=N))
    obs_data = fitsio.read(v_base + "{plan}-{loc}-observations-{n}.fits".format(plan=plan, loc=loc, n=N))
    cads = fitsio.read(rs_base + "{plan}/rsCadences-{plan}-{loc}.fits".format(plan=plan, loc=loc))

    counts = list()
    planned = list()
    cadences = list()

    for f in sim_data:
        obs_idx = f["observations"][:int(f["nobservations"])]
        mjds = [obs_data[i]["mjd"] for i in obs_idx]
        cad = f["cadence"]
        this_cad = cads[cads["CADENCE"] == cad][0]
        obs_epochs, last_epoch = epochs_completed(mjds, 0.6)

        counts.append(obs_epochs)
        planned.append(this_cad["NEPOCHS"])
        cadences.append(cad)

    return counts, planned, cadences


def tabulate(counts, planned, cadence):
    completion = dict()
    visits = dict()
    plan_count = dict()

    # print(counts)
    # print(len(counts), len(planned), len(cadence))

    for n, p, c in zip(counts, planned, cadence):
        if c in completion:
            completion[c].append(n/p)
            visits[c].append(n)
            plan_count[c].append(p)
        else:
            completion[c] = [n/p]
            visits[c] = [n]
            plan_count[c] = [p]

    return completion, visits, plan_count


def convertCadence(cad):
    if "_v" in cad:
        cad = cad[:cad.index("_v")]
    split = cad.split("_")
    nums = split[-1]
    name = "".join([str(n) + "_" for n in split[:-1]])
    try:
        epochs, exps = nums.split("x")
    except ValueError:
        return cad
    if int(exps) > 6:
        return cad
    epochs = int(epochs)

    name += "{}x" + exps
    if epochs < 6:
        return name.format("1:5")
    elif epochs < 11:
        return name.format("6:10")
    else:
        return name.format("10+")


def doHist(res_base, rs_base, plan, version=None, loc="apo", level=0.95):
    """Create old histograms by cadence
    """
    args = getCounts(res_base, rs_base, plan, version=version, loc=loc)
    completion, vis_count, plan_count = tabulate(*args)
    args = countEpochs(res_base, rs_base, plan, version=version, loc=loc)
    epoch_completion, epoch_vis_count, epoch_plan_count = tabulate(*args)
    print("found {} cadences".format(len(completion.keys())))

    print(f"N exp {len(vis_count)}, N epoch {len(epoch_vis_count)}")

    orig_keys = completion.keys()
    new_keys = [convertCadence(k) for k in orig_keys]
    sort_keys = np.sort(np.unique(new_keys))
    cad_to_indx = {k: i for i, k in enumerate(sort_keys)}

    reduced_keys = {k: [] for k in sort_keys}
    reduced_vis = {k: [] for k in sort_keys}
    reduced_plan = {k: [] for k in sort_keys}

    for key, v in completion.items():
        k = convertCadence(key)
        reduced_keys[k].extend(v)
        reduced_vis[k].extend(vis_count[key])
        reduced_plan[k].extend(plan_count[key])

    epoch_reduced_keys = {k: [] for k in sort_keys}
    epoch_reduced_vis = {k: [] for k in sort_keys}
    epoch_reduced_plan = {k: [] for k in sort_keys}

    for key, v in epoch_completion.items():
        k = convertCadence(key)
        epoch_reduced_keys[k].extend(v)
        # epoch_reduced_vis[k].extend(epoch_vis_count[key])
        # epoch_reduced_plan[k].extend(epoch_plan_count[key])

    plt.figure(figsize=(12, 12))

    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    x = list()
    mean = list()
    median = list()
    label = list()
    total = list()

    vis = list()
    planned = list()

    emean = list()

    for k, v in reduced_keys.items():
        x.append(cad_to_indx[k])

        mean.append(np.mean(v) * 100)
        median.append(np.median(v) * 100)

        label.append(k.strip())
        total.append(len(v))

        vis_sum = np.sum(reduced_vis[k])
        vis.append(vis_sum)
        plan_sum = np.sum(reduced_plan[k])
        planned.append(plan_sum)

    for k, v in epoch_reduced_keys.items():
        emean.append(np.mean(v) * 100)

    width = 0.3

    ax1.bar(np.array(x) - width, mean, width, color="b", label="exp")
    ax1.bar(np.array(x), emean, width, color="r", label="epoch")
    ax1.set_xticks(x, minor=False)
    ax1.set_xticklabels(label, rotation='vertical')
    ax1.set_ylabel("% comp".format(int(level*100)))
    ax1.set_ylim([0, 105])
    ax1.legend()

    ax3 = ax1.twiny()

    labels = ax1.get_xticks()
    ax3.set_xlim(ax1.get_xlim())
    ax3.set_xticks(labels)
    ax3.set_xticklabels(total, rotation='vertical')

    ax2.bar(np.array(x) - width, vis, width, color="c", label="visits", log=True)
    ax2.bar(np.array(x), planned, width, color="coral", label="planned visits", log=True)
    ax2.set_xticks(x, minor=False)
    ax2.set_xticklabels(label, rotation='vertical')
    ax2.set_ylabel("# visits")
    ax2.legend()
    # ax2.set_yscale("log")

    ax1.set_title("{}: {}".format(plan, loc))

    plt.tight_layout()
    if version is None:
        version = plan
    plt.savefig(os.path.join(res_base, version)+"/{}-{}-cadence_bar.pdf".format(plan, loc))
    plt.savefig(os.path.join(res_base, version)+"/{}-{}-cadence_bar.png".format(plan, loc))


def combineProgramMjds(base, plan, rs_base, version=None, loc="apo", N=0):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    obs_data = fitsio.read(v_base + "obsTargets-{plan}-{loc}-0.fits".format(plan=plan, loc=loc),
                           columns=["obs_mjd", "catalogid", "carton_to_target_pk"])
    comp_data = fitsio.read(rs_base + 
                            "/{plan}/final/rsCompletenessFinal-{plan}-both.fits".format(plan=plan), 
                            columns=["catalogid", "program", "carton", 
                                     f"assigned_{loc}", f"nexps_{loc}",
                                     "carton_to_target_pk"])

    progs = np.unique(comp_data["program"])

    prog_mjds = dict()
    for p in progs:
        comp_prog = comp_data[comp_data["program"] == p]
        w_targs = np.in1d(obs_data['catalogid'], comp_prog['catalogid'])
        prog_mjds[p] = obs_data[w_targs]["obs_mjd"]

    # assigned = comp_data[np.where(comp_data[f"assigned_{loc}"])]

    cartons = np.unique(comp_data["carton"])

    n_exp_total_carton = dict()
    carton_mjds = dict()
    for c in cartons:
        targs = comp_data[np.where(comp_data["carton"] == c)]
        n_exp_total_carton[c] = np.sum(targs[f"nexps_{loc}"])

        w_targs = np.in1d(obs_data['carton_to_target_pk'], targs['carton_to_target_pk'])
        carton_mjds[c] = obs_data[w_targs]["obs_mjd"]

    return prog_mjds, carton_mjds, n_exp_total_carton


def plannedProgExps(plan, rs_base, loc="apo"):
    comp_data = fitsio.read(rs_base + "/{plan}/final/rsCompletenessFinal-{plan}-both.fits".format(plan=plan), 
                            columns=["program", "cadence", "assigned", 
                                     "category", f"assigned_{loc}"])
    cadences = fitsio.read(rs_base + "/{plan}/rsCadences-{plan}-{loc}.fits".format(plan=plan, loc=loc),
                           columns=["CADENCE", "NEXP"])

    comp_data = comp_data[np.where(comp_data[f"assigned_{loc}"])]

    progs = {i: 0 for i in np.unique(comp_data["program"])}
    for c in comp_data:
        if c["category"] != "science":
            continue
        cad = cadences[np.where(cadences["CADENCE"] == c["cadence"])]
        progs[c["program"]] += np.sum(cad["NEXP"])

    return progs


def cumulativePlot(base, plan, rs_base, version=None, loc="apo"):
    prog_mjds, carton_mjds, c_totals = combineProgramMjds(base, plan, rs_base,
                                                          version=version, loc=loc)
    # new_prog = [convertCadence(k) for k in prog_mjds.keys()]
    new_prog = [k for k in prog_mjds.keys()]
    all_mjds = list()
    for k, i in prog_mjds.items():
        all_mjds.extend(list(i))
    min_mjd, max_mjd = int(np.min(all_mjds)), int(np.max(all_mjds))

    sort_keys = np.sort(np.unique(new_prog))

    mjds = np.arange(min_mjd, max_mjd, 1)
    plot_progs = {k: [] for k in sort_keys}

    sort_cartons = np.sort(np.unique([k for k in carton_mjds.keys()]))
    plot_cartons = {k: [] for k in sort_cartons}

    for m in mjds:
        today = dict()
        for k, v in prog_mjds.items():
            count = len(np.where(v < m)[0])
            if k in today:
                today[k].append(count)
            else:
                today[k] = [count]

        for k, v in today.items():
            plot_progs[k].append(np.sum(v))
        
        today = dict()
        for k, v in carton_mjds.items():
            count = len(np.where(v < m)[0])
            if k in today:
                today[k].append(count)
            else:
                today[k] = [count]

        for k, v in today.items():
            plot_cartons[k].append(np.sum(v))

    rows = int(np.ceil(len(np.unique(new_prog))/2))
    fig, axes = plt.subplots(rows, 2, sharex=True)
    fig.set_size_inches(12, rows*1.5)

    prog_to_indx = {k: i for i, k in enumerate(sort_keys)}

    flat_axes = axes.flatten()

    plt.title("{}-{}".format(plan, loc))

    for k, v in plot_progs.items():
        ax = flat_axes[prog_to_indx[k]]
        ax.plot(mjds, v)
        ax.set_title(k)

    for a in axes[:, 0]:
        a.set_ylabel("# visits", fontsize=16)

    axes[-1, 0].set_xlabel("mjd", fontsize=16)
    axes[-1, 1].set_xlabel("mjd", fontsize=16)
    plt.tight_layout()
    if version is None:
        version = plan
    plt.savefig(os.path.join(base, version)+"/{}-{}-cumulative.pdf".format(plan, loc))

    plt.close()

    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    years_fmt = mdates.DateFormatter('%Y')

    programs_to_highlight = ["bhm_rm", "bhm_spiders", "bhm_csc",
                             "mwm_galactic", "mwm_rv", "mwm_planet"]

    year1 = list()
    year2 = list()
    year3 = list()
    year4 = list()
    year5 = list()
    year5 = list()

    used_progs = list()  # ugh why are the observatories not the same...

    byYearCsv = ("carton, planned, year1_N, year1_frac, year2_N, year2_frac, "
                 "year3_N, year3_frac, year4_N, year4_frac, year5_N, year5_frac\n")

    plannedCounts = plannedProgExps(plan, rs_base, loc=loc)
    for k, v in plot_progs.items():
        if "ops" in k:
            continue
        elif not k in plannedCounts:
            print(f"{k} not in plannedCounts!!")
            continue
        plt.figure(figsize=(8,5))
        ax = plt.subplot(111)
        cal_days = daycnv(mjds+2400000.5, mode="dt")
        ax.plot(cal_days, v)
        ax.set_title(k)
        ax.set_ylabel("# visits", fontsize=16)
        ax.set_xlabel("date", fontsize=16)
        ax.axhline(plannedCounts[k])

        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_major_formatter(years_fmt)
        ax.xaxis.set_minor_locator(months)

        if version is None:
            version = plan
        plt.savefig(os.path.join(base, version)+"/{}-{}-{}-cumulative.png".format(plan, loc, k))
        plt.close()

        if k in programs_to_highlight:
            used_progs.append(k)
            bars = yearBars(mjds, v)
            year1.append(bars[0])
            year2.append(bars[1])
            year3.append(bars[2])
            year4.append(bars[3])
            year5.append(bars[4])

    for k, v in plot_cartons.items():
        if "ops" in k:
            continue

        bars = yearBars(mjds, v)
        byYearCsv += (f"{k}, {c_totals[k]}, {bars[0]}, {bars[0]/c_totals[k]:.2f}, "
                      f"{bars[1]}, {bars[1]/c_totals[k]:.2f}, "
                      f"{bars[2]}, {bars[2]/c_totals[k]:.2f}, "
                      f"{bars[3]}, {bars[3]/c_totals[k]:.2f}, "
                      f"{bars[4]}, {bars[4]/c_totals[k]:.2f}\n")

    csvPath = os.path.join(base, version) +\
              "/{}-{}-cartonsByYear.csv".format(plan, loc)
    
    with open(csvPath, "w") as csvFile:
        print(byYearCsv, file=csvFile)

    cartons = np.genfromtxt(csvPath, names=True,
                            dtype=None, encoding='utf-8', delimiter=',')

    figBase = os.path.join(base, version)

    for prog in ["bhm", "mwm"]:
        sub = cartons[np.where([prog in l for l in cartons["carton"]])]

        plt.figure(figsize=(9, 10))
        ax1 = plt.subplot(111)
        width = 0.1
        x_cor = np.arange(0, len(sub))

        ax1.barh(x_cor - 0.2, np.clip(sub["year1_frac"], 0, 2), width, color="c", label="year 1")
        ax1.barh(x_cor - 0.1, np.clip(sub["year2_frac"], 0, 2), width, color="b", label="year 2")
        ax1.barh(x_cor + 0.0, np.clip(sub["year3_frac"], 0, 2), width, color="m", label="year 3")
        ax1.barh(x_cor + 0.1, np.clip(sub["year4_frac"], 0, 2), width, color="r", label="year 4")
        rects = ax1.barh(x_cor + 0.2, np.clip(sub["year5_frac"], 0, 2), width, color="y", label="year 5")
        ax1.bar_label(rects, labels=[c["planned"] for c in sub], padding=1)

        ax1.set_yticks(x_cor, minor=False)
        ax1.set_yticklabels(sub["carton"])
        ax1.set_xlabel("Observed/Planned")
        ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        ax1.tick_params(labeltop=True)

        ax1.set_title("Fractional Carton Completion")

        ax1.set_xlim(0, 2.3)

        ax1.axvline(1, linewidth=1, linestyle="--", color="k", alpha=0.5)

        plt.tight_layout()
        plt.savefig(os.path.join(base, version)+f"/{plan}-{loc}-allYearCartonFrac_{prog}.png")
        plt.savefig(os.path.join(base, version)+f"/{plan}-{loc}-allYearCartonFrac_{prog}.pdf")
        plt.close()

    # make bar plots of programs by year
    plt.figure(figsize=(8, 5))
    ax1 = plt.subplot(111)
    width = 0.1
    x_cor = np.arange(0, len(used_progs))

    ax1.bar(x_cor - 0.2, year1, width, color="c", label="year 1", log=True)
    ax1.bar(x_cor - 0.1, year2, width, color="b", label="year 2", log=True)
    ax1.bar(x_cor + 0.0, year3, width, color="m", label="year 3", log=True)
    ax1.bar(x_cor + 0.1, year4, width, color="r", label="year 4", log=True)
    ax1.bar(x_cor + 0.2, year5, width, color="y", label="year 5", log=True)

    ax1.set_xticks(x_cor, minor=False)
    ax1.set_xticklabels(used_progs, rotation='vertical')
    ax1.set_ylabel("visits")
    ax1.legend()

    plt.tight_layout()

    plt.savefig(os.path.join(base, version)+"/{}-{}-by_year.png".format(plan, loc))
    plt.close()   


def yearBars(mjds, v):
    jan_2022 = 59580
    cuts = jan_2022 + np.arange(365, 365.25*6, 365.25)
    perYear = list()

    for c in cuts:
        use = np.where(mjds <= c)
        doneAsOf = np.array(v)[use]
        if len(doneAsOf) > 0:
            perYear.append(np.max(doneAsOf))
        else:
            perYear.append(0)
    return perYear


def passesCadence(targ_mjds, cadence=None):
    for m in targ_mjds:
        if m > 5e4:
            # dumb for now, get back to this
            return True


def countPlanned(program, programs, gots):
    """since utah has an old version of fitsio, strings are a pain.
    assuming things are the same length and maintain their order,
    this works. those assumptions should be good because python.
    """
    assert len(programs) == len(gots)
    where_prog = np.where(programs == program)
    prog_targs = gots[where_prog]
    got = np.where(prog_targs == 1)

    return len(prog_targs), len(got[0])


def plotTargMetric(base, rs_base, plan, version=None, reqs_file=None):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    apo_targs = fitsio.read(v_base + "obsTargets-{plan}-{loc}-0.fits".format(plan=plan, loc="apo"),
                            columns=["program", "catalogid"])
    lco_targs = fitsio.read(v_base + "obsTargets-{plan}-{loc}-0.fits".format(plan=plan, loc="lco"),
                            columns=["program", "catalogid"])

    all_comp = fitsio.read(rs_base + "/{plan}/final/rsCompletenessFinal-{plan}-both.fits".format(plan=plan),
                           columns=["catalogid", "program", "assigned_lco", "assigned_apo"])

    apo_comp = all_comp[np.where(all_comp["assigned_apo"])]
    apo_comp_prog = np.array([p.strip() for p in apo_comp["program"]])

    lco_comp = all_comp[np.where(all_comp["assigned_lco"])]
    lco_comp_prog = np.array([p.strip() for p in lco_comp["program"]])

    programs = [p for p in np.unique(apo_comp_prog) if p.lower() != "sky"]
    apo_done_counts = list()
    lco_done_counts = list()
    for p in programs:
        # count observesim done
        apo_prog = apo_comp[apo_comp["program"] == p]
        lco_prog = apo_comp[apo_comp["program"] == p]
        apo_done_counts.append(np.sum(np.in1d(apo_prog['catalogid'], apo_targs['catalogid'])))
        lco_done_counts.append(np.sum(np.in1d(lco_prog['catalogid'], lco_targs['catalogid'])))

    plt.figure(figsize=(16, 10))

    ax1 = plt.subplot(111)

    x = np.arange(len(programs))

    width = 0.3

    ax1.bar(x - width, lco_done_counts, width, color="c", label="lco", log=True)
    ax1.bar(x, apo_done_counts, width, color="r", label="apo", log=True)
    ax1.set_xticks(x, minor=False)
    ax1.set_xticklabels(programs, rotation='vertical')
    ax1.set_ylabel("# targs")
    # ax1.set_yscale("log")

    ax3 = ax1.twiny()

    labels = ax1.get_xticks()
    ax3.set_xlim(ax1.get_xlim())
    ax3.set_xticks(labels)
    ax3.set_xticklabels([x + y for x, y in zip(apo_done_counts, lco_done_counts)], rotation='vertical')

    ax1.legend()

    if reqs_file is None:
        reqs_file = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + "/etc/cadence_reqs.yml"

    req_by_cad = yaml.load(open(reqs_file))

    sum_text = ("{cad:30s}, {req:9s}, {input:8s}, {assign:9s}, {total:8s}, "
                "{apo:8s}, {lco:8s}, {assign_apo:10s}, {assign_lco:10s}\n").format(
                cad="carton", req="required", input="input", total="total",
                apo="apo", lco="lco", assign_apo="assign_apo",
                assign_lco="assign_lco", assign="assigned")

    for i, k in enumerate(programs):
        # check robostrategy planned vs assigned
        apo = apo_done_counts[i]
        apo_plan, apo_assign = countPlanned(k.strip(), apo_comp_prog, apo_comp["assigned_apo"])

        lco = lco_done_counts[i]
        lco_plan, lco_assign = countPlanned(k.strip(), lco_comp_prog, lco_comp["assigned_lco"])

        if k.strip() in req_by_cad:
            req = req_by_cad[k.strip()]
        else:
            req = ""

        sum_text += ("{cad:30s}, {req:9s}, {input:8d}, {assign:9d}, {total:8d}, "
                     "{apo:8d}, {lco:8d}, {assign_apo:10d}, {assign_lco:10d}\n").format(
                     cad=k.strip(), req=str(req), input=apo_plan+lco_plan,
                     assign=apo_assign+lco_assign, total=apo+lco, apo=apo,
                     lco=lco, assign_apo=apo_assign, assign_lco=lco_assign)

    res_base = v_base + plan

    plt.savefig(res_base + "-target_summary.pdf")
    plt.savefig(res_base + "-target_summary.png")
    with open(res_base + "-target_summary.txt", "w") as sum_file:
        print(sum_text, file=sum_file)


def compute_area_above_threshold(targets, obs_targets, threshold, nside, loc="apo"):
    '''
     Computes the sky area above a given completeness threshold
     Counts sky pixels on a HEALPix gris with NSIDE=nside

     inputs:
       targets - numpy recarray-like object with columns 'ra', 'dec', 'got'
                 corresponding to planned targets

       obs_targets - numpy recarray-like object with columns 'ra' and 'dec',
                     corresponding to actually observed targets

       threshold - floatingpoint number in range [0,1]

       nside - HEALPix parameter controlling pixel size

    '''

    assert threshold >= 0.0 and threshold <= 1.0
    assert len(targets) >= 1
    assert (nside >= 1) and (nside < 65536)

    npix = hp.nside2npix(nside)
    pixarea = hp.nside2pixarea(nside, degrees=True)

    hpx = hp.ang2pix(nside, targets["ra"], targets["dec"], lonlat=True)
    all_map = np.bincount(hpx, minlength=npix)
    got_map = np.bincount(hpx, weights=np.where(targets[f"assigned_{loc}"] > 0, 1, 0), minlength=npix)

    obs_hpx = hp.ang2pix(nside, obs_targets["ra"], obs_targets["dec"], lonlat=True)
    obs_map = np.bincount(obs_hpx, minlength=npix)

    with np.errstate(divide='ignore', invalid='ignore'):
        frac_map = np.divide(got_map, all_map)
        map_completed = np.where(frac_map >= threshold, 1, 0)

        obs_frac_map = np.divide(obs_map, all_map)
        obs_completed = np.where(obs_frac_map >= threshold, 1, 0)

    planned_area_completed = pixarea * np.sum(map_completed)
    obs_area_completed = pixarea * np.sum(obs_completed)

    return planned_area_completed, obs_area_completed


def grab_summary_files(base, rs_base, plan, version=None, loc="apo"):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    obs_file = v_base + "obsTargets-{plan}-{loc}-0.fits".format(plan=plan, loc=loc)
    comp_file = rs_base + "/{plan}/final/rsCompletenessFinal-{plan}-both.fits".format(plan=plan)

    all_targets = fitsio.read(comp_file,
                              columns=["catalogid", "covered", "program", "ra", "dec", f"assigned_{loc}"])

    loc_targs = np.extract(all_targets[f"assigned_{loc}"] > 0, all_targets)

    progs = np.unique(loc_targs["program"])

    # since mike is naming these with version numbers now
    # need to find the right version
    spiders_names = [p for p in progs if "spiders" in p]
    assert len(spiders_names) != 0, "didn't find an appropriate spiders_agn carton!"
    # prog_name = spiders_names[0]

    # targets = np.extract(loc_targs['carton'] == prog_name, loc_targs)

    w_spiders = [l in spiders_names for l in loc_targs['program']]
    targets = np.extract(w_spiders, loc_targs)

    obs_targets = fitsio.read(obs_file,
                              columns=["catalogid", "program", "obs_mjd", "ra", "dec"])

    # obs_targets = np.extract(obs_targets['carton'] == prog_name, obs_targets)

    # w_spiders = [l in spiders_names for l in obs_targets['program']]
    # obs_targets = np.extract(w_spiders, obs_targets)

    mask = np.in1d(obs_targets['catalogid'], targets['catalogid'])
    obs_targets = obs_targets[mask]

    # for time domain, ensure we get the earliest version
    obs_targets = np.sort(obs_targets, order="obs_mjd")

    pks, idxs = np.unique(obs_targets["catalogid"], return_index=True)

    # since obs_targets has an entry for each observation!!
    obs_targets = obs_targets[idxs]

    # we'll just leave this debug hint...
    # print(f"{loc}, {len(obs_targets)}, {len(targets)}")

    return targets, obs_targets


def spiders_area_for_program(base, rs_base, plan, version=None, loc="apo"):
    targets, obs_targets = grab_summary_files(base, rs_base, plan, version=version, loc=loc)

    planned, obs = compute_area_above_threshold(targets, obs_targets, threshold=0.8, nside=64, loc=loc)

    print(f"{plan} planned: Ntargets={len(targets)}, Area={planned:.2f} deg^2")
    print(f"{plan} obs: Ntargets={len(obs_targets)}, Area={obs:.2f} deg^2")

    return planned, obs


def spiders_area_for_program_time(base, rs_base, plan, version=None, loc="apo"):
    targets, obs_targets = grab_summary_files(base, rs_base, plan, version=version, loc=loc)

    start_mjd = int(np.min(obs_targets["obs_mjd"]))
    end_mjd = int(np.max(obs_targets["obs_mjd"])) + 1

    # force it to be multiples of 10 so I can print progress
    mjds = np.arange(start_mjd//10*10, end_mjd//10*10, 10)
    area_comp = list()

    threshold = 0.8
    nside = 64
    assert threshold >= 0.0 and threshold <= 1.0
    assert len(targets) >= 1
    assert (nside >= 1) and (nside < 65536)

    npix = hp.nside2npix(nside)
    pixarea = hp.nside2pixarea(nside, degrees=True)

    hpx = hp.ang2pix(nside, targets["ra"], targets["dec"], lonlat=True)
    all_map = np.bincount(hpx, minlength=npix)

    start = time.time()
    for m in mjds:
        mjd_targs = np.extract(obs_targets['obs_mjd'] < m, obs_targets)
        obs_hpx = hp.ang2pix(nside, mjd_targs["ra"], mjd_targs["dec"], lonlat=True)
        obs_map = np.bincount(obs_hpx, minlength=npix)

        with np.errstate(divide='ignore', invalid='ignore'):
            obs_frac_map = np.divide(obs_map, all_map)
            obs_completed = np.where(obs_frac_map >= threshold, 1, 0)

        obs_area_completed = pixarea * np.sum(obs_completed)
        # print(m, len(mjd_targs), obs_area_completed)

        area_comp.append(obs_area_completed)

    plt.figure(figsize=(16, 10))
    ax = plt.subplot(111)

    ax.plot(mjds, area_comp)
    ax.set_xlabel("MJD")
    ax.set_ylabel(r"area (deg$^2$)")
    ax.set_title(f"SPIDERS sky coverage vs time, {loc}")

    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    res_base = v_base + plan

    plt.savefig(res_base + f"-{loc}-spiders_v_time.png")
    plt.savefig(res_base + f"-{loc}-spiders_v_time.pdf")
    plt.close()


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