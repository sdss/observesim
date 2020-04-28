import time
import os
import numpy as np
import matplotlib.pyplot as plt
import fitsio
import yaml


def read_field(field_id, exp_to_mjd, assign):
    # fetch info, match mjd to exp
    field_targs = assign[np.where(assign["fieldid"] == field_id)]

    if len(np.unique(field_targs["exposure"])) > len(exp_to_mjd):
        exp_to_mjd.extend([-1 for i in range(len(np.unique(field_targs["exposure"])) - len(exp_to_mjd))])

    mjds = [exp_to_mjd[i] for i in field_targs["exposure"]]
    return field_targs["pk"], field_targs["cadence"], mjds


def countFields(res_base, rs_base, plan, version=None, loc="apo", N=0, save=True):
    if version is not None:
        v_base = os.path.join(res_base, version)
        v_base += "/"
    else:
        v_base = os.path.join(res_base, plan)
        v_base += "/"

    assign = fitsio.read(rs_base + "{plan}/rsAssignments-{plan}-{loc}.fits".format(plan=plan, loc=loc))
    assign = assign[np.where(assign["pk"] != -1)]

    completeness = fitsio.read(rs_base + "gamma-0/" + "rsCompleteness-gamma-0-apo.fits")
    completeness = completeness[np.where(completeness["pk"] != -1)]

    comp_dict = dict()
    for c in completeness:
        comp_dict[c["pk"]] = (c["targetid"], c["program"], c["got"], c["cadence"])

    allocation = fitsio.read(rs_base + "{plan}/rsAllocation-{plan}-{loc}.fits".format(plan=plan, loc=loc))

    sim_data = fitsio.read(v_base + "{plan}-{loc}-fields-{n}.fits".format(plan=plan, loc=loc, n=N))
    obs_data = fitsio.read(v_base + "{plan}-{loc}-observations-{n}.fits".format(plan=plan, loc=loc, n=N))

    # prep out struct
    all_targs = list()
    all_cads = list()
    all_mjds = list()
    all_fields = list()

    start_time = time.time()

    for f in sim_data:
        obs_idx = f["observations"][:int(f["nobservations"])]
        mjds = [obs_data[i]["mjd"] for i in obs_idx]
        cad = f["cadence"]
        real_fid = allocation[f["fieldid"]]["fieldid"]

        ids, cadences, targ_mjds = read_field(real_fid, mjds, assign)
        all_targs.extend(ids)
        all_cads.extend(cadences)
        all_mjds.extend(targ_mjds)
        all_fields.extend([real_fid for f in ids])

    assert len(all_targs) == len(all_cads), "targ != cad!!"
    assert len(all_targs) == len(all_mjds), "targ != mjd!!"
    assert len(all_targs) == len(all_fields), "targ != field!!"

    targ_ids = list()
    programs = list()
    gots = list()

    for k, c in zip(all_targs, all_cads):
        targetid, program, got, cadence = comp_dict[k]
        targ_ids.append(targetid)
        programs.append(program)
        gots.append(got)
        if cadence != c:
            print("field cad: {}, targ cad: {}".format(c, cadence))
            print("pk: {}, targ id: {}".format(k, targetid))
    dtype = [('pk', np.int32),
             ('target_id', np.int32),
             ('cadence', np.dtype('a40')),
             ('program', np.dtype('a40')),
             ('field_id', np.int32),
             ('got', np.int32),
             ('obs_mjd', np.float64)]
    obs_targs = np.zeros(len(all_targs), dtype=dtype)

    obs_targs["pk"] = all_targs
    obs_targs["cadence"] = all_cads
    obs_targs["field_id"] = all_fields
    obs_targs["obs_mjd"] = all_mjds
    obs_targs["got"] = gots
    obs_targs["program"] = programs
    obs_targs["target_id"] = targ_ids

    if save:
        fitsio.write(v_base + "obsTargets-{plan}-{loc}.fits".format(plan=plan, loc=loc), obs_targs,
                 clobber=True)

    return obs_targs

header = """<html><head><meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
<style>
    table, th, tr, td {{border: 1px solid black}}
</style>
</head><body>

<h1>Observesim: {plan}</h1>
<p>This page summarizes the observesim results for the {plan} robostrategy runs <a href="https://data.sdss.org/sas/sdss5/sandbox/robostrategy/allocations/{plan}">found here</a></p>

<h2>Summary</h2>
<p> The results for each observatory are summarized below. The average completion % for each field cadence class is shown, as well as the number of exposures falling into each class (note the log scale). </p>"""

bar_plots = """<a href="{plan}-apo-cadence_bar.pdf"><img src="{plan}-apo-cadence_bar.png" width="600px/"> </a>
<a href="{plan}-lco-cadence_bar.pdf"><img src="{plan}-lco-cadence_bar.png" width="600px/"> </a>"""

fits_out = """<h2> Output files: </h2>
<p> Each input field from the rsAllocation file is recorded in a fields file for  <a href="{plan}-apo-fields-0.fits">APO</a> and <a href="{plan}-lco-fields-0.fits">LCO</a>.
The fields files contain indices into an observations to record info for each observation. <a href="{plan}-apo-observations-0.fits"> APO observations</a>, and <a href="{plan}-lco-observations-0.fits"> LCO observations</a>. </p>"""

target_table = """<h2>Summary of Observed Targets </h2>
<p>For each observation of a field, the targets planned to be observed are specified in an
rsAssignments file. The targets can therefor be matched to field observations. The results of
this matching are here for <a href="obsTargets-{plan}-apo.fits">APO</a>
and <a href="obsTargets-{plan}-apo.fits">LCO</a>. They can then be sorted by cadence to determine how many
targets where observed total. The plot below shows a summary of total targets per cadence, note the log
scale.</p>

<a href="{plan}-target_summary.pdf"><img src="{plan}-target_summary.png" width="900px/"> </a>

<p>Below the total number of targets observed, as well as the number observed at each observatory,
is shown. A csv file containing this information is also <a href="{plan}-target_summary.txt">available</a>.</p>

<table><tbody>
<tr> <th colspan=2></th> <th colspan=4>robostrategy</th>
<th colspan=3>observesim</th> </tr>

<tr> <th>program</th> <th>required</th> <th>input</th>
<th>assigned</th> <th>assign_apo</th> <th>assign_lco</th>
<th>total</th> <th>apo</th> <th>lco</th> </tr>
"""

targ_table_row = """<tr><td>{prog}</td> <td>{req}</td> <td>{input}</td>
<td>{assign}</td> <td>{assign_apo}</td> <td>{assign_lco}</td>
<td>{total}</td> <td>{apo}</td> <td>{lco}</td> </tr>"""

table_heads = """</tbody></table>
<h2>Cumulative Exposures</h2>
<p> The plots below show the cumulative exposures for each field cadence class over time. A pdf showing all the cumulative plots is available for <a href="{plan}-apo-cumulative.pdf">APO</a> and <a href="{plan}-lco-cumulative.pdf">LCO</a></p>

<table><tbody>
<tr><td><h4>APO</h4></td> <td><h4>LCO </h4></td></tr>"""

table_row ="""<tr><td><a href="{apo_png}"><img src="{apo_png}" width="450/"></a></td>
<td><a href="{lco_png}"><img src="{lco_png}" width="450/"></a></td></tr>"""

tail = """</tbody></table>
</body></html>"""


def writeWebPage(base, plan, version=None):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    html = header + "\n" + bar_plots + "\n" + fits_out + "\n" + target_table + "\n"
    html = html.format(plan=plan)

    targ_sum_file = v_base + plan + "-target_summary.txt"

    targ_sum = np.genfromtxt(targ_sum_file, names=True, delimiter=",", dtype=None, encoding=None)

    for t in targ_sum:
        html += targ_table_row.format(prog=t["program"], req=t["required"],
                                total=t["total"], apo=t["apo"], lco=t["lco"],
                                assign=t["assigned"], assign_apo=t["assign_apo"],
                                assign_lco=t["assign_lco"], input=t["input"])

    html += "\n" + table_heads + "\n"

    files = os.listdir(v_base)
    cum_pngs = [f for f in files if "cumulative.png" in f]
    lco = [f for f in cum_pngs if "lco" in f and "none" not in f]
    apo = [f for f in cum_pngs if "apo" in f and "none" not in f]
    lco.sort()
    apo.sort()

    len_l = len(lco)
    len_a = len(apo)

    for i in range(max(len_l, len_a)):
        if i <= len_l - 1:
            l_file = lco[i]
        else:
            l_file = ""
        if i <= len_a - 1:
            a_file = apo[i]
        else:
            a_file = ""
        html += table_row.format(apo_png=a_file, lco_png=l_file) + "\n"

    return html + tail


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

    fields = fitsio.read(rs_base+"{plan}/rsAllocation-{plan}-{loc}.fits".format(plan=plan, loc=loc))

    planned = fields['needed']

    counts = np.mean(np.array(counts), axis=0)

    return counts, planned, [c.decode() for c in fields["cadence"]]


def tabulate(counts, planned, cadence):
    completion = dict()
    visits = dict()
    plan_count = dict()

    print(counts)
    print(len(counts), len(planned), len(cadence))

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
    if not "-" in cad:
        return cad.strip()
    name, num = cad.split("-")
    num = int(num)
    if num < 6:
        return name + "-1:5"
    elif num < 11:
        return name + "-6:10"
    else:
        return name + "-10+"


def doHist(res_base, rs_base, plan, version=None, loc="apo", level=0.95):
    args = getCounts(res_base, rs_base, plan, version=version, loc=loc)
    completion, vis_count, plan_count = tabulate(*args)
    print("found {} cadences".format(len(completion.keys())))

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

    width = 0.3

    ax1.bar(np.array(x) - width, mean, width, color="b", label="mean")
    ax1.bar(np.array(x), median, width, color="r", label="median")
    ax1.set_xticks(x, minor=False)
    ax1.set_xticklabels(label, rotation='vertical')
    ax1.set_ylabel("% comp".format(int(level*100)))
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


def getEpochsFromFile(base, plan, version=None, loc="apo", N=0):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    sim_data = fitsio.read(v_base + "{plan}-{loc}-fields-{n}.fits".format(plan=plan, loc=loc, n=N))
    obs_data = fitsio.read(v_base + "{plan}-{loc}-observations-{n}.fits".format(plan=plan, loc=loc, n=N))
    cad_mjds = dict()
    for f in sim_data:
        obs_idx = f["observations"][:int(f["nobservations"])]
        mjds = [obs_data[i]["mjd"] for i in obs_idx]
        cad = f["cadence"].decode()
        if cad in cad_mjds:
            cad_mjds[cad].extend(mjds)
        else:
            cad_mjds[cad] = mjds

    cad_mjds = {k: np.array(v) for k, v in cad_mjds.items()}
    return cad_mjds


def cumulativePlot(base, plan, version=None, loc="apo"):
    cad_mjds = getEpochsFromFile(base, plan, version=version, loc=loc)
    new_cad = [convertCadence(k) for k in cad_mjds.keys()]
    all_mjds = list()
    for k, i in cad_mjds.items():
        all_mjds.extend(list(i))
    min_mjd, max_mjd = int(np.min(all_mjds)), int(np.max(all_mjds))

    sort_keys = np.sort(np.unique(new_cad))

    mjds = np.arange(min_mjd, max_mjd, 1)
    plot_cads = {k: [] for k in sort_keys}

    for m in mjds:
        today = dict()
        for k, v in cad_mjds.items():
            count = len(np.where(v < m)[0])
            if convertCadence(k) in today:
                today[convertCadence(k)].append(count)
            else:
                today[convertCadence(k)] = [count]

        for k, v in today.items():
            plot_cads[k].append(np.sum(v))

    rows = int(np.ceil(len(np.unique(new_cad))/2))
    fig, axes = plt.subplots(rows, 2, sharex=True)
    fig.set_size_inches(12, rows*1.5)

    cad_to_indx = {k: i for i, k in enumerate(sort_keys)}

    flat_axes = axes.flatten()

    plt.title("{}-{}".format(plan, loc))

    for k, v in plot_cads.items():
        ax = flat_axes[cad_to_indx[k]]
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

    for k, v in plot_cads.items():
        plt.figure(figsize=(8,5))
        ax = plt.subplot(111)
        ax.plot(mjds, v)
        ax.set_title(k)
        ax.set_ylabel("# visits", fontsize=16)
        ax.set_xlabel("mjd", fontsize=16)
        if version is None:
            version = plan
        plt.savefig(os.path.join(base, version)+"/{}-{}-{}-cumulative.png".format(plan, loc, k))
        plt.close()


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

    apo_targs = fitsio.read(v_base + "obsTargets-{plan}-{loc}.fits".format(plan=plan, loc="apo"))
    lco_targs = fitsio.read(v_base + "obsTargets-{plan}-{loc}.fits".format(plan=plan, loc="lco"))
    lco_targs_mjds = {i: [] for i in np.unique(lco_targs["pk"])}
    apo_targs_mjds = {i: [] for i in np.unique(apo_targs["pk"])}

    targ_to_prog = dict()

    for t in lco_targs:
        lco_targs_mjds[t["pk"]].append(t["obs_mjd"])
        targ_to_prog[t["pk"]] = t["program"].decode()

    for t in apo_targs:
        apo_targs_mjds[t["pk"]].append(t["obs_mjd"])
        targ_to_prog[t["pk"]] = t["program"].decode()

    # #####################
    # #####################

    all_progs = np.sort(np.unique([v for k, v in targ_to_prog.items()]))
    lco_progs = {c: [] for c in all_progs}
    apo_progs = {c: [] for c in all_progs}

    for t, v in lco_targs_mjds.items():
        if passesCadence(v):
            lco_progs[targ_to_prog[t]].append(t)

    for t, v in apo_targs_mjds.items():
        if passesCadence(v):
            apo_progs[targ_to_prog[t]].append(t)

    # #####################
    # #####################

    plt.figure(figsize=(16, 10))

    ax1 = plt.subplot(111)

    apo_counts = [len(v) for k, v in apo_progs.items()]
    lco_counts = [len(v) for k, v in lco_progs.items()]
    names = [k.strip() for k in apo_progs.keys()] # same as lco_prog.keys()
    x = range(len(names))

    ax1.bar(x, lco_counts, color="c", label="lco", log=True)
    ax1.bar(x, apo_counts, bottom=lco_counts, color="r", label="apo", log=True)
    ax1.set_xticks(x, minor=False)
    ax1.set_xticklabels(names, rotation='vertical')
    ax1.set_ylabel("# targs")
    # ax1.set_yscale("log")

    ax3 = ax1.twiny()

    labels = ax1.get_xticks()
    ax3.set_xlim(ax1.get_xlim())
    ax3.set_xticks(labels)
    ax3.set_xticklabels([x + y for x, y in zip(apo_counts, lco_counts)], rotation='vertical')

    ax1.legend()

    if reqs_file is None:
        reqs_file = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + "/etc/cadence_reqs.yml"

    req_by_cad = yaml.load(open(reqs_file))

    apo_comp = fitsio.read(rs_base + "/{plan}/rsCompleteness-{plan}-{loc}.fits".format(plan=plan, loc="apo"),
                           columns=["pk", "program", "got"])
    apo_comp_prog = np.array([p.decode().strip() for p in apo_comp["program"]])

    lco_comp = fitsio.read(rs_base + "/{plan}/rsCompleteness-{plan}-{loc}.fits".format(plan=plan, loc="lco"),
                           columns=["pk", "program", "got"])
    lco_comp_prog = np.array([p.decode().strip() for p in lco_comp["program"]])

    sum_text = ("{cad:30s}, {req:9s}, {input:8s}, {assign:9s}, {total:8s}, "
                "{apo:8s}, {lco:8s}, {assign_apo:10s}, {assign_lco:10s}\n").format(
                cad="program", req="required", input="input", total="total",
                apo="apo", lco="lco", assign_apo="assign_apo",
                assign_lco="assign_lco", assign="assigned")

    for k in apo_progs.keys():
        apo = len(apo_progs[k])
        apo_plan, apo_assign = countPlanned(k.strip(), apo_comp_prog, apo_comp["got"])

        lco = len(lco_progs[k])
        lco_plan, lco_assign = countPlanned(k.strip(), lco_comp_prog, lco_comp["got"])

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
