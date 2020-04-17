import time, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import fitsio

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
        
#         if real_fid % 100 == 0:
#             print(real_fid, (time.time() - start_time)/60.)
        
        ids, cadences, targ_mjds = read_field(real_fid, mjds, assign)
        all_targs.extend(ids)
        all_cads.extend(cadences)
        all_mjds.extend(targ_mjds)
        all_fields.extend([real_fid for f in ids])
    
    assert len(all_targs) == len(all_cads), "targ != cad!!"
    assert len(all_targs) == len(all_mjds), "targ != mjd!!"
    assert len(all_targs) == len(all_fields), "targ != field!!"
    
    dtype = [('pk', np.int32),
             ('cadence', np.dtype('a40')),
             ('field_id', np.int32),
             ('obs_mjd', np.float64)]
    obs_targs = np.zeros(len(all_targs), dtype=dtype)
    
    obs_targs["pk"] = all_targs
    obs_targs["cadence"] = all_cads
    obs_targs["field_id"] = all_fields
    obs_targs["obs_mjd"] = all_mjds
    
    if save:
        fitsio.write(v_base + "obsTargets-{plan}-{loc}.fits".format(plan=plan, loc=loc), obs_targs,
                 clobber=True)

    return obs_targs


if __name__ == "__main__":
    usage = "sdss5_simulate"
    description = "Simulate the SDSS-V schedule"
    parser = argparse.ArgumentParser(description=description, usage=usage)
    parser.add_argument("-b", "--base", dest="base", type=str,
                        required=True, help="output FITS base name")
    parser.add_argument("-r", "--rsbase", dest="rs_base", type=str,
                        required=False, help="location of rs files if not specified in env",
                        default=None)
    parser.add_argument("-p", "--plan", dest="plan", type=str,
                        required=False, help="design plan",
                        default='plan-0')
    parser.add_argument("-v", "--version", dest="version", type=str,
                        required=False, help="use versioned directory for output",
                        default=None)
    args = parser.parse_args()
    base = args.base
    plan = args.plan
    version = args.version
    rs_base = args.rs_base

    if rs_base is None:
        rs_base = os.getenv('OBSERVING_PLAN_DIR') + "/"

    countFields(base, rs_base, plan, version=version, loc="apo", N=0)
    countFields(base, rs_base, plan, version=version, loc="lco", N=0)
