import os
import numpy as np
import fitsio


def countCartons(v_base, plan, rs_base):
    obs_data_apo = fitsio.read(v_base + f"obsTargets-{plan}-apo-0.fits",
                           columns=["carton"])
    obs_data_lco = fitsio.read(v_base + f"obsTargets-{plan}-lco-0.fits",
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

def fieldCounts(v_base, plan, rs_base, loc="apo"):
    fields = fitsio.read(v_base + f"{plan}-{loc}-fields-0.fits")

    tabulated = {
        "dark_10x4_4yr_v1": [0, 0],
        "dark_174x8_v1": [0, 0],
        "dark_100x8_v1": [0, 0],
        "dark_2x1_v1": [0, 0],
        "dark_2x2_v1": [0, 0],
        "dark_2x4_v1": [0, 0],
        "bright_x1": [0, 0],
        "bright_x2": [0, 0],
        "bright_x4": [0, 0]
    }

    for f in fields:
        cad = f["cadence"]
        if not cad in tabulated.keys():
            if "x1" in cad:
                cad = "bright_x1"
            elif "x2" in cad:
                cad = "bright_x2"
            elif "x4" in cad:
                cad = "bright_x4"
            else:
                raise Exception(f"missing cadence! {cad}")

        tabulated[cad][0] += f["nobservations"]
    
    allocation = fitsio.read(rs_base + 
                 f"/{plan}/final/rsAllocationFinal-{plan}-{loc}.fits")
     
    for f in allocation:
        cad = f["cadence"]
        if not cad in tabulated.keys():
            if "x1" in cad:
                cad = "bright_x1"
            elif "x2" in cad:
                cad = "bright_x2"
            elif "x4" in cad:
                cad = "bright_x4"
            else:
                raise Exception(f"missing cadence!  {cad}")

           
        tabulated[cad][1] += f["nfilled"]
    
    return tabulated

def quickSummary(base, plan, rs_base, version=None):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"
    header = """<html><head><meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
    <style>
        table, th, tr, td {{border: 2px solid black}}
    </style>
    </head><body>

    <table><tbody>

    <tr> <th>carton</th> <th>RS planned APO targs</th> <th>Done APO targs</th>
    <th>RS planned LCO targs</th> <th>Done LCO targs</th></tr>
    """

    carton_table_row = """<tr><td>{carton}</td> <td>{rs_apo}</td> <td>{done_apo}</td>
    <td>{rs_lco}</td> <td>{done_lco}</td> </tr>"""

    next_tab = """</tbody></table>

    <table><tbody>

    <tr> <th>cadence</th> <th>RS planned APO field_exps</th> <th>Done APO field_exps</th>
    <th>RS planned LCO field_exps</th> <th>Done LCO fiefield_expslds</th></tr>
    """

    field_table_row = """<tr><td>{cadence}</td> <td>{rs_apo}</td> <td>{done_apo}</td>
    <td>{rs_lco}</td> <td>{done_lco}</td> </tr>"""

    tail = """</tbody></table>
    </body></html>"""

    carton_counts = countCartons(v_base, plan, rs_base)
    tabulated_apo = fieldCounts(v_base, plan, rs_base, loc="apo")
    tabulated_lco = fieldCounts(v_base, plan, rs_base, loc="lco")

    for c in carton_counts:
        header += carton_table_row.format(**c)
    
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
    
    with open(f"{v_base}/obsSimStats-{plan}.html", "w") as webPage:
        print(header + tail, file=webPage)
