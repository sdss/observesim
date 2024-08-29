#!/usr/bin/env python
# coding: utf-8
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colorbar as cb
from matplotlib.patches import Wedge, Ellipse, Rectangle
import fitsio

import astropy.coordinates as coord
import astropy.units as u
from astropy.time import Time

from roboscheduler.moonphase import moonphase2


def countFieldMjd(obs=None, field_ids=None, mjd=None):
    done = np.extract(obs["mjd"] < mjd, obs)
    fields, counts = np.unique(done["field_pk"], return_counts=True)

    field_counts = {f: c for f, c in zip(fields, counts)}

    return np.array([field_counts.get(f, 0) for f in field_ids])


def drawMoonish(phase, ax):
    # 0 new, 4 full, 8 would be new again
    if phase == 0:
        el_col = "None"
        wed_2_col = "None"
        wed_1_col = "None"
    elif phase == 1:
        el_col = "w"
        wed_2_col = "None"
        wed_1_col = "k"
    elif phase == 2:
        el_col = "None"
        wed_2_col = "None"
        wed_1_col = "k"
    elif phase == 3:
        el_col = "k"
        wed_2_col = "None"
        wed_1_col = "k"
    elif phase == 4:
        el_col = "None"
        wed_2_col = "k"
        wed_1_col = "k"
    elif phase == 5:
        el_col = "k"
        wed_2_col = "k"
        wed_1_col = "None"
    elif phase == 6:
        el_col = "None"
        wed_2_col = "k"
        wed_1_col = "None"
    elif phase == 7:
        el_col = "w"
        wed_2_col = "k"
        wed_1_col = "None"

    wed1 = Wedge((0.0, 0.0), r=1, theta1=90, theta2=270, fc=wed_1_col)
    wed2 = Wedge((0.0, 0.0), r=1, theta1=270, theta2=90, fc=wed_2_col)
    el = Ellipse((0.0, 0.0), width=1, height=2, fc=el_col)
    ax.add_artist(wed1)
    ax.add_artist(wed2)
    ax.add_artist(el)

    ax.axis([-1, 1, -1, 1])


def CountFramesAllSky(base, plan, version=None):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    movieFrameDir = v_base + "moviePngsAllSky"

    try:
        os.makedirs(movieFrameDir)
    except:
        pass

    doneFrames = os.listdir(movieFrameDir)

    # start with fresh slate
    if len(doneFrames) > 0:
        for f in doneFrames:
            os.remove(movieFrameDir + "/" + f)

    apo_observations = fitsio.read(v_base + f"{plan}-apo-observations-0.fits",
                                   columns=["field_pk", "mjd", "skybrightness"])
    apo_observe_fields = fitsio.read(v_base + f"{plan}-apo-fields-0.fits")

    lco_observations = fitsio.read(v_base + f"{plan}-lco-observations-0.fits",
                                   columns=["field_pk", "mjd", "skybrightness"])
    lco_observe_fields = fitsio.read(v_base + f"{plan}-lco-fields-0.fits")

    # store the ra & decs of the fields as sky coordinate objects
    apo_ra = coord.Angle(-(apo_observe_fields['racen']+90)*u.degree)
    apo_both_ra = apo_ra.wrap_at(180*u.degree)
    apo_both_dec = coord.Angle(apo_observe_fields['deccen']*u.degree)

    lco_ra = coord.Angle(-(lco_observe_fields['racen']+90)*u.degree)
    lco_both_ra = lco_ra.wrap_at(180*u.degree)
    lco_both_dec = coord.Angle(lco_observe_fields['deccen']*u.degree)

    apo_field_ids = apo_observe_fields['pk']
    lco_field_ids = lco_observe_fields['pk']

    # find the min and max mjds for this set of observations
    min_mjd = np.floor(min(np.min(apo_observations['mjd']), np.min(lco_observations['mjd'])))
    max_mjd = np.ceil(max(np.max(apo_observations['mjd']), np.max(lco_observations['mjd'])))

    survey_start = Time(min_mjd, format="mjd")
    start_month = survey_start.datetime.strftime("%B %Y")
    survey_end = Time(max_mjd, format="mjd")
    end_month = survey_end.datetime.strftime("%B %Y")

    survey_length = (max_mjd - min_mjd)

    print(f"counting fields obs from {min_mjd} to {max_mjd}; {survey_length/365.} years")

    # count the number of nights & fields in this simulation
    n_nights = np.ceil(max_mjd) - np.floor(min_mjd)
    n_apo_fields = len(apo_observe_fields)
    n_lco_fields = len(lco_observe_fields)

    # make an array of all the nights
    mjds = np.arange(min_mjd, min_mjd+n_nights, 1)

    scaleFunc = np.arctan

    # set a ceiling to allow the visit maps to saturate, so that RM fields don't break things
    ceil = 10.
    dark_ceil = scaleFunc(ceil)

    bright_apo = np.extract(apo_observations["skybrightness"] > 0.35, apo_observations)
    dark_apo = np.extract(apo_observations["skybrightness"] <= 0.35, apo_observations)
    bright_lco = np.extract(lco_observations["skybrightness"] > 0.35, lco_observations)
    dark_lco = np.extract(lco_observations["skybrightness"] <= 0.35, lco_observations)

    # calculate moon phases for everything
    days = np.arange(min_mjd - 1, max_mjd + 1, 0.1)
    phase = moonphase2(days)

    previous = moonphase2(days[0] - 1)

    mjd_to_phase = dict()
    for p, d in zip(phase, days):
        if p - previous < 0:
            waxing = False
        else:
            waxing = True
        if waxing:
            mapped = p*4
        else:
            mapped = 8 - (p*4)

        previous = p
        mjd_to_phase[int(d)] = np.floor(mapped)

    # min_mjd = max_mjd - 45
    for i, mjd in enumerate(np.arange(min_mjd, max_mjd, 2)):
        w, h = 16, 9
        if i % 50 == 0:
            print(f"working on frame: {i}, mjd: {mjd}")
        f = plt.figure(figsize=(w, h))
        nrows = 5
        ax1 = plt.subplot2grid((nrows, 1), (0, 0), projection='mollweide', rowspan=nrows-1)

        ax3 = plt.subplot2grid((nrows, 1), (nrows-1, 0))

        # count visits
        apo_visits = countFieldMjd(obs=apo_observations, field_ids=apo_field_ids, mjd=mjd)
        lco_visits = countFieldMjd(obs=lco_observations, field_ids=lco_field_ids, mjd=mjd)

        cmap = "hot_r"
        mk = "D"

        im = ax1.scatter(apo_both_ra.radian, apo_both_dec.radian, s=2.7, c=scaleFunc(apo_visits),
                         cmap=cmap, vmin=0, vmax=dark_ceil, marker=mk)
        ax1.scatter(lco_both_ra.radian, lco_both_dec.radian, s=0.8, c=scaleFunc(lco_visits),
                    cmap=cmap, vmin=0, vmax=dark_ceil, marker=mk)

        from_edge = 0.1
        width = 0.1
        # [x. y, width, heigt]
        axin1 = ax1.inset_axes([1-from_edge, 1-from_edge, width, width*(w/h)])
        axin1.get_xaxis().set_visible(False)
        axin1.get_yaxis().set_visible(False)
        axin1.axis("off")
        drawMoonish(mjd_to_phase[int(mjd)], axin1)

        axin2 = ax1.inset_axes([0.35, 1.05, 0.3, 0.05])
        axin2.get_yaxis().set_visible(False)
        curr = (mjd-min_mjd)/survey_length
        bar = Rectangle((0, 0), curr, 1, fc="DarkOrange")
        axin2.add_artist(bar)
        axin2.set_xticks([0, 1])
        axin2.set_xticklabels([start_month, end_month])

        ax3.get_xaxis().set_visible(False)
        ax3.get_yaxis().set_visible(False)
        ax3.axis("off")

        f.suptitle(f"MJD = {mjd:8.2f}", fontsize=15)

        # Fine-tune figure; make subplots farther from each other.
        f.subplots_adjust(hspace=0.3)

        cbar = f.colorbar(im, ax=[ax3], location="bottom", ticks=[0, dark_ceil])

        cbar.ax.set_xticklabels(['Less Observations', 'More Observations'])

    #     f.savefig(plotname, dpi = 300)
        plt.savefig(f"{movieFrameDir}/frame{i:03d}.png", dpi=150)

        plt.close()

    return movieFrameDir


def CountFrames(base, plan, version=None, variant=""):
    if version is not None:
        v_base = os.path.join(base, version)
        v_base += "/"
    else:
        v_base = os.path.join(base, plan)
        v_base += "/"

    movieFrameDir = v_base + "moviePngs" + variant

    try:
        os.makedirs(movieFrameDir)
    except:
        pass

    doneFrames = os.listdir(movieFrameDir)

    # start with fresh slate
    if len(doneFrames) > 0:
        for f in doneFrames:
            os.remove(movieFrameDir + "/" + f)

    apo_observations = fitsio.read(v_base + f"{plan}-apo-observations-0.fits",
                                   columns=["field_pk", "mjd", "skybrightness"])
    apo_observe_fields = fitsio.read(v_base + f"{plan}-apo-fields-0.fits")

    lco_observations = fitsio.read(v_base + f"{plan}-lco-observations-0.fits",
                                   columns=["field_pk", "mjd", "skybrightness"])
    lco_observe_fields = fitsio.read(v_base + f"{plan}-lco-fields-0.fits")

    # store the ra & decs of the fields as sky coordinate objects
    apo_ra = coord.Angle(apo_observe_fields['racen']*u.degree)
    apo_both_ra = apo_ra.wrap_at(180*u.degree)
    apo_both_dec = coord.Angle(apo_observe_fields['deccen']*u.degree)

    lco_ra = coord.Angle(lco_observe_fields['racen']*u.degree)
    lco_both_ra = lco_ra.wrap_at(180*u.degree)
    lco_both_dec = coord.Angle(lco_observe_fields['deccen']*u.degree)

    apo_field_ids = apo_observe_fields['field_pk']
    lco_field_ids = lco_observe_fields['field_pk']

    # find the min and max mjds for this set of observations
    min_mjd = min(np.min(apo_observations['mjd']), np.min(lco_observations['mjd']))
    max_mjd = max(np.max(apo_observations['mjd']), np.max(lco_observations['mjd']))

    yrs = (np.ceil(max_mjd) - np.floor(min_mjd))/365.

    print("counting fields obs from {min_mjd} to {max_mjd}; {yrs} years")

    # count the number of nights & fields in this simulation
    n_nights = np.ceil(max_mjd) - np.floor(min_mjd)  # math.ceil(max_mjd - min_mjd)
    n_apo_fields = len(apo_observe_fields)
    n_lco_fields = len(lco_observe_fields)

    # make an array of all the nights
    mjds = np.arange(min_mjd, min_mjd+n_nights, 1)

    scaleFunc = np.arctan

    # set a ceiling to allow the visit maps to saturate, so that RM fields don't break things
    bright_ceil = 20.
    dark_ceil = scaleFunc(20.)

    bright_apo = np.extract(apo_observations["skybrightness"] > 0.35, apo_observations)
    dark_apo = np.extract(apo_observations["skybrightness"] <= 0.35, apo_observations)
    bright_lco = np.extract(lco_observations["skybrightness"] > 0.35, lco_observations)
    dark_lco = np.extract(lco_observations["skybrightness"] <= 0.35, lco_observations)

    for i, mjd in enumerate(np.arange(min_mjd, max_mjd, 7)):
        if i % 50 == 0:
            print(f"working on frame: {i}, mjd: {mjd}")
        if variant == "square":
            f = plt.figure(figsize=(20, 20))
        else:
            f = plt.figure(figsize=(20, 10))
        nrows = 3
        ax1 = plt.subplot2grid((nrows, 2), (0, 0), projection='mollweide', rowspan=nrows-1)
        ax2 = plt.subplot2grid((nrows, 2), (0, 1), projection='mollweide', rowspan=nrows-1)
        ax3 = plt.subplot2grid((nrows, 2), (nrows-1, 0))
        ax4 = plt.subplot2grid((nrows, 2), (nrows-1, 1))

        # count visits
        apo_bright_visits = countFieldMjd(obs=bright_apo, field_ids=apo_field_ids, mjd=mjd)
        apo_dark_visits = countFieldMjd(obs=dark_apo, field_ids=apo_field_ids, mjd=mjd)
        lco_bright_visits = countFieldMjd(obs=bright_lco, field_ids=lco_field_ids, mjd=mjd)
        lco_dark_visits = countFieldMjd(obs=dark_lco, field_ids=lco_field_ids, mjd=mjd)

        cmap = "hot_r"
        mk = "D"

        ax1.scatter(apo_both_ra.radian, apo_both_dec.radian, s=2.7, c=scaleFunc(apo_dark_visits),
                    cmap=cmap, vmin=0, vmax=dark_ceil, marker=mk)
        ax1.scatter(lco_both_ra.radian, lco_both_dec.radian, s=0.8, c=scaleFunc(lco_dark_visits),
                    cmap=cmap, vmin=0, vmax=dark_ceil, marker=mk)
        ax1.set_title('Dark Time', fontsize=10)

        im = ax2.scatter(apo_both_ra.radian, apo_both_dec.radian, s=2.7, c=scaleFunc(apo_bright_visits),
                         cmap=cmap, vmin=0, vmax=dark_ceil, marker=mk)
        ax2.scatter(lco_both_ra.radian, lco_both_dec.radian, s=0.8, c=scaleFunc(lco_bright_visits),
                    cmap=cmap, vmin=0, vmax=dark_ceil, marker=mk)
        ax2.set_title('Bright Time', fontsize=10)

        for ax in [ax3, ax4]:
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.axis("off")

        f.suptitle(f"MJD = {mjd:8.2f}", fontsize=15)

        # Fine-tune figure; make subplots farther from each other.
        f.subplots_adjust(hspace=0.3)

        cbar = f.colorbar(im, ax=[ax3, ax4], location="bottom", ticks=[0, dark_ceil])

        cbar.ax.set_xticklabels(['Less Observations', 'More Observations'])

    #     f.savefig(plotname, dpi = 300)
        plt.savefig(f"{movieFrameDir}/frame{i:03d}.png", dpi=150)

        plt.close()

    return movieFrameDir


if __name__ == "__main__":
    pass
    # plan = "gamma-0"
    # base = "/Users/jdonor/software/sdss/simCache/"

    # # movieFrameDir = CountFrames(base, plan)
    # movieFrameDir = CountFramesAllSky(base, plan)

    # # now makes mp4 file from frames in FramesForMP4 dir
    # os.chdir(movieFrameDir)
    # subprocess.call('ffmpeg -r 10 -i frame%03d.png -vcodec mpeg4 -y sdss5_sim.mp4', shell=True)

    # print(f'Get you mp4 video in {movieFrameDir}!')
