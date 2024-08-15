#!/usr/bin/env python
# encoding: utf-8

import os

import numpy as np
from astroplan import Observer
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import yaml

import observesim.weather
import roboscheduler.scheduler
import observesim.observe
from observesim.progress import fieldsFromDB, doneForObs


def sortFields(fieldids, nexps, exp, maxTime=0):

    for f, n in zip(fieldids, nexps):
        if exp * n < maxTime:
            return f, n

    return -1, -1


def apoCheck(alt, az):
    enc = [a > 45 or z > 100 or z < 80 for a, z in zip(alt, az)]
    altc = [a < 86 and a > 30 for a in alt]

    return np.logical_and(enc, altc)


def lcoCheck(alt, az):
    return alt > 30


def accSlewTime(degrees):
    # compute time for Du Pont with acc/decceleration considered
    # starting with 1/2*a*t_1^2 + a*t_1*t_2 + 1/2*d*t_2^2 = dist
    # and a*t_1 = d*t_2, its trivial to solve

    acc = 0.022

    return 3 * np.sqrt(2 * degrees / (7 * acc))


def decTime(degrees):
    # time for Du Pont to move on dec axis
    if degrees < 61.6:
        return accSlewTime(degrees)
    else:
        return (degrees - 61.6)/0.63 + 84.4


def raTime(degrees):
    # time for Du Pont to move on RA axis
    if degrees < 37.2:
        return accSlewTime(degrees)
    else:
        return (degrees - 37.2)/0.49 + 65.6


class Simulation(object):
    """A class to encapsulate an SDSS-5 simulation
    """

    def __init__(self, plan, observatory, idx=1, schedule="normal", redo_exp=True,
                 oldWeather=False, with_hist=False, rsFinal=True,
                 hist_plan="eta-9"):

        out_path = os.getenv('RS_OUTDIR')
        cfg_file = os.path.join(out_path, "sim_cfg.yml")
        if os.path.isfile(cfg_file):
            print(f"found sim config file, applying: \n {cfg_file}")
            cfg = yaml.load(open(cfg_file), Loader=yaml.FullLoader)
        else:
            cfg_file = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + "/etc/nominal_cfg.yml"
            cfg = yaml.load(open(cfg_file), Loader=yaml.FullLoader)

        if(observatory == 'apo'):
            timezone = "US/Mountain"
            fclear = 0.5
            elev = 2788
            self.telescope = {"alt": 30, "az": 90, "par_angle": 0,
                              "alt_slew": 1.5, "az_slew": 2.0, "rot_slew": 2.0}
            self.obsCheck = apoCheck
            self.moveTelescope = self.moveSloanTelescope
            self.nom_duration = np.float32(cfg["nom_duration_apo"] / 60. / 24.)
            # self.cals = np.float32(3. / 60. / 24.)
            self.field_overhead = np.float32(cfg["field_overhead_apo"] / 60. / 24.)
            self.design_overhead = np.float32(cfg["design_overhead_apo"] / 60. / 24.)
        if(observatory == 'lco'):
            timezone = "US/Eastern"
            fclear = 0.7
            elev = 2134
            self.telescope = {"ra": 0, "dec": -30}
            self.obsCheck = lcoCheck
            self.moveTelescope = self.moveDuPontTelescope
            self.nom_duration = np.float32(cfg["nom_duration_lco"] / 60. / 24.)
            # self.cals = np.float32(3. / 60. / 24.)
            # extra 90s or so for cals at LCO?
            self.field_overhead = np.float32(cfg["field_overhead_lco"] / 60. / 24.)
            self.design_overhead = np.float32(cfg["design_overhead_lco"] / 60. / 24.)

        self.redo_exp = redo_exp
        self.snr_b = cfg["snr_b"]
        self.snr_r = cfg["snr_r"]
        self.snr_ap = cfg["snr_ap"]

        self.obsHist = {"lst": list(),
                        "ra": list(),
                        "bright": list(),
                        "field_pk": list(),
                        "weather": list(),
                        "mjd": list()}

        out_path = os.getenv('RS_OUTDIR')
        priority_file = os.path.join(out_path, "priority.yml")
        if os.path.isfile(priority_file):
            print(f"found priority file, applying: \n {priority_file}")
            priorities = yaml.load(open(priority_file), Loader=yaml.FullLoader)
        else:
            priorities = dict()

        self.scheduler = roboscheduler.scheduler.Scheduler(observatory=observatory,
                                                           schedule=schedule,
                                                           priorities=priorities)

        if observatory == "lco" and not oldWeather:
            base = os.getenv("OBSERVESIM_OUTPUT_BASE")
            modelsDir = os.path.join(base, "weather_models")
            fname = os.path.join(modelsDir, f"saved_model_{observatory}_{idx}.csv")
            self.weather = observesim.weather.Weather3(mjd_start=self.scheduler.start,
                                                       mjd_end=self.scheduler.end,
                                                       model_fname=fname)
        elif not oldWeather:
            self.weather = observesim.weather.Weather2(mjd_start=self.scheduler.start,
                                                       mjd_end=self.scheduler.end,
                                                       seed=idx, loc=observatory)

        else:
            self.weather = observesim.weather.Weather(mjd_start=self.scheduler.start,
                                                      mjd_end=self.scheduler.end,
                                                      seed=idx, fclear=fclear)

        self.observatory = Observer(longitude=self.scheduler.longitude * u.deg,
                                    latitude=self.scheduler.latitude*u.deg,
                                    elevation=elev*u.m, name=observatory, timezone=timezone)
        if with_hist:
            fields_sum, all_designs = fieldsFromDB(obs=observatory.upper(),
                                                   plan=hist_plan)
            self.scheduler.initdb(designbase=hist_plan, 
                                  fieldsArray=fields_sum,
                                  realDesigns=all_designs,
                                  fromFits=False)
        else:
            self.scheduler.initdb(designbase=plan, rsFinal=rsFinal)
        self.field_ra = self.scheduler.fields.racen
        self.field_dec = self.scheduler.fields.deccen
        self.field_pk = self.scheduler.fields.pk

        self.observe = observesim.observe.Observe(defaultExp=self.nom_duration,
                                                  loc=observatory)
        # self.bossReadout = np.float32(70. / 60. / 60. / 24.)

        self.curr_mjd = np.float32(1e9)
        self.nextchange_weather = 5e9

        self.coord = SkyCoord(self.field_ra * u.deg, self.field_dec * u.deg)

        self.slews = list()
        self.slew_mjds = list()
        self.slew_alt = list()
        self.slew_az = list()
        self.slew_rot = list()
        self.slew_ra = list()
        self.slew_dec = list()

        self.hit_lims = 0

        self.redo_apg = 0
        self.redo_r = 0
        self.redo_b = 0

        if with_hist:
            field_mjds = doneForObs(obs=observatory.upper(),
                                    plan=hist_plan)

            result = {"mjd":-99999,
                      "duration": 900,
                      "apgSN2": -1.0,
                      "rSN2": -1.0,
                      "bSN2": -1.0}
            for f in field_mjds:
                result["mjd"] = f["mjd"]
                self.scheduler.update(field_pk=f["field_pk"], result=result,
                                      finish=True)

    def moveDuPontTelescope(self, mjd, fieldidx):
        next_ra, next_dec = self.field_ra[fieldidx], self.field_dec[fieldidx]

        ra_slew = np.abs(next_ra-self.telescope["ra"])
        dec_slew = np.abs(next_dec-self.telescope["dec"])

        if ra_slew > 180:
            ra_slew = 360 - ra_slew
            assert ra_slew > 0, "forgot circular math? ra"

        dec_time = decTime(dec_slew)
        ra_time = raTime(ra_slew)

        self.telescope["ra"] = self.field_ra[fieldidx]
        self.telescope["dec"] = self.field_dec[fieldidx]

        return max([dec_time, ra_time]), ra_slew, dec_slew

    def moveSloanTelescope(self, mjd, fieldidx):
        altaz = self.observatory.altaz(Time(mjd, format="mjd"), self.coord[fieldidx])
        alt = altaz.alt.deg
        az = altaz.az.deg
        angle = self.observatory.parallactic_angle(Time(mjd, format="mjd"), self.coord[fieldidx]).deg

        alt_slew = np.abs(alt-self.telescope["alt"])
        az_slew = np.abs(az-self.telescope["az"])
        if az_slew > 180:
            az_slew = 360 - az_slew
            assert az_slew > 0, "forgot circular math?  az"
        rot_slew = np.abs(angle-self.telescope["par_angle"])
        if rot_slew > 180:
            rot_slew = 360 - rot_slew
            assert rot_slew > 0, "forgot circular math? rot"

        alt_time = alt_slew / self.telescope["alt_slew"]
        az_time = az_slew / self.telescope["az_slew"]
        rot_time = rot_slew / self.telescope["rot_slew"]

        self.telescope["alt"] = alt
        self.telescope["az"] = az
        self.telescope["par_angle"] = angle

        return max([alt_time, az_time, rot_time]), alt_slew, az_slew, rot_slew

    def siteObs(self, fieldidx, mjd):
        """Check observability issues at site, e.g. zenith at APO
           or enclosure, etc
           for any number of mjds, e.g. for a whole observing window
        """

        try:
            len(mjd)
        except TypeError:
            mjd = np.array([mjd])

        try:
            len(fieldidx)
        except TypeError:
            fieldidx = np.array([fieldidx])

        altaz = self.observatory.altaz(Time(mjd, format="mjd"), self.coord[fieldidx],
                                       grid_times_targets=True)
        # altaz shape = (fields x mjds)
        alt = altaz.alt.deg.flatten()
        az = altaz.az.deg.flatten()
        res = self.obsCheck(alt, az)
        good = res.reshape((len(fieldidx), len(mjd)))

        # axis 1 is along fields, I guess...
        return np.all(good, axis=1)

    def bright(self, mjd=None):
        if mjd is None:
            mjd = self.curr_mjd
        skybrightness = self.scheduler.skybrightness(mjd)
        return skybrightness > 0.35

    def nextField(self):
        # dark time or brighttime? to guess at how long we need for obs
        if not self.bright():
            airmass_weight = 1.05
        else:
            airmass_weight = 0.05
        # integer division floors; no partial exposures
        maxExp = int((self.nextchange - self.curr_mjd)//(self.nom_duration * 1.3 ** airmass_weight))
        if maxExp == 0:
            # self.curr_mjd = self.curr_mjd + self.nom_duration
            return -1, 1, True
        field_pk, nexposures = self.scheduler.nextfield(mjd=self.curr_mjd,
                                                        maxExp=maxExp)
        # assert fieldid is not None, f"can't schedule {self.curr_mjd}, {self.bright()}"
        if(field_pk is not None):
            fieldidx = np.where(self.field_pk == field_pk)[0]
            site_check = self.siteObs(fieldidx, [self.curr_mjd + n*(self.nom_duration) for n in range(nexposures)])
            # maxTime = self.nextchange - self.curr_mjd
            # maxTime = maxExp * self.nom_duration
            maxTime = maxExp * (self.nom_duration)

            if not site_check:
                field_idxs, nexps = self.scheduler.nextfield(mjd=self.curr_mjd,
                                                             maxExp=maxExp,
                                                             returnAll=True)

                obs_fields = self.siteObs(field_idxs, [self.curr_mjd + n*(self.nom_duration) for n in range(nexposures)])
                field_idxs = field_idxs[obs_fields]
                nexps = nexps[obs_fields]
                if len(field_idxs) == 0:
                    # print("all fields collide with something :( ")
                    # print(obs_fields)
                    self.hit_lims += 1./20
                    return -1, 1./20, False

                fieldidx, nexposures = sortFields(field_idxs, nexps, self.nom_duration, maxTime=maxTime)
                if fieldidx == -1:
                    # print("baawaaaaaahhhahahaa :( ")
                    # self.curr_mjd = self.curr_mjd + self.nom_duration/20
                    return -1, 1./20, False
            field_pk = int(self.field_pk[fieldidx])

            return field_pk, nexposures, False
        else:
            # if not self.bright():
            #     assert False, f"{self.curr_mjd} ugh"
            return -1, 1, False

    def bookKeeping(self, fieldidx, i=-1, cloudy=False):
        """figure out SN and keep track, etc
        """
        alt, az = self.scheduler.radec2altaz(mjd=self.curr_mjd,
                                             ra=self.field_ra[fieldidx],
                                             dec=self.field_dec[fieldidx])
        airmass = 1/np.cos(np.pi * (90-alt) / 180.)
        if alt < 20:
            print(i, alt, az, self.curr_mjd, fieldidx, "TOO LOW!!")
            if alt < 0:
                print("booooooooo")
                # assert False, "ugh"

        result = self.observe.result(mjd=self.curr_mjd, field_pk=self.field_pk[fieldidx],
                                     airmass=airmass, cloudy=cloudy,
                                     epochidx=self.scheduler.fields.icadence[fieldidx])
        duration = result["duration"]
        if duration < 0 or np.isnan(duration):
            print("HOOOWWWOWOWOWOWW")
            print(i, alt, az, self.curr_mjd, field_pk)

        # self.curr_mjd = self.curr_mjd + duration + self.bossReadout
        self.curr_mjd = self.curr_mjd + duration + self.design_overhead
        # print(f"{float(duration):0.4f}, {float(self.design_overhead):0.6f}, {float(self.curr_mjd):0.7f}")

        # move telescope for tracking
        self.moveTelescope(self.curr_mjd, fieldidx)

        self.obsHist["lst"].append(self.scheduler.lst(self.curr_mjd)[0])
        self.obsHist["ra"].append(self.field_ra[fieldidx])
        self.obsHist["bright"].append(self.scheduler.skybrightness(self.curr_mjd))
        self.obsHist["field_pk"].append(self.field_pk[fieldidx])
        self.obsHist["weather"].append(False)
        self.obsHist["mjd"].append(float(self.curr_mjd))

        return result

    def observeField(self, field_pk, nexposures, cloudy=False):
        fieldidx = int(np.where(self.field_pk == field_pk)[0])

        slewtime, *axes = self.moveTelescope(self.curr_mjd, fieldidx)

        self.slews.append(int(slewtime))
        self.slew_mjds.append(int(self.curr_mjd))

        if len(axes) == 3:
            self.slew_alt.append(float(axes[0]))
            self.slew_az.append(float(axes[1]))
            self.slew_rot.append(float(axes[2]))

            self.slew_ra.append(np.nan)
            self.slew_dec.append(np.nan)
        else:
            self.slew_ra.append(float(axes[0]))
            self.slew_dec.append(float(axes[1]))

            self.slew_alt.append(np.nan)
            self.slew_az.append(np.nan)
            self.slew_rot.append(np.nan)

        # slewtime is in seconds...
        # self.curr_mjd = self.curr_mjd + self.cals + np.float32(slewtime / 60. / 60. / 24.)
        self.curr_mjd = self.curr_mjd + self.field_overhead
        # print(f"{fieldidx}, {float(self.field_overhead):0.6f}, {float(self.curr_mjd):0.7f}")

        field_exp_count = nexposures
        for i in range(nexposures):
            # each "exposure" is a design

            if(self.curr_mjd > self.nextchange):
                oops = (self.curr_mjd - self.nextchange) * 24 * 60
                if oops > 5:
                    print("NOOOO! BAD!", oops)
                    # print(i, nexposures, self.telescope)
                continue

            res = self.bookKeeping(fieldidx, i=i, cloudy=cloudy)

            if self.bright():
                if res["apgSN2"] < self.snr_ap and self.redo_exp:
                    field_exp_count += 1
                    self.redo_apg += 1
                    self.scheduler.update(field_pk=field_pk, result=res,
                                          finish=False)
                    res = self.bookKeeping(fieldidx, i=i, cloudy=cloudy)
                    self.scheduler.update(field_pk=field_pk, result=res,
                                          finish=True)
                else:
                    self.scheduler.update(field_pk=field_pk, result=res,
                                          finish=True)
            else:
                if (res["rSN2"] < self.snr_r or res["bSN2"] < self.snr_b) and self.redo_exp:
                    field_exp_count += 1
                    if res["rSN2"] < self.snr_r:
                        self.redo_r += 1
                    else:
                        self.redo_b += 1
                    self.scheduler.update(field_pk=field_pk, result=res,
                                          finish=False)
                    res = self.bookKeeping(fieldidx, i=i, cloudy=cloudy)
                    self.scheduler.update(field_pk=field_pk, result=res,
                                          finish=True)
                else:
                    self.scheduler.update(field_pk=field_pk, result=res,
                                          finish=True)
        # if self.bright():
        #     ap_tot = np.sum(self.scheduler.observations.apgSN2[-1*field_exp_count:])
        #     # print(f"{nexposures} {field_exp_count} {len(self.scheduler.observations.apgSN2[-1*field_exp_count:])}")
        #     # print(f"AP SN {ap_tot:7.1f} VS {300 * nexposures}")
        #     if ap_tot < 2025 * nexposures and self.redo_exp:
        #         self.redo_apg += 1
        #         self.scheduler.update(field_pk=field_pk, result=res,
        #                               finish=False)
        #         res = self.bookKeeping(fieldidx, i=i, cloudy=cloudy)
        #         self.scheduler.update(field_pk=field_pk, result=res,
        #                               finish=False)
        # else:
        #     r_tot = np.sum(self.scheduler.observations.rSN2[-1*field_exp_count:])
        #     b_tot = np.sum(self.scheduler.observations.bSN2[-1*field_exp_count:])
        #     # print(f"{nexposures} {field_exp_count} {len(self.scheduler.observations.bSN2[-1*field_exp_count:])}")
        #     # print(f"B  SN {b_tot:7.1f} VS {2.5 * nexposures} \nR  SN {r_tot:7.1f} VS {5 * nexposures}")
        #     if (b_tot < 2 * nexposures or r_tot < 4 * nexposures) and self.redo_exp:
        #         if r_tot < 4 * nexposures:
        #             self.redo_r += 1
        #         else:
        #             self.redo_b += 1
        #         self.scheduler.update(field_pk=field_pk, result=res,
        #                               finish=False)
        #         res = self.bookKeeping(fieldidx, i=i, cloudy=cloudy)
        #         self.scheduler.update(field_pk=field_pk, result=res,
        #                               finish=False)

    def observeMJD(self, mjd):
        mjd_evening_twilight = self.scheduler.evening_twilight(mjd)
        mjd_morning_twilight = self.scheduler.morning_twilight(mjd)
        self.curr_mjd = mjd_evening_twilight
        # int_mjd = int(self.curr_mjd)

        surveyGoal = np.sum(self.scheduler.fields.slots)
        surveyDone = np.sum([len(self.scheduler.fields.hist[i]) for i in self.scheduler.fields.pk])

        self.scheduler.surveyComplete = surveyDone / surveyGoal

        if mjd % 100 == 0:
            print("!!!!", mjd, self.scheduler.surveyComplete, " done")

        # guesses = np.arange(0, 1, 0.05)

        self.nextchange = mjd_morning_twilight

        this_moon = 0

        isclear, self.nextchange_weather, cloudy = self.weather.clear(now=self.curr_mjd,
                                                              until=mjd_morning_twilight)
        # carries on until a change or it hits morning

        while(self.curr_mjd < mjd_morning_twilight and
              self.curr_mjd < self.scheduler.end_mjd()):

            if self.nextchange_weather <= self.curr_mjd:
                isclear, self.nextchange_weather, cloudy = self.weather.clear(now=self.curr_mjd,
                                                                      until=mjd_morning_twilight)

            onoff, nextchange_on = self.scheduler.on(mjd=self.curr_mjd)

            nextchange = np.min(np.array([self.nextchange_weather, nextchange_on, mjd_morning_twilight]))

            if not isclear:
                # count = 0
                # dur = float(nextchange - self.curr_mjd)
                while self.curr_mjd < nextchange:
                    if nextchange - self.curr_mjd < self.nom_duration:
                        self.curr_mjd = nextchange
                        continue
                    self.obsHist["lst"].append(self.scheduler.lst(self.curr_mjd)[0])
                    self.obsHist["ra"].append(-1)
                    self.obsHist["bright"].append(float(self.scheduler.skybrightness(self.curr_mjd)))
                    self.obsHist["field_pk"].append(-1)
                    self.obsHist["weather"].append(True)
                    self.obsHist["mjd"].append(float(self.curr_mjd))
                    self.curr_mjd += self.nom_duration
                    # count += 1
                # print("WEATHER ", self.curr_mjd, f"night {night_len*24:.1f}, weather {dur*24:.1f}", count)
            elif (onoff != 'on'):
                self.curr_mjd = nextchange

            if self.nextchange - self.curr_mjd < self.nom_duration:
                self.curr_mjd = self.nextchange
                continue

            if self.scheduler.skybrightness(self.curr_mjd) > this_moon:
                this_moon = self.scheduler.skybrightness(self.curr_mjd)

            field_pk, nexposures, noTime = self.nextField()
            if field_pk == -1:
                if noTime:
                    self.curr_mjd = self.curr_mjd + self.nom_duration
                    continue
                if this_moon < 0.98:
                    print("skipped ", self.curr_mjd, self.scheduler.skybrightness(self.curr_mjd), this_moon)
                    # raise Exception()
                # print("skipped ", self.curr_mjd)
                self.obsHist["lst"].append(self.scheduler.lst(self.curr_mjd)[0])
                self.obsHist["ra"].append(np.nan)
                self.obsHist["bright"].append(float(this_moon))
                self.obsHist["field_pk"].append(-1)
                self.obsHist["weather"].append(False)
                self.obsHist["mjd"].append(float(self.curr_mjd))
                self.curr_mjd = self.curr_mjd + self.nom_duration
                print("nothing scheduled", self.curr_mjd, self.obsHist["mjd"][-1]-self.curr_mjd)
                continue
            self.observeField(field_pk, nexposures, cloudy=cloudy)

        # if mjd % 10 == 0:
        #     self.scheduler.priorityLogger.write(name=str(mjd) + "-" + self.observatory.name)

    def lstToArray(self):
        assert len(self.obsHist["weather"]) == len(self.obsHist["lst"]), "lst tracking bad!"
        dtype = [('lst', np.float64),
                 ('ra', np.float64),
                 ('bright', np.float32),
                 ('field_pk', np.int32),
                 ('weather', np.bool_),
                 ('mjd', np.float64)]
        lstOut = np.zeros(len(self.obsHist["lst"]), dtype=dtype)
        lstOut["lst"] = np.array(self.obsHist["lst"])
        lstOut["ra"] = np.array(self.obsHist["ra"])
        lstOut["bright"] = np.array(self.obsHist["bright"])
        lstOut["field_pk"] = np.array(self.obsHist["field_pk"])
        lstOut["weather"] = np.array(self.obsHist["weather"])
        lstOut["mjd"] = np.array(self.obsHist["mjd"])
        return(lstOut)

    def slewsToArray(self):
        dtype = [('time', np.int32),
                 ('mjd', np.int32),
                 ('alt', np.float64),
                 ('az', np.float64),
                 ('rot', np.float64),
                 ('ra', np.float64),
                 ('dec', np.float64)]
        arrayOut = np.zeros(len(self.slews), dtype=dtype)
        arrayOut["time"] = np.array(self.slews)
        arrayOut["mjd"] = np.array(self.slew_mjds)
        arrayOut["alt"] = np.array(self.slew_alt)
        arrayOut["az"] = np.array(self.slew_az)
        arrayOut["rot"] = np.array(self.slew_rot)
        arrayOut["ra"] = np.array(self.slew_ra)
        arrayOut["dec"] = np.array(self.slew_dec)

        return(arrayOut)
