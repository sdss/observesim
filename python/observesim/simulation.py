#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import numpy as np
from astroplan import Observer
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

import observesim.weather
import roboscheduler.scheduler
import observesim.observe


def sortFields(fieldids, nexps, priorities, exp, maxTime=0):
    dtype = [('field', int), ('nexp', int), ('priority', float), ("exp", float)]
    values = [(f, n, p, a) for f, n, p, a in zip(fieldids, nexps, priorities, exp)]
    arrified = np.array(values, dtype=dtype)
    sortedFields = np.sort(arrified, order='priority')

    for f in sortedFields:
        if f["exp"] * f["nexp"] < maxTime:
            return f["field"], f["nexp"]

    return -1, -1


def apoCheck(alt, az):
    enc = [a > 45 or z > 100 or z < 80 for a, z in zip(alt, az)]
    alt = [a < 86 and a > 30 for a in alt]

    return np.logical_and(enc, alt)


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

    def __init__(self, plan, observatory, idx=1, schedule="normal"):
        if(observatory == 'apo'):
            timezone = "US/Mountain"
            fclear = 0.5
            elev = 2788
            self.telescope = {"alt": 30, "az": 90, "par_angle": 0,
                              "alt_slew": 1.5, "az_slew": 2.0, "rot_slew": 2.0}
            self.obsCheck = apoCheck
            self.moveTelescope = self.moveSloanTelescope
        if(observatory == 'lco'):
            timezone = "US/Eastern"
            fclear = 0.7
            elev = 2134
            self.telescope = {"ra": 0, "dec": -30}
            self.obsCheck = lcoCheck
            self.moveTelescope = self. moveDuPontTelescope

        self.obsHist = {"lst": list(),
                        "ra": list(),
                        "bright": list(),
                        "fieldid": list()}

        self.scheduler = roboscheduler.scheduler.Scheduler(observatory=observatory,
                                                           schedule=schedule)
        self.weather = observesim.weather.Weather(mjd_start=self.scheduler.start,
                                                  mjd_end=self.scheduler.end,
                                                  seed=idx, fclear=fclear)
        self.observatory = Observer(longitude=self.scheduler.longitude * u.deg,
                                    latitude=self.scheduler.latitude*u.deg,
                                    elevation=elev*u.m, name=observatory, timezone=timezone)
        self.scheduler.initdb(designbase=plan)
        self.field_ra = self.scheduler.fields.racen
        self.field_dec = self.scheduler.fields.deccen
        self.fieldid = self.scheduler.fields.fieldid

        cadencelist = self.scheduler.fields.cadencelist.cadences
        cadences = self.scheduler.fields.cadence

        self.nom_duration = np.float32(15. / 60. / 24.)
        self.cals = np.float32(3. / 60. / 24.)
        self.observe = observesim.observe.Observe(defaultExp=self.nom_duration,
                                                  cadencelist=cadencelist, cadences=cadences)
        self.bossReadout = np.float32(70. / 60. / 60. / 24.)

        self.curr_mjd = np.float32(1e9)

        self.coord = SkyCoord(self.field_ra * u.deg, self.field_dec * u.deg)

        self.hit_lims = 0

    def moveDuPontTelescope(self, mjd, fieldid):
        next_ra, next_dec = self.field_ra[fieldid], self.field_dec[fieldid]

        dec_slew = np.abs(next_ra-self.telescope["ra"])
        ra_slew = np.abs(next_dec-self.telescope["dec"])

        dec_time = decTime(dec_slew)
        ra_time = raTime(ra_slew)

        self.telescope["ra"] = self.field_ra[fieldid]
        self.telescope["dec"] = self.field_dec[fieldid]

        return max([dec_time, ra_time])


    def moveSloanTelescope(self, mjd, fieldid):
        altaz = self.observatory.altaz(Time(mjd, format="mjd"), self.coord[fieldid])
        alt = altaz.alt.deg
        az = altaz.az.deg
        angle = self.observatory.parallactic_angle(Time(mjd, format="mjd"), self.coord[fieldid]).deg

        alt_slew = np.abs(alt-self.telescope["alt"])
        az_slew = np.abs(az-self.telescope["az"])
        rot_slew = np.abs(angle-self.telescope["par_angle"])

        alt_time = alt_slew / self.telescope["alt_slew"]
        az_time = az_slew / self.telescope["az_slew"]
        rot_time = rot_slew / self.telescope["rot_slew"]

        self.telescope["alt"] = alt
        self.telescope["az"] = az
        self.telescope["par_angle"] = angle

        return max([alt_time, az_time, rot_time])


    def siteObs(self, fieldid, mjd):
        """Check observability issues at site, e.g. zenith at APO
           or enclosure, etc
           for any number of mjds, e.g. for a whole observing window
        """
        # try:
        #     len(mjd)
        # except TypeError:
        #     mjd = [mjd]
        # good = list()
        # for m in mjd:
        #     # print(m, fieldid)
        #     # print(self.field_ra[fieldid], self.field_dec[fieldid])
        #     altaz = self.observatory.altaz(Time(0m, format="mjd"), self.coord[fieldid])
        #     alt = altaz.alt.deg
        #     az = altaz.az.deg
        #     good.append(self.obsCheck(alt, az))

        # return np.all(good, axis=0)

        # this way seems to be at least 30% faster?
        try:
            len(mjd)
        except TypeError:
            mjd = np.array([mjd])

        try:
            len(fieldid)
        except TypeError:
            fieldid = np.array([fieldid])

        altaz = self.observatory.altaz(Time(mjd, format="mjd"), self.coord[fieldid],
                                       grid_times_targets=True)
        # altaz shape = (fields x mjds)
        alt = altaz.alt.deg.flatten()
        az = altaz.az.deg.flatten()
        res = self.obsCheck(alt, az)
        good = res.reshape((len(fieldid), len(mjd)))

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
            return -1, 1
        fieldid, nexposures = self.scheduler.nextfield(mjd=self.curr_mjd,
                                                       maxExp=maxExp)
        if(fieldid is not None):
            fieldidx = np.where(self.fieldid == fieldid)[0]
            new_alt, new_az = self.scheduler.radec2altaz(mjd=self.curr_mjd,
                                            ra=self.field_ra[fieldidx],
                                            dec=self.field_dec[fieldidx])
            new_duration = self.nom_duration * (1/np.cos(np.pi *
                                                (90-new_alt) / 180.)) ** airmass_weight
            # final_alt, final_az = self.scheduler.radec2altaz(
            #                         mjd=self.curr_mjd + new_duration * nexposures,
            #                         ra=self.field_ra[fieldid],
            #                         dec=self.field_dec[fieldid])

            site_check = self.siteObs(fieldidx, [self.curr_mjd + n*(new_duration) for n in range(nexposures)])
            maxTime = self.nextchange - self.curr_mjd
            if maxTime < new_duration * nexposures or \
                not site_check:
                # print("first field: ")
                # print(fieldid, site_check, new_alt, new_az)
                field_idxs, nexps, priorities = self.scheduler.nextfield(mjd=self.curr_mjd,
                                                                         maxExp=maxExp,
                                                                         returnAll=True)
                obs_fields = self.siteObs(field_idxs, [self.curr_mjd + n*(new_duration) for n in range(nexposures)])
                field_idxs = field_idxs[obs_fields]
                nexps = nexps[obs_fields]
                if len(field_idxs) == 0:
                    # print("all fields collide with something :( ")
                    # print(obs_fields)
                    self.hit_lims += 1./20
                    # self.curr_mjd = self.curr_mjd + self.nom_duration/20
                    return -1, 1./20
                alts, azs = self.scheduler.radec2altaz(mjd=self.curr_mjd,
                                                       ra=self.field_ra[field_idxs],
                                                       dec=self.field_dec[field_idxs])
                adj_exp = self.nom_duration * np.power((1/np.cos(np.pi *
                                                       (90 - alts) / 180.)), airmass_weight)
                fieldidx, nexposures = sortFields(field_idxs, nexps, priorities, adj_exp,
                                                  maxTime=maxTime)
                if fieldidx == -1:
                    # print("baawaaaaaahhhahahaa :( ")
                    # self.curr_mjd = self.curr_mjd + self.nom_duration/20
                    return -1, 1./20
            fieldid = int(self.fieldid[fieldidx])
            return fieldid, nexposures
        else:
            return -1, 1
        # else:
        #     lsts.append(self.scheduler.lst(self.curr_mjd)[0])
        #     observed.append(-1)
        #     self.curr_mjd = self.curr_mjd + duration
        #     exp_tonight += duration

    def observeField(self, fieldid, nexposures):

        # slewtime = np.float32(2. / 60. / 24.) # times some function of angular distance?
        slewtime = self.moveTelescope(self.curr_mjd, fieldid)
        # print("moved to {} in {} s".format(fieldid, slewtime))

        # slewtime is in seconds...
        self.curr_mjd = self.curr_mjd + self.cals + np.float32(slewtime / 60. / 60. / 24.)

        fieldidx = int(np.where(self.fieldid == fieldid)[0])

        for i in range(nexposures):
            # add each exposure
            alt, az = self.scheduler.radec2altaz(mjd=self.curr_mjd,
                                                 ra=self.field_ra[fieldidx],
                                                 dec=self.field_dec[fieldidx])
            airmass = 1/np.cos(np.pi * (90-alt) / 180.)
            # observed.append(self.scheduler.fields.racen[fieldid])
            if alt < 20:
                print(i, alt, az, self.curr_mjd, fieldid)
                if alt < 0:
                    print("booooooooo")
                    # assert False, "ugh"

            result = self.observe.result(mjd=self.curr_mjd, fieldid=fieldidx,
                                         airmass=airmass,
                                         epochidx=self.scheduler.fields.icadence[fieldidx])
            duration = result["duration"]
            if duration < 0 or np.isnan(duration):
                print("HOOOWWWOWOWOWOWW")
                print(i, alt, az, self.curr_mjd, fieldid)

            self.curr_mjd = self.curr_mjd + duration + self.bossReadout
            # lsts.append(self.scheduler.lst(self.curr_mjd)[0])

            self.obsHist["lst"].append(self.scheduler.lst(self.curr_mjd)[0])
            self.obsHist["ra"].append(self.field_ra[fieldidx])
            self.obsHist["bright"].append(self.bright())
            self.obsHist["fieldid"].append(fieldid)

            if(self.curr_mjd > self.nextchange):
                oops = (self.curr_mjd - self.nextchange) * 24 * 60
                if oops > 5:
                    print("NOOOO! BAD!", oops)
                    print(i, nexposures, alt)
            # print("observed: ", i,  fieldid, result)
            self.scheduler.update(fieldid=fieldid, result=result)

            # exp_tonight += duration

    def observeMJD(self, mjd):
        # uncomment to do a quick check
        exp_tonight = 0
        mjd_evening_twilight = self.scheduler.evening_twilight(mjd)
        mjd_morning_twilight = self.scheduler.morning_twilight(mjd)
        night_len = mjd_morning_twilight - mjd_evening_twilight
        self.curr_mjd = mjd_evening_twilight
        int_mjd = int(self.curr_mjd)
        if int_mjd % 100 == 0:
            print("!!!!", int_mjd)
        # print(int_mjd)
        while(self.curr_mjd < mjd_morning_twilight and
              self.curr_mjd < self.scheduler.end_mjd()):
            # should we do this now?
            isclear, nextchange_weather = self.weather.clear(mjd=self.curr_mjd)
            onoff, nextchange_on = self.scheduler.on(mjd=self.curr_mjd)
            nextchange_all = np.array([nextchange_weather, nextchange_on, mjd_morning_twilight])
            self.nextchange = np.min(nextchange_all)
            if not ((isclear == True) and (onoff == 'on')):
                if(self.nextchange > self.curr_mjd):
                    self.curr_mjd = self.nextchange
                continue
            fieldid, nexposures = self.nextField()
            # print(self.curr_mjd, fieldid, nexposures)
            if fieldid == -1:
                self.curr_mjd = self.curr_mjd + self.nom_duration * nexposures
                continue
            self.observeField(fieldid, nexposures)

        # weather_used[mjd] = {"length": night_len, "observed": exp_tonight}

    def lstToArray(self):
        dtype = [('lst', np.float64),
                   ('ra', np.float64),
                   ('bright', np.bool_),
                   ('fieldid', np.int32)]
        lstOut = np.zeros(len(self.obsHist["lst"]), dtype=dtype)
        lstOut["lst"] = np.array(self.obsHist["lst"])
        lstOut["ra"] = np.array(self.obsHist["ra"])
        lstOut["bright"] = np.array(self.obsHist["bright"])
        lstOut["fieldid"] = np.array(self.obsHist["fieldid"])
        return(lstOut)
