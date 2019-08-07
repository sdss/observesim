#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import numpy as np

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


class Simulation(object):
    """A class to encapsulate an SDSS-5 simulation
    """

    def __init__(self, base, plan, observatory, idx=1):
        if(observatory == 'apo'):
            fclear = 0.5
        if(observatory == 'lco'):
            fclear = 0.7
        
        self.scheduler = roboscheduler.scheduler.Scheduler(observatory=observatory)
        self.weather = observesim.weather.Weather(mjd_start=self.scheduler.start,
                                             mjd_end=self.scheduler.end,
                                             seed=idx, fclear=fclear)
        self.scheduler.initdb(designbase=plan)
        self.field_ra = self.scheduler.fields.racen
        self.field_dec = self.scheduler.fields.deccen
        cadencelist = self.scheduler.fields.cadencelist.cadences
        cadences = self.scheduler.fields.cadence

        self.nom_duration = np.float32(15. / 60. / 24.)
        self.cals = np.float32(3. / 60. / 24.)
        self.observe = observesim.observe.Observe(defaultExp=self.nom_duration, 
                                             cadencelist=cadencelist, cadences=cadences)

        self.curr_mjd = np.float32(1e9)


    def nextField(self, mjd_morning_twilight):
        # should we do this now?
        isclear, nextchange_weather = self.weather.clear(mjd=self.curr_mjd)
        onoff, nextchange_on = self.scheduler.on(mjd=self.curr_mjd)
        nextchange_all = np.array([nextchange_weather, nextchange_on, mjd_morning_twilight])
        self.nextchange = np.min(nextchange_all)
        if not ((isclear == True) and (onoff == 'on')):
            if(self.nextchange > self.curr_mjd):
                self.curr_mjd = self.nextchange
            return -1, -1
        # dark time or brighttime? to guess at how long we need for obs
        skybrightness = self.scheduler.skybrightness(self.curr_mjd)
        if skybrightness < 0.3:
            airmass_weight = 1.05
        else:
            airmass_weight = 0.05
        # integer division floors; no partial exposures
        maxExp = int((self.nextchange - self.curr_mjd)//(self.nom_duration * 1.3 ** airmass_weight))
        if maxExp == 0:
            self.curr_mjd = self.curr_mjd + self.nom_duration
            return -1, -1
        fieldid, nexposures = self.scheduler.nextfield(mjd=self.curr_mjd,
                                                  maxExp=maxExp)
        if(fieldid is not None):
            new_alt, new_az = self.scheduler.radec2altaz(mjd=self.curr_mjd, 
                                            ra=self.field_ra[fieldid],
                                            dec=self.field_dec[fieldid])
            new_duration = self.nom_duration * (1/np.cos(np.pi *\
                             (90-new_alt) / 180.)) ** airmass_weight
            final_alt, final_az = self.scheduler.radec2altaz(
                                    mjd=self.curr_mjd + new_duration * nexposures, 
                                    ra=self.field_ra[fieldid],
                                    dec=self.field_dec[fieldid])
            maxTime = self.nextchange - self.curr_mjd
            if  maxTime < new_duration * nexposures or \
                final_alt < 20:
                fieldids, nexps, priorities = self.scheduler.nextfield(mjd=self.curr_mjd,
                                                  maxExp=maxExp, returnAll=True)
                alts, azs = self.scheduler.radec2altaz(mjd=self.curr_mjd, 
                                            ra=self.field_ra[fieldids],
                                            dec=self.field_dec[fieldids])
                adj_exp = self.nom_duration * np.power((1/np.cos(np.pi * \
                                    (90 - alts) / 180.)), airmass_weight)
                fieldid, nexposures = sortFields(fieldids, nexps, priorities, adj_exp, maxTime=maxTime)
                if fieldid == -1:
                    print("baawaaaaaahhhahahaa :( ")
                    self.curr_mjd = self.curr_mjd + self.nom_duration/20
                    return -1, -1
            return fieldid, nexposures
        # else:
        #     lsts.append(self.scheduler.lst(self.curr_mjd)[0])
        #     observed.append(-1)
        #     self.curr_mjd = self.curr_mjd + duration
        #     exp_tonight += duration


    def observeField(self, fieldid, nexposures):

        slewtime = np.float32(2. / 60. / 24.) # times some function of angular distance?

        self.curr_mjd = self.curr_mjd + self.cals + slewtime

        for i in range(nexposures):
            # add each exposure
            alt, az = self.scheduler.radec2altaz(mjd=self.curr_mjd, 
                                        ra=self.field_ra[fieldid],
                                        dec=self.field_dec[fieldid])
            airmass = 1/np.cos(np.pi * (90-alt) / 180.)
            # observed.append(self.scheduler.fields.racen[fieldid])

            result = self.observe.result(mjd=self.curr_mjd, fieldid=fieldid,
                                    airmass=airmass)
            duration = result["duration"]
            if duration < 0:
                print("HOOOWWWOWOWOWOWW")
                print(alt, az, self.curr_mjd, fieldid)

            self.curr_mjd = self.curr_mjd + duration
            # lsts.append(self.scheduler.lst(self.curr_mjd)[0])

            if(self.curr_mjd > self.nextchange):
                oops = (self.curr_mjd - self.nextchange) * 24 * 60
                if oops > 5:
                    print("NOOOO! BAD!", oops)
                    print(i, nexposures, alt)
            
            self.scheduler.update(fieldid=fieldid, result=result)

            # exp_tonight += duration


    def observeMJD(self, mjd):
        # uncomment to do a quick check
        # if mjd > 59200:
        #     continue
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
            fieldid, nexposures = self.nextField(mjd_morning_twilight)
            self.observeField(fieldid, nexposures)

            
        # weather_used[mjd] = {"length": night_len, "observed": exp_tonight}
