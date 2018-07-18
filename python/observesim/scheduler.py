import os
import numpy as np
import scipy.optimize as optimize
import PyAstronomy.pyasl as pyasl
import astropy.units as units
import astropy.time as atime
import pydl.pydlutils.yanny as yanny
import observesim.fields
import observesim.observations
import observesim.cadence
from observesim.moonphase import moonphase2
from observesim.sunpos2 import sunpos2

"""Scheduler module class.

Dependencies:

 numpy
 scipy
 PyAstronomy
 astropy
 pydl

"""


def dateandtime2mjd(date=None, time='12:00', to_tai=7):
    """Utility to calculate an MJD"""
    if((type(date) is list) | (type(date) is np.ndarray)):
        isotimes = ["{date} {time}".format(date=cdate, time=ctime)
                    for cdate, ctime in zip(date, time)]
    else:
        isotimes = "{date} {time}".format(date=date, time=time)
    times = atime.Time(isotimes, format='iso', scale='tai')
    times = times + np.int32(to_tai) * units.hour
    return(times.mjd)


class SchedulerBase(object):
    """Scheduler base class with generic utilities.

    Parameters:
    ----------

    Attributes:
    ----------

    Methods:
    -------

    ralst2ha(ra=, lst=) : convert RA and LST to hour angle
    hadec2altaz(ha=, dec=, lat=) : convert HA, Dec, latitude to alt, az
    alt2airmass(alt=) : convert altitude to airmass

"""
    def __init__(self):
        return

    def _arrayify(self, quantity=None):
        """Cast quantity as ndarray of numpy.float64"""
        try:
            length = len(quantity)
        except TypeError:
            length = 1
        return np.zeros(length, dtype=np.float64) + quantity

    def _mjd2jd(self, mjd=None):
        """Convert MJD to JD"""
        return (self._arrayify(mjd) + np.float64(2400000.5))

    def ralst2ha(self, ra=None, lst=None):
        """Return HA (degrees) given RA and LST

        Parameters:
        ----------

        ra : np.float64
            right ascension (deg)

        lst : np.float64
            local sidereal time (deg)

        Returns:
        -------

        ha : np.float64
            hour angle (deg)

        """
        ha = (((self._arrayify(lst) - self._arrayify(ra) + 360. + 180.)
               % 360.) - 180.)
        return(ha)

    def hadec2altaz(self, ha=None, dec=None, lat=None):
        """Return (alt, az) (degrees) of given HA and Dec and latitude

        Parameters:
        ----------

        ha : np.float64
            hour angle (deg)

        dec : np.float64
            declination (deg)

        lat : np.float64
            latitude (deg)

        Returns:
        -------

        alt : np.float64
            altitude (deg)

        az : np.float64
            azimuth (deg E of N)
        """
        aha = self._arrayify(ha)
        adec = self._arrayify(dec)
        (alt, az) = pyasl.hadec2altaz(aha, adec,
                                      np.float64(lat) + np.zeros(len(aha)))
        return (alt, az)

    def alt2airmass(self, alt):
        """Return airmass given altitude

        Parameters:
        ----------

        alt : np.float64
            altitude (deg)

        Returns:
        -------

        airmass : np.float64
            airmass (1/sin(altitude))

        """
        airmass = 1. / np.sin(np.pi / 180. * self._arrayify(alt))
        return(airmass)


class Observer(SchedulerBase):
    """Observer class to define different observatories.

    Parameters:
    ----------

    observatory : str
        Name of observatory to use (must be in observatory file)
        (default 'apo')

    observatoryfile : str
        Name of Yanny-format observatory file to read
        (default $OBSERVESIM_DIR/data/observatories.par)

    Attributes:
    ----------

    observatory : str
        Name of observatory

    latitude : numpy.float64
        Latitude of observatory

    longitude : numpy.float64
        Longitude (E of Greenwich) of observatory

    Methods:
    -------

    ralst2ha(ra=, lst=) : convert RA and LST to hour angle
    hadec2altaz(ha=, dec=, lat=) : convert HA, Dec, latitude to alt, az
    alt2airmass(alt=) : convert altitude to airmass
    lst(mjd=) : return LST in degrees for observer at given MJD (days)
    radec2altaz(mjd=, ra=, dec=) : return alt/az for ra/dec at given MJD
    sun_radec(mjd=) : return position of Sun in Ra/Dec
    sun_altaz(mjd=) : return position of Sun in Alt/AZ
    moon_radec(mjd=) : return position of Moon in Ra/Dec
    moon_altaz(mjd=) : return position of Moon in Alt/AZ
    moon_illumination(mjd=) : return illumination of Moon at given MJD
    evening_twilight(mjd=): return evening twilight on MJD
    morning_twilight(mjd=): return morning twilight on MJD
"""
    def __init__(self, observatory='apo', observatoryfile=None):
        """Create Observer object"""
        super().__init__()
        self.observatory = observatory
        if(observatoryfile is None):
            observatoryfile = os.path.join(os.getenv('OBSERVESIM_DIR'),
                                           'data', 'observatories.par')
        self._file = observatoryfile
        self._data = yanny.yanny(self._file)
        observatories = np.array([observatory.decode()
                                  for observatory in
                                  self._data['OBSERVATORY']['observatory']])
        indx = np.where(observatories == self.observatory)[0][0]
        self.latitude = self._data['OBSERVATORY']['latitude'][indx]
        self.longitude = self._data['OBSERVATORY']['longitude'][indx]

    def lst(self, mjd=None):
        """Return LST (degrees) given MJD for observer

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        lst : np.float64
            local sidereal time (deg)
        """
        mjds = self._arrayify(mjd)
        lst = (np.float64(15.) *
               pyasl.ct2lst(self._mjd2jd(mjds),
                            np.zeros(len(mjds)) + self.longitude))
        return (lst)

    def sun_radec(self, mjd=None):
        """Return (ra, dec) in deg J2000 for Sun at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        ra : np.float64
            right ascension, J2000 (deg)

        dec : np.float64
            declination, J2000 (deg)
        """
        jd = self._mjd2jd(mjd=self._arrayify(mjd))
        (tmp_jd, ra, dec) = sunpos2(jd)
        return (ra, dec)

    def moon_radec(self, mjd=None):
        """Return (ra, dec) in deg J2000 for Moon at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        ra : np.float64
            right ascension, J2000 (deg)

        dec : np.float64
            declination, J2000 (deg)
        """
        jd = self._mjd2jd(mjd=self._arrayify(mjd))
        ra, dec, dist, geolon, geolat = pyasl.moonpos(jd)
        return (ra, dec)

    def radec2altaz(self, mjd=None, ra=None, dec=None):
        """Return (alt, az) for (ra, dec) in deg J2000 at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        ra : np.float64
            right ascension, J2000 (deg)

        dec : np.float64
            declination, J2000 (deg)

        Returns:
        -------

        alt : np.float64
            altitude (deg)

        az : np.float64
            azimuth (deg E of N)
        """
        lst = self.lst(mjd=mjd)
        ha = self.ralst2ha(ra=ra, lst=lst)
        (alt, az) = self.hadec2altaz(ha=ha, dec=dec, lat=self.latitude)
        return (alt, az)

    def moon_illumination(self, mjd=None):
        """Return Moon illumination at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        illumination : np.float64
            fraction of Moon illuminated
        """
        jd = self._mjd2jd(mjd=self._arrayify(mjd))
        return (moonphase2(jd))

    def sun_altaz(self, mjd=None):
        """Return (alt, az) for Sun at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        alt : np.float64
            altitude (deg)

        az : np.float64
            azimuth (deg E of N)
        """
        (ra, dec) = self.sun_radec(mjd=mjd)
        (alt, az) = self.radec2altaz(mjd=mjd, ra=ra, dec=dec)
        return (alt, az)

    def moon_altaz(self, mjd=None):
        """Return (alt, az) for Moon at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        alt : np.float64
            altitude (deg)

        az : np.float64
            azimuth (deg E of N)
        """
        (ra, dec) = self.moon_radec(mjd=mjd)
        (alt, az) = self.radec2altaz(mjd=mjd, ra=ra, dec=dec)
        return (alt, az)

    def lunation(self, mjd=None):
        """Return Moon illumination, or zero if Moon at alt<0"""
        (moon_alt, moon_az) = self.moon_altaz(mjd=mjd)
        if(moon_alt < 0):
            return(0.)
        else:
            return(self.moon_illumination(mjd=mjd))

    def _twilight_function(self, mjd=None, twilight=-8.):
        """Utility function for root-finding to get twilight times"""
        (alt, az) = self.sun_altaz(mjd=mjd)
        return (alt - twilight)

    def evening_twilight(self, mjd=None, twilight=-8.):
        """Return MJD (days) of evening twilight for MJD

        Parameters:
        ----------

        mjd : np.int32, int
            Modified Julian Day (days)

        Returns:
        -------

        evening_twilight : np.float64
            time of twilight in MJD (days)
        """
        if(np.floor(np.float64(mjd)) != np.float64(mjd)):
            raise ValueError("MJD should be an integer")
        noon_ish = (np.float64(mjd) -
                    self.longitude / 15. / 24. - 0.5)
        midnight_ish = noon_ish + 0.5
        twi = optimize.brenth(self._twilight_function,
                              noon_ish, midnight_ish,
                              args=(twilight))
        return(np.float64(twi))

    def morning_twilight(self, mjd=None, twilight=-8.):
        """Return MJD (days) of morning twilight for MJD

        Parameters:
        ----------

        mjd : np.int32, int
            Modified Julian Day (days)

        Returns:
        -------

        morning_twilight : np.float64
            time of twilight in MJD (days)
        """
        if(np.floor(np.float64(mjd)) != np.float64(mjd)):
            raise ValueError("MJD should be an integer")
        midnight_ish = (np.float64(mjd) -
                        self.longitude / 15. / 24.)
        nextnoon_ish = midnight_ish + 0.5
        twi = optimize.brenth(self._twilight_function,
                              midnight_ish, nextnoon_ish,
                              args=(twilight))
        return(np.float64(twi))


class Master(Observer):
    """Master class to interpret master schedule as an observer

    Parameters:
    ----------
    schedulefile : str
        schedule file to use; default $OBSERVESIM_DIR/data/master_schedule.par

    Attributes:
    ----------
    start : np.int32
        MJD (days) of first night of survey
    end : np.int32
        MJD (days) of last night of survey
    mjds : ndarray of np.int32
        MJDs (days) when survey is potentially active
    events : ndarray of numpy.str_
        names of events of note
    event_dates : ndarray of numpy.str_
        list of dates in ISO format for events of note
    event_times : ndarray of numpy.str_
        list of times in ISO format for events of note
    event_mjd : ndarray of numpy.float64
        MJDs (days) of events of note

    Methods:
    -------
    on() : is the survey on
"""
    def __init__(self, schedulefile=None, observatory='apo',
                 observatoryfile=None):
        """Create Master object for schedule"""
        super().__init__(observatory=observatory,
                         observatoryfile=observatoryfile)
        if(schedulefile is None):
            masterfile = 'master_schedule_{observatory}.par'.format(observatory=observatory)
            schedulefile = os.path.join(os.getenv('OBSERVESIM_DIR'),
                                        'data', masterfile)
        self._schedulefile = schedulefile
        self.schedule = yanny.yanny(self._schedulefile)
        self._validate()
        self.event_dates = np.array([date.decode() for date
                                     in self.schedule['SCHEDULE']['date']])
        self.event_times = np.array([time.decode() for time
                                     in self.schedule['SCHEDULE']['time']])
        self.event_mjds = self._dateandtime2mjd()
        self.events = np.array([event.decode() for event
                                in self.schedule['SCHEDULE']['event']])
        self.start = self._start()
        self.end = self._end()
        self.mjds = self._mjds()
        return

    def _dateandtime2mjd(self):
        return(dateandtime2mjd(date=self.event_dates,
                               time=self.event_times,
                               to_tai=self.schedule['to_tai']))

    def _validate(self):
        # should make sure:
        #  one start (first event)
        #  one end (last event)
        #  start MJD is a daytime time 
        #  START_SURVEY is "on" 
        #  END_SURVEY is "off" 
        return

    def on(self, mjd=None):
        if(mjd < self.event_mjds[0]):
            return('off', self.event_mjds[0])
        if(mjd >= self.event_mjds[-1]):
            return('off', mjd + 1.)
        # Assumes there is only one
        indx = np.where((mjd >= self.event_mjds[0:-1]) &
                        (mjd < self.event_mjds[1:]))[0][0]
        return(self.schedule[self.events[indx]],
               self.event_mjds[indx + 1])

    def end_mjd(self):
        """Return end MJD

        Returns:

        end_mjd : np.float64
            MJD of last event (end of survey)
        """
        return(self.event_mjds[-1])

    def _start(self):
        # Assumes there is only one
        indx = np.where(self.events == 'START_SURVEY')[0][0]
        # Assumes START_SURVEY turns state on
        return(np.int32(np.floor(self.event_mjds[indx])))

    def _end(self):
        # Assumes there is only one
        indx = np.where(self.events == 'END_SURVEY')[0][0]
        # Assumes END_SURVEY turns state off
        return(np.int32(np.ceil(self.event_mjds[indx])))

    def _mjds(self):
        nmjd = self.end - self.start + 1
        mjds = self.start + np.arange(nmjd, dtype=np.int32)
        keep = np.zeros(nmjd, dtype=np.int32)
        for indx in np.arange(len(self.events) - 1):
            this_event = self.events[indx]
            if(self.schedule[this_event] == 'on'):
                keep_start = np.int32(np.floor(self.event_mjds[indx]))
                keep_end = np.int32(np.ceil(self.event_mjds[indx + 1]))
                ikeep = np.where((mjds >= keep_start) &
                                 (mjds <= keep_end))[0]
                keep[ikeep] = 1
        ikeep = np.where(keep)[0]
        return(mjds[ikeep])


class Scheduler(Master):
    """Scheduler class.

    Parameters:
    ----------

    airmass_limit : float, np.float32
        airmass limit for observations

    Attributes:
    ----------

    airmass_limit : float, np.float32
        airmass limit for observations

    master : Master object
        Master schedule to use for scheduling

    observer : Observer object
        Observer to use for scheduling

    fields : Fields object
        object for fields

    observations : Observations object
        object accessing list of observations

    Methods:
    -------

    initdb() : initialize field list and set to unobserved

    nextfield(mjd=mjd) : return field to observe at mjd

    observable(mjd=mjd) : return fieldids observable at mjd

    update(fieldid=fieldid, result=result) : update observations with result

    Comments:
    --------

    Scheduling proceeds conceptually as follows
         - fields are limited to set that are conceivably observable
         - A strategy to optimize completion

    In this default Scheduler, the strategy is a completely heuristic one
         - take lowest HA cases in bins of 5 deg
         - take lowest transit altitude case among those

    """
    def __init__(self, airmass_limit=2.,
                 schedulefile=None, observatory='apo', observatoryfile=None):
        """Return Scheduler object
        """
        super().__init__(schedulefile=schedulefile, observatory=observatory,
                         observatoryfile=observatoryfile)
        self.airmass_limit = airmass_limit
        return

    def initdb(self, designbase='plan-0'):
        """Initialize Scheduler fields and observation lists
        """
        filebase = os.path.join(os.getenv('OBSERVING_PLAN_DIR'),
                                designbase)
        self.cadencelist = observesim.cadence.CadenceList()
        self.cadencelist.fromfits(filename=filebase + "-cadences.fits")
        self.fields = observesim.fields.Fields()
        self.fields.fromfits(filename=filebase + "-fields-" +
                             self.observatory + ".fits")
        self.observations = observesim.observations.Observations(observatory=self.observatory)
        return

    def observable(self, mjd=None, check_lunation=True, check_cadence=True):
        """Return array of fields observable

        Parameters:
        ----------

        mjd : np.float64
            current MJD
"""
        (alt, az) = self.radec2altaz(mjd=mjd, ra=self.fields.racen,
                                     dec=self.fields.deccen)
        airmass = self.alt2airmass(alt)
        lunation = self.lunation(mjd)
        observable = (alt > 0.) & (airmass < self.airmass_limit)
        iobs = np.where(observable)[0]
        if(check_cadence):
            indxs = np.where(self.fields.nextmjd > mjd)[0]
            observable[indxs] = False
            iob = np.where(observable)[0]
            #if(len(iob) == 0):
                #fitsio.write('test.fits', self.fields.nextmjd, clobber=True)
                #fitsio.write('test.fits', np.int32((alt > 0.) & (airmass < self.airmass_limit)))
                #sys,exit(1)
            indxs = np.where(self.fields.nextmjd <= mjd)[0]
            for indx in indxs:
                if(observable[indx]):
                    cadence = self.cadencelist.cadences[self.fields.cadence[indx]]
                    iobservations = self.fields.observations[indx]
                    mjd_past = self.observations.mjd[iobservations]
                    observable[indx] = cadence.evaluate_next(mjd_past=mjd_past,
                                                             mjd_next=mjd,
                                                             lunation_next=lunation,
                                                             check_lunation=check_lunation)

        iobservable = np.where(observable)[0]
        return self.fields.fieldid[iobservable]

    def pick(self, mjd=None, fieldid=None):
        """Return the fieldid to pick from using heuristic strategy

        Parameters:
        ----------

        fieldid : ndarray  of np.int32
            array of available fieldid values

        Returns:
        -------

        pick_fieldid : ndarray of np.int32
            fieldid
"""
        lst = self.lst(mjd)
        ha = self.ralst2ha(ra=self.fields.racen[fieldid], lst=lst)
        iha = np.int32(np.floor(np.abs(ha) / 5.))
        minha = iha.min()
        iminha = np.where(iha == minha)[0]
        dec = self.fields.deccen[fieldid[iminha]]
        if(self.observatory == 'apo'):
            ipick = np.argmin(dec)
        else:
            ipick = np.argmax(dec)
        pick_fieldid = fieldid[iminha[ipick]]
        return(pick_fieldid)

    def nextfield(self, mjd=None):
        """Picks the next field to observe

        Parameters:
        ----------

        mjd : np.float64
            Current MJD (days)

        Returns:
        --------

        fieldid : np.int32, int
            ID of field to observe
"""
        observable_fieldid = self.observable(mjd=mjd)
        if(len(observable_fieldid) == 0):
            # print("Nothing observable")
            observable_fieldid = self.observable(mjd=mjd,
                                                 check_cadence=False)
        fieldid = self.pick(fieldid=observable_fieldid, mjd=mjd)
        return(fieldid)

    def update(self, fieldid=None, result=None):
        """Update the observation list with result of observations

        Parameters:
        -----------

        fieldid : np.int32, int
            ID of field

        result : ndarray
            One element, contains 'mjd', 'duration', 'sn2'

        Comments:
        ---------

"""
        (alt, az) = self.radec2altaz(mjd=result['mjd'],
                                     ra=self.fields.racen[fieldid],
                                     dec=self.fields.deccen[fieldid])
        airmass = self.alt2airmass(alt)
        lunation = self.lunation(result['mjd'])
        lst = self.lst(result['mjd'])
        iobs = self.observations.add(fieldid=fieldid,
                                     mjd=result['mjd'],
                                     duration=result['duration'],
                                     sn2=result['sn2'],
                                     lunation=lunation,
                                     airmass=airmass,
                                     lst=lst)
        self.fields.add_observations(result['mjd'], fieldid, iobs)
        return
