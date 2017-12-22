import numpy as np
import PyAstronomy.pyasl as pyasl
import observesim.master
import observesim.fields
import observesim.observations

"""Scheduler module class.

Dependencies:

 numpy
 scipy

"""


class Scheduler(object):
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
        object accessing list of fields

    observations : Observations object
        object accessing list of observations

    Methods:
    -------

    initdb() : initialize field list and set to unobserved
    field(mjd=mjd) : return field to observe at mjd
    observable(mjd=mjd) : return fieldids observable at mjd
    set_priority_all(mjd=mjd) : reset all priorities
    set_priority(fieldid=fieldid, mjd=mjd) : reset priority for one field
    update(fieldid=fieldid, result=result) : update observations with result

    Comments:
    --------

    Scheduling proceeds conceptually as follows:

      * fields are limited to set that are conceivably observable
      * highest priority fields to observe are selected among the
      * A strategy to optimize completion

    In this default Scheduler, the strategy is a completely heuristic one:
         - take lowest HA cases in bins of 5 deg
         - take lowest transit altitude case among those

    """
    def __init__(self, airmass_limit=2.):
        """Return Scheduler object

        Parameters:
        ----------

        observer : Observer object
            Observer object to use for calculations

        airmass_limit : float, np.float32
            airmass limit for observations (default 2)
        """
        self.master = observesim.master.Master()
        self.observer = observesim.master.Observer()
        self.airmass_limit = airmass_limit
        return

    def initdb(self):
        """Initialize Scheduler fields and observation lists
        """
        self.fields = observesim.fields.Fields(observatory=self.observer.observatory)
        self.observations = observesim.observations.Observations(observatory=self.observer.observatory)
        return

    def set_priority(self, fieldid=None, mjd=None):
        """Set priority for a single field

        Parameters:
        ----------

        fieldid : np.int32, int
            fieldid for field to set priority for

        mjd : np.float64
            current MJD (only use observations prior to MJD)

        Comments:
        --------

        Sets Scheduler.fields.priority for fieldid
        """
        observations = self.observations.forfield(fieldid=fieldid, mjd=mjd)
        if(self.fields.fieldtype[fieldid] == 'standard'):
            tsn2 = observations['sn2'].sum()
            if(tsn2 > 2500.):
                self.fields.priority[fieldid] = self.fields._limit
            else:
                self.fields.priority[fieldid] = 1
        else:
            print("Error.")
        return

    def set_priority_all(self, mjd=None):
        """Set status of all fields (typically for beginning of night)

        Parameters:
        ----------

        mjd : np.float64
            current MJD (only use observations prior to MJD)

        """
        for fieldid in self.fields.fieldid:
            self.set_priority(mjd=mjd, fieldid=fieldid)
        return

    def observable(self, mjd=None):
        """Return array of fields observable

        Parameters:
        ----------

        mjd : np.float64
            current MJD
        """
        lst = self.observer.lst(mjd)
        ha = (((lst - self.fields.racen) + 180. + 360.) % 360.) - 180.
        lat = self.observer.latitude + np.zeros(len(ha))
        (alt, az) = pyasl.hadec2altaz(ha, self.fields.deccen, lat)
        airmass = 1. / np.sin(np.pi / 180. * alt)
        iobservable = np.where((alt > 0.) &
                               (airmass < self.airmass_limit))[0]
        return self.fields.fieldid[iobservable]

    def highest_priority(self, fieldid=None):
        """Return the fieldids which are in the highest priority class

        Parameters:
        ----------

        fieldid : ndarray  of np.int32
            array of fieldid values

        Returns:
        -------

        highest_fieldid : ndarray of np.int32
            array of fieldid values in highest priority class
        """
        priority = self.fields.priority[fieldid]
        highest_priority = priority.min()
        ihighest_priority = np.where(priority == highest_priority)[0]
        priority_fieldid = fieldid[ihighest_priority]
        return(priority_fieldid)

    def pick(self, mjd=None, fieldid=None):
        """Return the fieldid to pick from using heuristic strategy

        Parameters:
        ----------

        fieldid : ndarray  of np.int32
            array of fieldid values

        Returns:
        -------

        pick_fieldid : ndarray of np.int32
            fieldids
        """
        lst = self.observer.lst(mjd)
        ha = (((lst - self.fields.racen[fieldid]) + 180. + 360.) % 360.) - 180.
        iha = np.int32(np.floor(np.abs(ha) / 5.))
        minha = iha.min()
        iminha = np.where(iha == minha)[0]
        dec = self.fields.deccen[fieldid[iminha]]
        ipick = np.argmin(dec)
        pick_fieldid = fieldid[iminha[ipick]]
        return(pick_fieldid)

    def field(self, mjd=None):
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
        highest_fieldid = self.highest_priority(fieldid=observable_fieldid)
        fieldid = self.pick(fieldid=highest_fieldid, mjd=mjd)
        return(fieldid)

    def update(self, fieldid=None, result=None):
        """Update the observation list with result of observations

        Parameters:
        ----------

        fieldid : np.int32, int
            ID of field

        result : ndarray
            One element, contains 'mjd', 'duration', 'sn2'

        Comments:
        --------

        Updates the Scheduler.observations object

        """
        self.observations.add(fieldid=fieldid,
                              mjd=result['mjd'],
                              duration=result['duration'],
                              sn2=result['sn2'])
        self.set_priority(fieldid=fieldid, mjd=result['mjd'])
        return
