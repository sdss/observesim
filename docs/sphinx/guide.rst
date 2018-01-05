.. role:: header_no_toc
  :class: class_header_no_toc

.. title:: Guide for Observing Simulations for SDSS-V

Guide for Observing Simulations for SDSS-V
==========================================

Introduction
------------

This Python product simulates the SDSS-V observations.

It is designed specifically with the fiber robot observations in mind
but it contains some generic utilities for simulating any
observational situation.

The simulation aspects of this product consists of utilities to
simulate weather and observational results. The Weather and Observe
classes handle these functions.

A set of nested classes handle the scheduling aspects of the
simulations. These classes are Observer (which has methods returning
results associated with a specific observer on the Earth), its
subclass Master (which adds the master schedule), and its subclass
Scheduler (which finally adds the decision-making functionality). The
Scheduler uses a simple algorithm for scheduling and is meant as a
base class for more sophisticated approaches.

A pair of classes handles tracking of the fields and
observations. Fields contains methods that will search and return
information on the list of fields. Observations contains methods that
will search and return information on the observations taken so
far. These classes hold the information internally in memory. 

The simulation software depends on Scheduler, Weather, and
Observe. Scheduler depends on Fields and Observations (in fact, it
contains an object of each class as an attribute.

Please see the :doc:`API for reference </api.rst>`.


Future Anticipated Changes in Design
------------------------------------

The Scheduler class should at some later point be refactored into its
own software product.

The Fields and Observations classes currently allow access to
properties of the fields and observations simply by looking at ndarray
attributes (like racen, deccen). This access should be abstracted so
this access by the Scheduler is not dependent on the internal storage
mechanism for the information.

The Fields class currently assumes a particular distribution on the
sky. There needs to be a process that generates the set of desired
fields. Whether this is part of the Fields class or a different
product altogether remains to be seen. 

The Fields and Observations should also at some later point be
refactored into either one or two different products, and in survey
operations need to be accessing databases rather than holding results
in memory. The methods for Fields and Observations will ideally not
change much when that change is made.

The "results" passed around in the scheduler and the simulator is a
simple dictionary. I'm not sure how this should evolve so was keeping
it simple for now.

None of this yet incorporates any functionality associated with
assigning targets to fields. 

Running a Simulation
--------------------

To run a specific simulation, call from the Unix command line:

.. code-block:: sh

   sdss5_simulate -o <output base> [-s <seed>] 
 
This will run through the entire SDSS-V schedule, simulate the
observations, and output two files:

.. code-block:: sh

   <output base>-fields.fits
   <output base>-observations.fits

which contain the list of fields used and the list of observations
made and their results. 

Perusing the sdss5_simulate code will inform one's understanding of
the classes described below. 

Observer
--------

The :ref:`Observer <Observer>` class defines the properties of the
observer and has methods for checking ephemerides and observability
for the specific observer.

This class gets information about the observatory from the Yanny file:

::

     $OBSERVESIM_DIR/data/observatories.par

Example content of this file is here:

.. code-block:: tcl

   typedef struct {
      char observatory[100];
      double longitude;
      double latitude;
   } OBSERVATORY;

   OBSERVATORY apo -105.82027778 32.7797556

This file can be read with tools in the pydl.pydlutils.yanny module.

Master Schedule
---------------

The :ref:`Master <Master>` class is a subclass of :ref:`Observer
<Observer>` that incorporates information about the master schedule
for the observatory in question within SDSS-V.

::

    import observesim.scheduler as scheduler

    schedule = scheduler.Master()
    template = "mjd={mjd}, illumination={illumination}"
    for mjd in schedule.mjds:
        illumination = schedule.moon_illumination(mjd=mjd)
        print(template.format(mjd=mjd, illumination=illumination))

It also contains the method "on()", which returns whether the survey
is on at any given (floating point) MJD, and the next change in on/off
status. 

The master schedule itself is kept as a Yanny file at:

::

     $OBSERVESIM_DIR/data/master_schedule.par

This file can be read with tools in the pydl.pydlutils.yanny module.

Example content of the master schedule file is here:

::

    # Timezone offset in hours to apply to get to TAI
    # (i.e. Greenwich time)
    to_tai 7  # Mountain Standard Time

    # Whether events start ("on") or stop ("off") observing
    START_SURVEY on
    END_SURVEY off
    START_SHUTDOWN off
    END_SHUTDOWN on

    typedef enum {
      START_SURVEY,
      END_SURVEY,
      START_SHUTDOWN,
      END_SHUTDOWN
    } EVENT;

    typedef struct {
      char date[10];
      char time[5];
      EVENT event;
    } SCHEDULE;

    SCHEDULE 2020-07-01 12:00 START_SURVEY
    SCHEDULE 2020-07-10 12:00 START_SHUTDOWN
    SCHEDULE 2020-08-20 12:00 END_SHUTDOWN
    SCHEDULE 2021-07-01 12:00 END_SURVEY

Scheduler
---------

The :ref:`Scheduler <Scheduler>` class is a subclass of :ref:`Master
<Master>` that has tools necessary to schedule observations.

First, it contains attributes "fields" and "observations" which are
expected to be objects of the Fields and Observations classes
described below. These are used to access and update information about
the fields and observations.

It then provides the method "field()" which returns the next field to
observe given a (floating point) MJD. 

Finally, it provides the method "update()" which updates the status of
a field given a result. The result is currently assumed to be a
dictionary with 'sn2', 'mjd', and 'duration', as described below in
the simulation section. 

Fields
------

The :ref:`Fields <Fields>` class creates and stores the information
about the survey fields.

In its current form the Fields class initializes itself from the
Sloane distribution of tiles on the sky.

In its current form, one accesses the fields from a set of ndarrays
containing the RA, Dec, type of field, priority, etc. The fields are
indexed by fieldid, which is just the zero-indexed position of each
field in the arrays.

Observations
------------

The :ref:`Observations <Observations>` class creates and stores the
information about the observations.

In its current form, one access observations from the Observations
class through its attributes, which are just ndarrays of durations,
MJDs, and signal-to-noise squared.  The "add()" method adds an
observations. The "forfield()" method returns all observations for a
specific field.

Weather
-------

The :ref:`Weather <Weather>` class handles simulated weather
conditions. 

It is initialized with a starting and ending MJD, an optional random
seed, and optional other parameters including the fraction of clear
weather.

It has one method "clear()" which returns whether a specific
floating-point MJD is clear, and what the MJD of the next change in
state is.

The weather patterns have a power law frequency spectrum cut off at
high frequency with a Gaussian. This characteristic scale is by
default two days (but this leaves plenty of variation below that
scale). It does not characterize the weather beyond clear or not
clear.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
