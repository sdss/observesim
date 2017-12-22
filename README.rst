Observing Simulation for SDSS-V
==============================

This product is for observing simulation code for SDSS-V. This
location on GitHub is temporary before being transferred to the main
SDSS repositories.


Master Schedule
---------------

The ``master`` module allows the interpretation of the master schedule
and the calculation of basic astronomical parameters for each night
within it.

The master schedule itself is kept as a Yanny file at:

::

     $OBSERVESIM_DIR/data/master_schedule.par

and an example is given below. It contains a list of events labeled by
local time and date, each of which either turns the survey on or turns
it off.

The list of observatories is at

::

     $OBSERVESIM_DIR/data/observatories.par

and an example is given below.

An example of how to use the master schedule module is as follows:

::

    import observesim.master as master

    schedule = master.Master()
    template = "mjd={mjd}, illumination={illumination}"
    for mjd in schedule.mjds:
        night = master.Night(mjd=mjd)
        illumination = night.moon.illumination(mjd=mjd)
        print(template.format(mjd=mjd, illumination=illumination))

The ``Night`` class also has attributes which are the floating point
MJDs of evening and morning twilight. Thus, a simulation can use this
module to loop through all of the nights, and then simulate each night.
The ``sun`` and ``moon`` attributes of ``Night``, and the ``lst()``
method, will give appropriate information for planning the night.

Simulation Code
---------------

The ``sdss5_simulate`` executable performs an actual simulation. You
must specify ``-o [outbase]``, and then the output goes in:

::

     [outbase]-fields.fits
     [outbase]-observations.fits

which give the list of fields and observations. You can also set a
specific random seed with ``-s [integer]``.

The first critical object that this code uses is the ``Scheduler``
object, which contains the algorithms for picking the next field, and
can be used to access the current state of observations. Specifically,
it contains ``Fields``, ``Observations``, ``Master``, and ``Observer``
objects. I don't know that I've settled on the exact best factoring of
the code here.

The second critical object that is used is the ``Weather`` object,
which returns whether at any given MJD it is clear, and when the next
change of state of the weather will be. 

Dependencies:

::

     numpy
     scipy
     astropy
     pydl
     PyAstronomy

The reason for PyAstronomy is that it has Meeus-based routines for
astronomical calculations that are substantially faster than the more
accurate routines found in astropy. There are a couple of problems in it
with the Sun position calculation for large vectors. This part of our
code should be cleaned up

Example contents for master\_schedule.par
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Example contents for observatories.par
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    typedef struct {
      char observatory[100];
      double longitude;
      double latitude;
    } OBSERVATORY;

    OBSERVATORY APO -105.82027778 32.7797556


| |Build Status|
| |Coverage Status|

------------

.. |Build Status| image:: https://travis-ci.org/blanton144/observesim.svg?branch=master
   :target: https://travis-ci.org/blanton144/observesim

.. |Build Status| image:: https://travis-ci.org/blanton144/observesim.svg?branch=master
   :target: https://travis-ci.org/blanton144/observesim

.. |Coverage Status| image:: https://coveralls.io/repos/github/blanton144/observesim/badge.svg?branch=master
   :target: https://coveralls.io/github/blanton144/observesim?branch=master
