.. role:: header_no_toc
  :class: class_header_no_toc

.. title:: Observing Simulation for SDSS-V

Observing Simulation for SDSS-V
===============================

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

Future Anticipated Changes in Design
------------------------------------

The Scheduler class should at some later point be refactored into its
own software product. The Fields and Observations should at some later
point also be refactored into either one or two different products,
and in survey operations need to be accessing databases rather than
holding results in memory. The methods for Fields and Observations
will ideally not change much when that change is made.

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

Observer
--------

The :ref:`Observer <Observer>` class defines the properties of the
observer and has methods for checking ephemerides and observability
for the specific observer.

This class gets information about the observatory from the
$OBSERVESIM_DIR/data/observatories.par file. Example content of this
file is here:

.. code-block:: tcl

   typedef struct {
      char observatory[100];
      double longitude;
      double latitude;
   } OBSERVATORY;

   OBSERVATORY apo -105.82027778 32.7797556


Master Schedule
---------------

Scheduler
---------

Fields
------

Observations
------------

Weather
-------

Reference
---------

- :ref:`What's new in observesim? <changelog>`

.. toctree::
   :maxdepth: 2

   api


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
