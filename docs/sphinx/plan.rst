.. role:: header_no_toc
  :class: class_header_no_toc

.. title:: Plan for scheduling problem

What we need from the targeting teams
=====================================

* List of targets, and along with each one a number of epochs, exposures per epoch, and the name of a cadence type.

* Geometrical description of coverage of targets, in the form of several nside=64 healpix maps:
 * window: 1 indicates a region with targets and 0 indicates a region without.
 * epochs: maximum number of epochs for objects in this.
 * exposures: number of exposures per epoch

* Any list of specific field centers, if applicable.

Tiling
======

* We will define a set of appropriate field centers with the desired coverage of the sky.

* We will define a set of tiles covering those field centers that cover the region

* We will choose the 
	
Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
