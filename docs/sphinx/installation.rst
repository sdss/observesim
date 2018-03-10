.. role:: header_no_toc
  :class: class_header_no_toc

.. title:: Installing observesim

.. _observesim-installation:
	
Installation
============

observesim has been tested on Mac OS X and Linux operating systems. 

**Installation from Source**

In bash (change appropriately for tcsh or other shell):

    export OBSERVESIM_DIR=[your desired location]
    git clone https://github.com/sdss/observesim $OBSERVESIM_DIR
    export PYTHONPATH=$OBSERVESIM_DIR/python:$PYTHONPATH
    export PATH=$OBSERVESIM_DIR/bin:$PATH

In tcsh (change appropriately for tcsh or other shell):

    setenv OBSERVESIM_DIR [your desired location]
    git clone https://github.com/sdss/observesim $OBSERVESIM_DIR
    setenv PYTHONPATH $OBSERVESIM_DIR/python:$PYTHONPATH
    setenv PATH $OBSERVESIM_DIR/bin:$PATH

**Modules**

There is a Modules file in:

    $OBSERVESIM_DIR/modules/observesim.module

that can be used to build this into a Modules distribution, which is
standard for SDSS software. 

**SDSS install**

observesim has not been tested under sdssinstall yet.

**pip**

Although observesim is designed to be pip-installable, it has not been
registered, so is not yet.

.. _observesim-dependencies:

Dependencies
============

 * numpy
 * scipy
 * astropy
 * PyAstronomy
 * pydl
 * pyyaml
 * pygments
 * click
 * fitsio
