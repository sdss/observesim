Observing Simulation for SDSS-V
==============================

Run an SDSS V simulation based on a robostrategy plan. 

Example usage:
- Suppose `~/data/simulations` points towards robostrategy outputs. To run a single simulation on the beta-1 RS plan (for APO, the default):  `sdss5_simulate -b ~/data/simulations/ -p beta-1`
    
- To run a suite of 10 simulations, using 4 processes (via a very simple implementation of python multiprocessing): `batch_sim -b ~/data/simulations/ -p beta-1 -n 10 -m 4`



Full documentation can be found at http://observesim.readthedocs.io/

| |Build Status|
| |Coverage Status|

------------

.. |Build Status| image:: https://travis-ci.org/blanton144/observesim.svg?branch=master
   :target: https://travis-ci.org/blanton144/observesim

.. |Build Status| image:: https://travis-ci.org/blanton144/observesim.svg?branch=master
   :target: https://travis-ci.org/blanton144/observesim

.. |Coverage Status| image:: https://coveralls.io/repos/github/blanton144/observesim/badge.svg?branch=master
   :target: https://coveralls.io/github/blanton144/observesim?branch=master
