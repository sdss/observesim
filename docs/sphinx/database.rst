.. role:: header_no_toc
  :class: class_header_no_toc

:tocdepth: 2


Using the database
==================

Schema
------

The data used for survey simulations is stored in the SDSS-V dedicated database ``sdss5db``, in particular in the ``targetdb`` schema. The following links give access to the :download:`graffle <../../schema/targetdb.graffle>` (:download:`PDF <../../schema/targetdb.pdf>`) and the :download:`SQL <../../schema/targetdb.sql>` files that implement the ``targetdb`` schema.

The current version of the schema is ``v0.3.3``. If you make modifications to the schema, please remember to update the version in both the graffle and SQL file.


What's in targetdb?
-------------------

``targetdb`` contains the tables and data strictly necessary for survey simulations (and, in the future, for on-sky survey scheduling). As such, it contains only the targets that *can* be scheduled, as opposed to the databases used for target selection, that contain the entirety of the parents catalogues used to select the sample. ``targetdb`` only contains the target information (e.g., cadence) necessary for scheduling, in order to maintain the size of the schema small enough for Utah-observatories syncing.

The main section of :download:`targetdb <../../schema/targetdb.pdf>` comprises the ``target`` table as well as multiple auxiliary tables (``target_cadence``, ``program``, ``target_type``, etc.) At this time all targets have ``target_type='Science'`` but in the future we will add standard, skies, and maybe guide star targets. Many targets have magnitude information (in a variety of bands) as well as stellar parameters. While it is unlikely that we will need that data for scheduling they can be used for a number of QA and sample checks. Refer to :ref:`loading-targetdb` for details on the data currently loaded.

Another section of ``targetdb`` contains tables with information about the fibres and actuators, according to the `latest focal plane system layout <https://internal.sdss.org/trac/as4/wiki/FPSLayout>`_. A ``tile`` table links position on the sky with targets and robot configuration.

A final group of tables stores information about exposures and spectra. In the operations database, those tables are likely to be moved to their own schema. The ``simulation`` table connects with ``tile`` and ``exposure`` and can be used to catalogue what exposures and tiles have been created for each simulation.

Several tables, such as ``weather`` and ``lunation`` are placeholder for future functionality.


Basic connection
----------------

The main copy of the ``sdss5db`` database lives at Utah and can be accessed on port 5432 of ``db.sdss.utah.edu``. The user ``sdssdb_admin`` has read/write access to the database. The following line should open a Postgresql terminal for ``sdss5db``

.. code-block:: none

    psql -h db.sdss.utah.edu -U sdssdb_admin -p 5432 sdss5db

The password for user ``sdssdb_admin`` is the same as for ``apodb`` and other SDSS databases. In general, it is recommended (and required for the :ref:`orm`) to store the database password in your ``~/.pgpass`` file (remember to set its permissions to 600). For instance

.. code-block:: none

    *:*:sdss5db:sdssdb_admin:<insert_password_here>

will give you access to ``sdss5db`` as the ``sdssdb_admin`` user in any host.


.. _tunnel:

Tunnelling
----------

It is possible to tunnel to the Postgresql server hosting the database. The easiest way is to SSH to a machine at Utah and forward the port 5432 in ``db.sdss.utah.edu`` to your local machine. For instance, this SSH configuration

.. code-block:: none

    Host utah
        User u0931042
        HostName manga.sdss.org
        LocalForward localhost:6666 db.sdss.utah.edu:5432
        ForwardAgent yes

will redirect the Postgresql server port 5432 at Utah to your port 6666 (change the user unit and the machine to which you prefer to connect). Then, doing

.. code-block:: none

    psql -U sdssdb_admin -p 6666 sdss5db

will log you to the DB at Utah. This assumes that Postgresql is installed in your computer.


Dumping and restoring
---------------------

If you prefer to get your own local copy of ``sdss5db`` you can dump the database to a SQL file and restore it in your machine. SSH to Utah and run

.. code-block:: none

     pg_dump sdss5db > sdss5db.sql

(you can add a ``-n targetdb`` if you only want to dump the ``targetdb`` schema).

Then, in your machine, create ``sdss5db`` if it does not exist

.. code-block:: none

    createdb -T template0 sdss5db

Log in to the new DB and set the correct permissions

.. code-block:: postgresql

    psql sdss5db
    CREATE USER sdssdb_admin WITH PASSWORD '<insert_password_here>';
    GRANT ALL ON DATABASE sdss5db to sdssdb_admin;
    GRANT SELECT ON ALL TABLES IN SCHEMA targetdb to sdssdb_admin;

(note that you only need to do the above steps once, even if you delete and recreate the database). Now exit Postgresql and restore the database by doing

.. code-block:: none

    psql -U sdssdb_admin sdss5db < sdss5db.sql


.. _orm:

ORM access to targetdb
----------------------

In general, you will want to interact with the database using `Object-relational mapping (ORM) <https://en.wikipedia.org/wiki/Object-relational_mapping>`__ classes and methods, such as those provided by `SQLAlchemy <http://www.sqlalchemy.org>`_ or `peewee <http://docs.peewee-orm.com/en/latest/>`_. observesim provides database connections and model classes for ``targetdb`` for both libraries.

observesim includes three connection profiles, ``local`` for a database served by the local Postgresql server, ``utah`` to connect to the DB server at ``db.sdss.utah.edu`` from a machine at Utah, and ``tunnel`` to use a connection as described in :ref:`tunnel`. The profiles are defined in the ``database`` section of the package configuration file, ``observesim.cfg``. The configuration can be overridden by creating a local configuration file, ``~/.observesim/observesim.cfg``, as explained `here <http://sdss-python-template.readthedocs.io/en/latest/#configuration-file-and-logging>`__.


Using peewee
~~~~~~~~~~~~

If you are using peewee 3, importing the model classes is quite simple. Just do ::

    from observesim.db.peewee import targetdb

which will try, in order, the local, Utah, and tunnel connections and connect to the first available one. A log info message will indicate the connection that is being used. If no connection is available, a null connection will be set. The model classes will still be available but no queries can be performed. You can use the `~observesim.db.peewee.DatabaseConnection` class to set a specific connection. For instance ::

    targetdb.database.connect_from_config('tunnel')

will change the connection to the tunnel profile. The model classes will, from that moment, query the corresponding tables via the tunnel connection. You can also pass custom configuration options ::

    targetdb.database.connect_from_parameters('mydb', user='me', password='1234')

Once the connection has been created, you can use the normal `peewee`_ syntax to interact with the database, e.g. ::

    first_target = targetdb.Target.select().first()
    print(first_target)
    >>> <Target: pk=6590939>
    print(first_target.target_type)
    >>> <TargetType: pk=0, label='Science'>


Using SQLAlchemy
~~~~~~~~~~~~~~~~

With `SQLAlchemy`_ you first need to import the connection you want to use, e.g. ::

    from observesim.db.sqlalchemy.connections import local as db

then import the model classes for targetdb ::

    from observesim.db.sqlalchemy.models import targetdb
    session = db.Session()
    first_target = session.query(targetdb.Target).first()
    print(first_target)
    >>> <Target: pk=6590939>
    print(first_target.target_type)
    >>> <TargetType: pk=0, label='Science'>

As usual, and in both the peewee and SQLAlchemy implementations, models have the same name as the table but changing snake style to camelcase (e.g., ``target_type`` --> ``TargetType``).


.. _loading-targetdb:

Loading targetdb
----------------

At this time, ``targetdb`` has been loaded using a series of (to a certain degree) mock catalogues for the different subprograms for MWM and BHM. The data is described in `this wiki page <https://internal.sdss.org/trac/as4/wiki/Targetlist>`__ and the files used for loading can be found `here <https://data.sdss.org/sas/mangawork/users/u0931042/sdss5_target_list/>`__. In particular, the file `target_list.dat <https://data.sdss.org/sas/mangawork/users/u0931042/sdss5_target_list/target_list.dat>`__ contains a table with the files to load and the mapping of column names for each one of them. This file can be used with `observesim loaddb` (see below) to ingest the data to targetdb.

The focal plane system layout (fibres and actuators) is described in a :download:`file <../../python/observesim/etc/fps_RTConfig.txt>` provided by Rick Pogge and expanded `here <https://internal.sdss.org/trac/as4/wiki/FPSLayout>`__.

The loading can be do using the CLI script ``observesim`` (located in the ``bin/`` directory) with the ``loaddb`` option. For instance

.. code-block:: none

    observesim -v loaddb --remove target_list.dat fps_RTConfig.txt

will load targets and actuators/fibres to ``targetdb`` using the information in files ``target_list.dat`` and ``fps_RTConfig.txt``, and will remove all previous information from the relevant tables. The command calls, in order, `~observesim.db.load.load_targets` and `~observesim.db.load.load_fibres` from the `observesim.db.load` module. Refer to ``observesim loaddb --help`` for information on additional options.


API Reference
-------------

.. automodapi:: observesim.db.peewee
   :no-inheritance-diagram:
   :headings: ~^

.. automodapi:: observesim.db.sqlalchemy.connections
   :no-inheritance-diagram:
   :headings: ~^

.. automodapi:: observesim.db.load
   :no-inheritance-diagram:
   :headings: ~^
