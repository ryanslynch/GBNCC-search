GBNCC-search: A complete python data reduction pipeline for the GBNCC
pulsar survey.

Authors: Kevin Stovall (Univeristy of New Mexico) and Ryan S. Lynch
(McGill University)

Repository Maintained by: Ryan S. Lynch (rlynch@physics.mcgill.ca)

This file last modified on: 2013-11-04

This pipeline is used for downloading data from the GBNCC data archive
(located at NRAO in Charlottesvile, VA), searching the data using a
standalone script or through the use of HPC clusters, and tracking
search jobs and verifying output.  Logistics are handled through
communication with a MySQL database that stores up-to-date information
about the raw data files and their status regarding searches.  The
pipeline currently supports usage on the Guillimin HPC cluster.

The main executables used for data reduction are:

download.py : Communicate with the MySQL database to download raw data
submit.py   : Submit jobs to a supported HPC cluster
search.py   : Search data with the use of an HPC cluster
trackjobs.py: Track the status of jobs in progress, verify results,
and update the MySQL database

In addition to standard python packages, the following packages are
required:

- PRESTO (www.cv.nrao.edu/~sransom/presto)
- MySQLdb 
- numpy
- scipy
- pylab
- matplotlib
- mpfit
- pyfits
- pyslalib

Before using this pipeline, you should edit pipeline/lib/python/config.py,
run setup.py, and update your PATH and PYTHONPATH.  Then find pulsars!

Contents: 

setup.py : Verify packages, test database connection, and set up directories
pipeline/ : Top level pipeline directory 
pipeline/bin/ : Executables
pipeline/bin/addObservations.py : Add raw data files to database
pipeline/bin/download.py : Download raw data
pipeline/bin/search.py : Search GBNCC data
pipeline/bin/submit.py : Submit jobs to a supported HPC cluster
pipeline/bin/trackjobs.py : Track HPC cluster jobs
pipeline/lib/ : Libraries
pipeline/lib/db_schema.txt : MySQL database description
pipeline/lib/for_tarball/ : Packages for Guillimin python tarball
pipeline/lib/GBNCC.zaplist : Standard Fourier zap list
pipeline/lib/make_python_tarball.py : Make tarball of python packages
pipeline/lib/python : Python modules
pipeline/lib/python/analyse_sp.py : Single pulse searching algorithm
pipeline/lib/python/clans2batch1.py : Script for running CLAS2 jobs
pipeline/lib/python/config.py : Pipeline configuration variables
pipeline/lib/python/database.py : Methods for communicating with MySQL 
database
pipeline/lib/python/diagnostics.py : Methods for calculating job 
diagnostics
pipeline/lib/python/group_sp.py : Methods for analyse_sp.py
pipeline/lib/python/handle_exit.py : Clean-up during unexpected exits
pipeline/lib/python/msub.py : Methods for communicating with moab que
pipeline/lib/python/profile_tools.py : Methods for analyzing pulsar profiles
pipeline/lib/python/ratings.py : Methods for calculating candidate ratings
pipeline/lib/python/sp_utils.py : Methods for analyse_sp.py
pipeline/lib/python/utils.py : General pipeline utilities
pipeline/lib/python-2.6.8-packages-base.tar.gz : Base tarball of common 
python packages for use on Guillimn
