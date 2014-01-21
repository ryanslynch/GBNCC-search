#!/usr/bin/env python
from __future__ import print_function
import os, sys
sys.path.insert(0, "pipeline/lib/python")
import config

print("Checking for the necessary python packages...", end="")
sys.stdout.flush()
packages = ["os", "sys", "shutil", "stat", "socket", "struct", "subprocess",
            "glob", "contextlib", "signal", "unittest", "re", "platform",
            "tarfile", "time", "datetime", "pytz", "MySQLdb", "xml",
            "numpy", "scipy", "pylab", "mpfit", "pyfits", "cPickle", "pickle",
            "psr_utils", "pyslalib", "presto", "prepfold", "pypsrcat",
            "sifting", "ratings", "diagnostics", "profile_tools", "infodata",
            "analyse_sp", "group_sp", "sp_utils",
            "config", "utils", "database", "handle_exit", "msub", "matplotlib",
            "__future__"]
for package in packages:
    try:
        __import__(package)
    except ImportError:
        print("\nERROR: Could not find %s"%package)
        print("Exiting")
        sys.exit(1)
print("sucess")

print("Testing database connection...", end="")
sys.stdout.flush()
import database
try:
    db = database.Database("observations")
except:
    print("\nWARNING: Could not connect to observations database")
    print("Continuing...")
else:
    db.close()
    print("success")

print("Attempting to create directory structure...", end="")
sys.stdout.flush()
for dir in [config.datadir,os.path.join(config.datadir,"holding"),
            config.jobsdir,config.logsdir,config.baseoutdir]:
    if os.path.exists(dir) and os.path.isdir(dir):
        pass
    else:
        try:
            os.makedirs(dir)
        except OSError:
            print("\nERROR: Could not create %s"%dir)
            print("Exiting")
            sys.exit(1)
print ("success")

print("To begin processing add\n"\
      "{0}/bin to your PATH and\n"\
      "{0}/lib/python to your PYTHONPATH\n"\
      "and run driver.sh".format(config.pipelinedir))
