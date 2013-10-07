# Script for running multiple single-node batch CLANS2 jobs
# John Zollweg, Cornell Center for Advanced Computing, April 20, 2009
# last modified April 21, 2009 - JAZ
# usage: clans2batch1.py num
# i.e.: clans2batch1.py 15

__author__ = "John Zollweg (zollweg@cac.cornell.edu)"
__version__ = "$Revision: 14 $"
__date__ = "$Date: 2009-06-11 15:22:36 -0400 (Thu, 11 Jun 2009) $"

import os,sys,shutil,string,time
from subprocess import *
import nsub

if len(sys.argv) < 2:
    print """clans2batch1.py runs a number of single-node batch jobs
    Usage: clans2batch1.py num
           where:  the first argument is the number of jobs to run
                       """
    sys.exit()

jno = []  # Id of job given by scheduler
try:
    njobs = int(sys.argv[1])
except ValueError, err:
    print "Could not read second argument as number of runs.", sys.argv[1]
    sys.exit()

batchfiles = []
for j in range(1,njobs+1):
    job = os.path.join("..","job%i.bat" % j)
    if os.path.isfile(job):
        batchfiles.append(job)
    else:
        print "Could not find batch file %s. Exiting." % job
        sys.exit()

anyFailed=False
for batchfile in batchfiles:
    (job_number,out,err) = nsub.submit_job(batchfile)
    if job_number:
        jno.append(job_number)
    else:
        print 'failed to submit job for batchfile',batchfile
        anyFailed=True
        break

if anyFailed:
    for jobid in jno:
        nsub.cancel_job(jobid)
else:
    while not nsub.are_jobs_done(jno):
        time.sleep(10)

