import os, subprocess

# Name of institution where pipeline is being run
institution = "McGill"
# Name of HPC machine where pipeline is being run
machine     = "guillimin"
# Timezone of processing site
timezone    = "Canada/Eastern"
# User name on 'machine'
user        = "rlynch"
# Email address where job notifications will be sent (if enabled)
email       = "rlynch+gbncc_jobs@physics.mcgill.ca"
# Walltime limit (hh:mm:ss)
walltimelim = "150:00:00"
# Maximum size of the 'pending' job queue
queuelim    = 30
# Time to wait between submitting a new job when there are no new files or the
# 'pending' queue is full
sleeptime   = 5*60
# Disk quota size of datadir (in bytes)
datadir_lim = 1024**4 - 12*1024**3 # 1 TB - 12 GB of overhead
# Top level analysis directory
topdir      = "/sb/project/bgf-180-ac/GBNCC"
# Base working directory for data reduction (should have at least 13 GB free)
baseworkdir = "/localscratch"
# Base temporary directory for data reduction (should have at least 2 GB free)
basetmpdir  = "/localscratch"
# Directory where pipeline scripts are stored
pipelinedir = os.path.join(topdir, "pipeline")
# Directory where raw data files are stored before being processed
datadir     = "/gs/scratch/rlynch"
# Directory where job submission files are stored
jobsdir     = os.path.join(topdir, "jobs")
# Directory wehre log files are stored
logsdir     = os.path.join(topdir, "logs")
# Directory where output files are permanently stored
baseoutdir  = os.path.join(topdir, "results")
# Location of FFT zaplist
zaplist     = os.path.join(pipelinedir, "lib", "GBNCC.zaplist")
# Pipeline version (as the git hash)
version     = subprocess.Popen("cd %s ; git rev-parse HEAD 2> /dev/null"%pipelinedir,shell=True,stdout=subprocess.PIPE).stdout.readline().strip()
# Databases dictionary
DATABASES = {
    "observations" : {
    "dbnm"   : "CONTACT RYAN FOR DB NAME",
    "hostnm" : "CONTACT RYAN FOR HOSTNM NAME",
    "usernm" : "CONTACT RYAN FOR USERNM NAME",
    "passwd" : "CONTACT RYAN FOR PASSWD NAME",
        },
    }


# Dictionary for holding job submission scripts
subscripts = {"guillimin": 
"""#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -N {jobnm}
#PBS -M {email}
#PBS -m ae
#PBS -q sw
#PBS -l nodes={nodenm}:ppn=1
#PBS -l walltime={walltimelim}

export CPATH={workdir}/python-2.6.8/include:{workdir}/python-2.6.8/wx-2.8.12/include:/software/compilers/intel/composerxe-2011.4.191/mkl/include:/software/compilers/intel/composerxe-2011.4.191/mkl/include
export LIBRARY_PATH={workdir}/python-2.6.8/include:/software/compilers/intel/composerxe-2011.4.191/compiler/lib/intel64:/software/compilers/intel/composerxe-2011.4.191/mkl/lib/intel64:/software/compilers/intel/composerxe-2011.4.191/compiler/lib/intel64:/software/compilers/intel/composerxe-2011.4.191/mkl/lib/intel64
export PYTHONPATH={workdir}/python-2.6.8/lib64/python2.6/site-packages:{workdir}/python-2.6.8/lib/python2.6/site-packages
export LD_LIBRARY_PATH={workdir}/python-2.6.8/GotoBLAS_LAPACK/shared:{workdir}/python-2.6.8/wx-2.8.12/lib:{workdir}/python-2.6.8/lib/lib64:{workdir}/python-2.6.8/lib/lib:/software/compilers/intel/composerxe-2011.4.191/compiler/lib/intel64:/software/compilers/intel/composerxe-2011.4.191/mkl/lib/intel64:/sb/software/libraries/MKL/10.3/lib/intel64:/software/libraries/PGPLOT/5.2:/sb/software/libraries/CFITSIO/3.28/lib:/usr/local/lib:/usr/lib/:/sb/project/bgf-180-aa/lib:/sb/project/bgf-180-aa/lib64:/sb/project/bgf-180-aa/src/presto_memaccel/presto/lib
export PATH={workdir}/python-2.6.8/bin:{workdir}/python-2.6.8/wx-2.8.12/bin:/opt/moab/bin:/software/tools/datamover:/software/tools/swig-2.0.4/bin:/sb/software/libraries/MKL/10.3/bin:/software/tools/clig/bin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/usr/lpp/mmfs/bin:/sb/software/tools/scripts:/usr/lpp/mmfs/bin:/home/rlynch/bin:/sb/project/bgf-180-aa/bin:/sb/project/bgf-180-aa/src/presto_memaccel/presto/bin

if [ {nodenm} == 1 ]
  then
    echo -e \"$HOSTNAME
{jobid}
0 0\" > {jobsdir}/{jobnm}.checkpoint
    mkdir -p {workdir}
    mv {filenm} {workdir}
    cp {zaplist} {workdir}
    tar -C {workdir} -xzf {pipelinedir}/lib/python-2.6.8-packages-GBNCC.tar.gz
  else
    set -- $({jobsdir}/{jobnm}.checkpoint)
    echo -e \"$HOSTNAME
{jobid}
$3 $4\" > {jobsdir}/{jobnm}.checkpoint
    mv {baseworkdir}/$2 {baseworkdir}/{jobid}
fi
cd {workdir}
python-2.6.8/bin/search.py -w {workdir} -i {hashnm} {basenm}.fits 
#rm -rf {workdir}
""",

"condor":
"""Executable={pipelinedir}/bin/search.py
Universe=vanilla
# New configuration
+Online_GBNCC = True
Requirements = TARGET.Online_GBNCC =?= True && Machine == "nemo-slave1081.nemo.phys.uwm.edu"
Arguments= -s {tmpdir} -w {workdir} -o {outdir} -i {basenm}.fits -z {zaplist}
getenv=True
Output={logsdir}/GBNCC.$(cluster).$(process).out
Error={logsdir}/GBNCC.$(cluster).$(process).err
Log={logsdir}/GBNCC.log
queue
"""
}

subscript = subscripts[machine]
