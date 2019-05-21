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
walltimelim = "120:00:00"                                        
# Maximum size of the 'pending' job queue                        
queuelim    = 30                                                 
# Time to wait between submitting a new job when there are no new files or the
# 'pending' queue is full                                                     
sleeptime   = 5*60                                                            
# Disk quota size of datadir (in bytes)                                       
datadir_lim = 1024**4 - 12*1024**3 # 1 TB - 12 GB of overhead                 
# Top level analysis directory                                                
topdir      = "/gs/project/bgf-180-ad/GBNCC"                                  
# Base working directory for data reduction (should have at least 13 GB free) 
baseworkdir = "/localscratch"                                                 
# Base temporary directory for data reduction (should have at least 2 GB free)
basetmpdir  = "/localscratch"                                                 
# Directory where pipeline scripts are stored                                 
pipelinedir = os.path.join(topdir, "pipeline")                                
# Directory where third party software is found
softwaredir = "/gs/project/bgf-180-ad/PALFA4/software"
# Directory where raw data files are stored before being processed            
datadir     = "/gs/scratch/rlynch/GBNCC"                                      
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
    "hostnm" : "CONTACT RYAN FOR HOST NAME",                                    
    "usernm" : "CONTACT RYAN FOR USERNAME",                                     
    "passwd" : "CONTACT RYAN FOR PASSWORD",                                     
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
#PBS -A bgf-180-ad

if [ {nodenm} == 1 ]
  then
    echo -e \"$HOSTNAME
{jobid}
0 0\" > {jobsdir}/{jobnm}.checkpoint
    mkdir -p {workdir}
    mv {filenm} {workdir}
    cp {zaplist} {workdir}
  else
    set -- $({jobsdir}/{jobnm}.checkpoint)
    echo -e \"$HOSTNAME
{jobid}
$3 $4\" > {jobsdir}/{jobnm}.checkpoint
    mv {baseworkdir}/$2 {baseworkdir}/{jobid}
fi
cd {workdir}
search.py -w {workdir} -i {hashnm} {basenm}.fits
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
""",

"nimrod":
"""#!/bin/bash
#PBS -m a
#PBS -u {user}
#PBS -V
#PBS -M {email}
#PBS -N {jobnm}
#PBS -l nodes=1:compute:ppn=1
#PBS -l walltime={walltimelim}

mkdir -p {workdir}
mv {filenm} {workdir}
cp {zaplist} {workdir}
cd {workdir}
search.py -w {workdir} -i {hashnm} {basenm}.fits
#rm -rf {workdir}
"""
}

subscript = subscripts[machine]
