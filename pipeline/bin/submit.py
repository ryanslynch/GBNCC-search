#!/usr/bin/env python
import os, shutil, glob, time, datetime, pytz, config, utils, database, msub

#checkpoints = glob.glob(os.path.join(config.jobsdir, "*.checkpoint"))
checkpoints = []

print("Starting GBNCC job submitter...")
while True:
    filenms = glob.glob(os.path.join(config.datadir, "guppi*GBNCC*fits"))
    nqueued = utils.getqueue(config.machine)

    while nqueued<config.queuelim and (len(filenms)>0 or len(checkpoints)>0):
        if len(checkpoints) > 0:
            checkpoint = checkpoints.pop()
            basenm,hashnm = checkpoint.split(".")[0:2]
            basenm = os.path.basename(basenm)
            filenm = basenm + ".fits"
            with open(checkpoint, "r") as f: nodenm = f.readline().strip()
        
        elif len(filenms) > 0:
            filenm = filenms.pop()
            shutil.move(filenm, os.path.join(config.datadir, "holding"))
            filenm = os.path.join(config.datadir, "holding",
                                  os.path.basename(filenm))
            basenm = os.path.basename(filenm).rstrip(".fits")
            hashnm = os.urandom(8).encode("hex")
            nodenm = "1"
        
        jobnm   = basenm + "." +  hashnm
        workdir = os.path.join(config.baseworkdir, basenm, hashnm)
        tmpdir  = os.path.join(config.basetmpdir, basenm, hashnm, "tmp")

        subfilenm = os.path.join(config.jobsdir, jobnm+".sh")
        subfile   = open(subfilenm, "w")
        subfile.write(config.subscript.format(filenm=filenm, basenm=basenm, 
                                              jobnm=jobnm, workdir=workdir,
                                              hashnm=hashnm,
                                              tmpdir=tmpdir, 
                                              outdir=config.baseoutdir,
                                              logsdir=config.logsdir,
                                              nodenm=nodenm, 
                                              zaplist=config.zaplist,
                                              pipelinedir=config.pipelinedir,
                                              walltimelim=config.walltimelim, 
                                              email=config.email))
        subfile.close()
        jobid,msg = utils.subjob(config.machine,subfilenm,options="-o {0} -e {0}".format(config.logsdir))
        if jobid is None: 
            print("ERROR: %s: %s"%(jobnm,msg))
        
        else:
            print("Submitted %s"%jobnm)
            while msub.get_all_jobs()[jobid]["State"] == "Idle":
                time.sleep(5)
            
            date = datetime.datetime.now()
            nodenm = msub.get_all_jobs()[jobid]["MasterHost"]
            checkpoint = os.path.join(config.jobsdir, jobnm+".checkpoint")
            with open(checkpoint, "w") as f:
                f.write(nodenm+"\n")
                f.write("0 0\n")
            db = database.Database("observations")
            query = "UPDATE GBNCC SET ProcessingStatus='p',"\
                    "ProcessingID='{jobid}',ProcessingSite='{site}',"\
                    "ProcessingAttempts=ProcessingAttempts+1,"\
                    "ProcessingDate='{date}',PipelineVersion='{version}' "\
                    "WHERE FileName='{filenm}'".format(jobid=jobid,
                                                     site=config.machine,
                                                     date=date.isoformat(),
                                                     version=config.version,
                                                     filenm=os.path.basename(filenm))
            db.execute(query)
            db.commit()
            db.close()
            time.sleep(5)
            nqueued = utils.getqueue(config.machine)
            
    else:
        print("Nothing to submit.  Sleeping...")
        time.sleep(config.sleeptime)
        
