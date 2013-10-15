#!/usr/bin/env python
import os, glob, time, datetime, pytz, config, utils
import database as DB

while True:
    checkpoints = glob.glob(os.path.join(config.jobsdir, "*.checkpoint"))
    filenms = glob.glob(os.path.join(config.datadir, "guppi*GBNCC*fits"))
    nqueued = utils.getqueue(config.machine)

    while nqueued<config.queuelim and (len(filenms)>0 or len(checkpoints)>0):
        if len(checkpoints) > 0:
            checkpoint = checkpoints.pop()
            basenm,hashnm = checkpoint.split(".")[0:2]
            filenm = basenm + ".fits"
            with open(checkpoint, "r") as f: nodenm = f.readline().strip()
        
        elif len(filenms) > 0:
            filenm = filenms.pop()
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
                                              tmpdir=tmpdir, 
                                              outdir=config.baseoutdir,
                                              logsdir=config.logsdir,
                                              nodenm=nodenm, 
                                              zaplist=config.zaplist,
                                              pipelinedir=config.pipelinedir,
                                              walltimelim=config.walltimelim, 
                                              email=config.email))
        subfile.close()
        jobid,msg = utils.subjob(config.machine, subfile)
        if jobid is None: 
            print("ERROR: %s: %s"%(jobnm,msg))
        
        else:
            date = datetime.datetime.now()
            db = DB.Database(DB.databases["observations"])
            query = "UBDATE GBNCC SET ProcessingStatus='p',"\
                    "ProcessingID={jobid},ProcessingSite='{site}'"\
                    "ProcessingAttempts=ProcessingAttemps+1,"\
                    "ProcessingDate={date},PipelineVersion={version} "\
                    "WHERE FileName={filenm}".format(jobid=jobid,
                                                     site=config.institution,
                                                     date=date.isoformat(),
                                                     version=config.version,
                                                     filenm=filenm)
            nodenm = msub.get_all_jobs[jobid]["MasterHost"]
            checkpoint = os.path.join(config.jobsdir, jobnm+".checkpoint")
            with open(checkpoint, "w") as f:
                f.write(nodenm+"\n")
                f.write("0 0\n")
            db.execute(query)
            db.commit()
            db.close()
            time.sleep(5)
            nqueued = utils.getqueue(config.machine)
            
    else: time.sleep(config.sleeptime)
        
