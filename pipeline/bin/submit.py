#!/usr/bin/env python
import os, subprocess, shutil, glob, time, datetime, pytz, config, utils, database, PBSQuery
prestodir = os.environ["PRESTO"]

#checkpoints = glob.glob(os.path.join(config.jobsdir, "*.checkpoint"))
checkpoints = []

queue = PBSQuery.PBSQuery()

print("Starting GBNCC job submitter...")
while True:
    db = database.Database("observations")
    query = "SELECT ProcessingID,FileName FROM GBNCC WHERE "\
            "ProcessingStatus='i'"
    db.execute(query)
    ret = db.fetchall()
    if len(ret) != 0:
        for jobid,filenm in ret:
            alljobs = queue.getjobs()
            if alljobs is not None:
                if alljobs.has_key(str(jobid)):
                    if alljobs[str(jobid)]["job_state"][0] == "R":
                        nodenm = alljobs[str(jobid)]["exec_host"][0]
                        jobnm  = alljobs[str(jobid)]["Job_Name"][0]
                        #checkpoint = os.path.join(config.jobsdir, jobnm+".checkpoint")
                        #with open(checkpoint, "w") as f:
                        #    f.write(nodenm+"\n")
                        #    f.write("0 0\n")
                        query = "UPDATE GBNCC SET ProcessingStatus='p' "\
                                "WHERE FileName='{filenm}'".format(filenm=filenm)
                        db.execute(query)
                    else:
                        pass
    db.close()
    
    filenms = glob.glob(os.path.join(config.datadir, "guppi*GBNCC*fits"))
    nqueued = utils.getqueue(config.machine,queue)

    while nqueued<config.queuelim and (len(filenms)>0 or len(checkpoints)>0):
        if len(checkpoints) > 0:
            checkpoint = checkpoints.pop()
            basenm,hashnm = checkpoint.split(".")[0:2]
            basenm = os.path.basename(basenm)
            filenm = basenm + ".fits"
            #with open(checkpoint, "r") as f:
            #    nodenm = f.readline().strip()
            #    jobid  = f.readline().strip()
        
        elif len(filenms) > 0:
            filenm = filenms.pop()
            shutil.move(filenm, os.path.join(config.datadir, "holding"))
            filenm = os.path.join(config.datadir, "holding",
                                  os.path.basename(filenm))
            basenm = os.path.basename(filenm).rstrip(".fits")
            hashnm = os.urandom(8).encode("hex")
            nodenm = "1"
            jobid  = "$PBS_JOBID"
        
        jobnm   = basenm + "." +  hashnm
        workdir = os.path.join(config.baseworkdir, jobid, basenm, hashnm)
        tmpdir  = os.path.join(config.basetmpdir, jobid, basenm, hashnm, "tmp")

        subfilenm = os.path.join(config.jobsdir, jobnm+".sh")
        subfile   = open(subfilenm, "w")
        subfile.write(config.subscript.format(filenm=filenm, basenm=basenm, 
                                              jobnm=jobnm, workdir=workdir,
                                              baseworkdir=config.baseworkdir,
                                              hashnm=hashnm, jobid=jobid,
                                              jobsdir=config.jobsdir,
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
            prestoversion = subprocess.Popen("cd %s ; git rev-parse HEAD 2> /dev/null"%prestodir,shell=True,stdout=subprocess.PIPE).stdout.readline().strip()
            jobid = jobid.strip()
            print("Submitted %s with ID %s"%(jobnm,jobid))
            time.sleep(30)
            alljobs = queue.getjobs()
            if alljobs is not None:
                if alljobs[jobid]["job_state"][0] == "Q":
                    status = "i"
                else:
                    status = "p"
                    nodenm = alljobs[jobid]["exec_host"][0]
                    checkpoint = os.path.join(config.jobsdir, jobnm+".checkpoint")
                    with open(checkpoint, "w") as f:
                        f.write(nodenm+"\n")
                        f.write("0 0\n")

                date = datetime.datetime.now()
                query = "UPDATE GBNCC SET ProcessingStatus='{status}',"\
                        "ProcessingID='{jobid}',ProcessingSite='{site}',"\
                        "ProcessingAttempts=ProcessingAttempts+1,"\
                        "ProcessingDate='{date}',"\
                        "PipelineVersion='{version}', "\
                        "PRESTOVersion='{prestoversion}' "\
                        "WHERE FileName='{filenm}'".format(status=status,
                                                           jobid=jobid,
                                                           site=config.machine,
                                                           date=date.isoformat(),
                                                           version=config.version,
                                                           prestoversion=prestoversion,
                                                           filenm=os.path.basename(filenm))
                
                db = database.Database("observations")
                db.execute(query)
                db.commit()
                db.close()
        nqueued = utils.getqueue(config.machine,queue)
            
    else:
        print("Nothing to submit.  Sleeping...")
        time.sleep(config.sleeptime)
        
