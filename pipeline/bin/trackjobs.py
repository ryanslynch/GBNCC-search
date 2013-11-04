#!/usr/bin/env python
import os, time, msub, utils, config
import database as DB

print("Starting GBNCC job tracker...")
while True:
    db    = DB.Database("observations")
    query = "SELECT ID,ProcessingID,FileName FROM GBNCC WHERE ProcessingStatus='p' AND ProcessingSite='%s'"%config.machine
    db.execute(query)
    ret   = db.fetchall()
    
    if len(ret) != 0:
        for ID,jobid,filenm in ret:
            if msub.is_job_done(jobid):
                MJD,beamid = filenm.split("_")[1:3]
                outdir = os.path.join(config.baseoutdir, MJD, beamid)
                status = utils.results_status(outdir)
                print("Job %s completed with status %s"%(jobid,status))
                query = "UPDATE GBNCC SET ProcessingStatus='%s' "\
                        "WHERE ProcessingID=%i"%(status,jobid)
                db.execute(query)
                
            else: pass
    
    db.close()
    time.sleep(30*60)
