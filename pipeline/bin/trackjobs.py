#!/usr/bin/env python
import time, msub, utils, config
import database as DB

while True:
    db    = DB.Database(DB.databases["observations"])
    query = "SELECT ProcessingID,FileName FROM GBNCC WHERE ProcessingStatus='p'"
    db.execute(query)
    ret   = db.cursor().fetchall()
    
    for ID,jobid,filenm in ret:
        if msub.is_job_done(jobid):
            MJD,beamid = filenm.split("_")[1:3]
            outdir = os.path.join(config.baseoutdir, MJD, beamid)
            status = results_status(outdir)
                query = "UPDATE GBNCC SET ProcessingStatus='%s' "\
                        "WHERE ProcessingID=%i"%(status,jobid)
                db.execute(query)
                
        else: pass
    
    db.close()
    time.sleep(30*60)
