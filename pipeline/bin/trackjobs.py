#!/usr/bin/env python
import os, time, utils, config, PBSQuery
import database as DB

print("Starting GBNCC job tracker...")
queue = PBSQuery.PBSQuery()
while True:
    db    = DB.Database("observations")
    query = "SELECT ID,ProcessingID,FileName FROM GBNCC WHERE (ProcessingStatus='p' OR ProcessingStatus='i') AND ProcessingSite='%s'"%config.machine
    db.execute(query)
    ret   = db.fetchall()
    alljobs = queue.getjobs()
    
    if len(ret) != 0 and alljobs is not None:
        for ID,jobid,filenm in ret:
            if not alljobs.has_key(str(jobid)):
                MJD,beamid = filenm.split("_")[1:3]
                basenm = filenm.strip(".fits")
                outdir = os.path.join(config.baseoutdir, MJD, beamid)
                status = utils.results_status(outdir, basenm)
                print("Job %s completed with status %s"%(jobid,status))
                query = "UPDATE GBNCC SET ProcessingStatus='%s' "\
                        "WHERE ID=%i"%(status,ID)
                db.execute(query)
                
            else: pass
    
    db.close()
    time.sleep(30*60)
