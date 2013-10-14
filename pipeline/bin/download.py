#!/usr/bin/python
import os, sys, subprocess, atexit
import database, config

def download(outdir):
     db = database.Database("observations")
     query  = "SELECT ID,FilePath,FileName FROM GBNCC WHERE ProcessingStatus='u' OR (ProcessingStatus='f' AND ProcessingAttempts < 3)"
     db.execute(query)
     ret     = db.cursor.fetchone()
     ID      = ret[0]
     filenm  = os.path.join(*ret[1:])
     cmd     = "rsync -aux pulsar.physics.mcgill.ca:%s %s"%(filenm,outdir)

     query   = "UPDATE GBNCC SET ProcessingStatus='d',ProcessingSite='%s' "\
               "WHERE ID=%i"%(config.machine,ID)
     db.execute(query)
     db.commit()
     retcode = subprocess.call(cmd, shell=True)
     if retcode == 0:
         query = "UPDATE GBNCC SET ProcessingStatus='q' WHERE ID=%i"%ID
         db.execute(query)
         db.commit()
         print("Successfully downloaded %s"%filenm)
     else:
         query = "UPDATE GBNCC SET ProcessingStatus='u' WHERE ID=%i"%ID
         db.execute(query)
         db.commit()
         print("ERROR: Failed to download %s"%filenm)
     
     db.close()
    
     return None

def cleanup():
    query = "UPDATE GBNCC SET ProcessingStatus='u',ProcessingSite=NULL "\
            "WHERE ProcessingSite='%s'"%config.machine
    db = database.Database("observations")
    db.execute(query)
    db.close()
    

def main():
    if len(sys.argv) == 2:
        outdir = sys.argv[1]
    else:
        outdir = config.datadir

    while True:
        download(outdir)
    

if __name__ == "__main__":
    with atexit.handle_exit(cleanup):
        main()
