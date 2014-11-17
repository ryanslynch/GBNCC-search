#!/usr/bin/env python
from __future__ import print_function
import os, sys, subprocess, time, handle_exit
import database, config, utils

def download(outdir):
     db = database.Database("observations")
     #query  = "SELECT ID,FilePath,FileName FROM GBNCC WHERE "\
     #         "ProcessingStatus='u' OR (ProcessingStatus='f' AND "\
     #         "ProcessingAttempts < 10)"
     query  = "SELECT ID,FilePath,FileName FROM GBNCC WHERE "\
              "ProcessingStatus='u'"
     db.execute(query)
     ret     = db.cursor.fetchone()
     if ret is not None:
          ID      = ret[0]
          filenm  = os.path.join(*ret[1:])
          cmd     = "rsync astro.cv.nrao.edu:%s %s"%(filenm,outdir)
          
          query   = "UPDATE GBNCC SET ProcessingStatus='d',"\
                    "ProcessingSite='%s' WHERE ID=%i"%(config.machine,ID)
          db.execute(query)
          db.commit()
          retcode = subprocess.call(cmd, shell=True)
          if retcode == 0:
               query = "UPDATE GBNCC SET ProcessingStatus='q' WHERE ID=%i"%ID
               db.execute(query)
               db.commit()
               print("Successfully downloaded %s"%filenm)
          else:
               query = "UPDATE GBNCC SET ProcessingStatus='u', "\
                       "ProcessingSite=NULL WHERE ID=%i"%ID
               db.execute(query)
               db.commit()
               print("ERROR: Failed to download %s"%filenm)
     
     else:
          print("No files to download. Sleeping...")
          time.sleep(30*60)
     
     db.close()
    
     return None

def cleanup():
    query = "UPDATE GBNCC SET ProcessingStatus='u',ProcessingSite=NULL "\
            "WHERE ProcessingStatus='d' AND ProcessingSite='%s'"%config.machine
    db = database.Database("observations")
    db.execute(query)
    db.close()
    

def main():
    if len(sys.argv) == 2:
        outdir = sys.argv[1]
    else:
        outdir = config.datadir

    while True:
        if utils.get_size(config.datadir) < config.datadir_lim:
            download(outdir)
        else:
            print("Disk full. Sleeping...")
            time.sleep(config.sleeptime)

if __name__ == "__main__":
    with handle_exit.handle_exit(cleanup):
         print("Starting GBNCC downloader...")
         main()
