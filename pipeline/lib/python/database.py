import os, sys, fnmatch, tempfile, optparse, MySQLdb, pyfits, datetime
from psr_utils import ra_to_rad, dec_to_rad, RADTODEG, DEGTORAD, SECPERDAY
from pyslalib.slalib import sla_eqgal
import config, subprocess

hostname = "pulsar.physics.mcgill.ca"
username = "gbtpsr"
password = "NRAO100-m"
dbname   = "GBTDataDB"


databases = {
    "observations" : {
        "dbnm"   : "GBTDataDB",
        "hostnm" : "pulsar.physics.mcgill.ca",
        "usernm" : "gbtpsr",
        "passwd" : "NRAO100-m",
        },
    }

class Database(object):
    def __init__(self, database):
        self.db = MySQLdb.connect(host=database["hostnm"],
                                  user=database["usernm"],
                                  db=database["dbname"], 
                                  passwd=databse["passwd"])
        self.cursor = self.db.cursor()

    def execute(self, query):
        try:
            self.cursor.execute(query)
        except Exception,e:
            print("ERROR: %i: %s"%(e.args[0], e.args[1]))
            self.db.rollback()
            self.db.close()
            exit(1)
    
    def __del__(self):
        self.db.close()

def _strfmt(x):
    if type(x) == str and x != "NULL":
        return "\"%s\""%x
    else:
        return str(x)




def parse_files(filenms):
    info = []
    for filenm in filenms:
        d = {}
        hdulist = pyfits.open(filenm)
        hdr0 = hdulist[0].header
        hdr1 = hdulist[1].header
        hdulist.close()
        d["BeamID"] = int(hdr0["SRC_NAME"].strip("GBNCC"))
        d["RightAscension"] = ra_to_rad(hdr0["RA"])*RADTODEG
        d["Declination"] = dec_to_rad(hdr0["DEC"])*RADTODEG
        l,b = sla_eqgal(d["RightAscension"]*DEGTORAD,d["Declination"]*DEGTORAD)
        d["GalacticLon"] = l*RADTODEG
        d["GalacticLat"] = b*RADTODEG
        d["DateObserved"] = hdr0["STT_IMJD"] + hdr0["STT_SMJD"]/SECPERDAY
        d["CenterFreq"] = hdr0["OBSFREQ"]
        d["BW"] = abs(hdr0["OBSBW"])
        d["SamplingTime"] = hdr1["TBIN"]*1000000
        d["IntTime"] = hdr1["NAXIS2"]*hdr1["NSBLK"]*hdr1["TBIN"]
        d["FileName"] = os.path.split(filenm)[1]
        d["FilePath"] = os.path.split(os.path.abspath(filenm))[0]
        d["ArchiveLocation"] = "NULL"
        d["FileSize"] = os.path.getsize(filenm)
        d["ProcessingStatus"] = "u"
        d["ProcessingAttempts"] = 0
        d["ProcessingSite"] = "NULL"
        d["ProcessingID"] = "NULL"
        d["ProcessingDate"] = "NULL"
        d["PipelineVersion"] = "NULL"
        d["UploadDate"] = datetime.datetime.fromtimestamp(os.path.getmtime(filenm)).strftime("%Y-%m-%dT%H:%M:%S")
        info.append(d)

    return info
        
        

def add_entries(filenms, dbname=dbname):
     db = Database(databases["observations"])
     info = parse_files(filenms)
     for i in info:
         query = "INSERT INTO GBNCC (%s) VALUES (%s)" % \
             (",".join(i.keys()), ",".join(map(_strfmt, i.values())))

         db.execute(query)

     db.commit()
     db.close()

     return True


 def download(dbname=dbname):
     db = Database(databases["observations"])
     query  = "SELECT ID,FilePath,FileName FROM GBNCC WHERE ProcessingStatus='u' OR (ProcessingStatus='f' AND ProcessingAttemps < 3)"
     db.execute(query)
     ret     = db.cursor.fetchone()
     ID      = ret[0]
     filenm  = os.path.join(*ret[1:])
     cmd     = "rsync -aux astro.cv.nrao.edu:%s %s"%(filenm,config.datadir)

     query   = "UPDATE GBNCC SET ProcessingStatus='d' WHERE ID=%i"%ID
     db.execute(query)
     db.commit()
     retcode = subprocess.call(cmd, shell=True)
     if retcode == 0:
         query = "UPDATE GBNCC SET ProcessingStatus='q' WHERE ID=%i"%ID
         db.execute(query)
         db.commit()
     else:
         query = "UPDATE GBNCC SET ProcessingStatus='u' WHERE ID=%i"%ID
         db.execute(query)
         db.commit()
     
     db.close()
    
    return None
            

if __name__ == "__main__":
    filenms = sys.argv[1:]
    add_entries(filenms)
