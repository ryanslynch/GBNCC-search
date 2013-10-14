#!/usr/bin/python
import os, sys, pyfits, datetime
import database
from psr_utils import ra_to_rad, dec_to_rad, RADTODEG, DEGTORAD, SECPERDAY
from pyslalib.slalib import sla_eqgal

def parse_files(filenms):
    info = []
    for filenm in filenms:
        d = {}
        hdulist = pyfits.open(filenm, ignore_missing_end=True)
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
        

def upload(filenms):
     db = database.Database("observations")
     info = parse_files(filenms)
     for i in info:
         query = "INSERT INTO GBNCC (%s) VALUES (%s)" % \
             (",".join(i.keys()), ",".join(map(database._strfmt, i.values())))
         db.execute(query)

     db.commit()
     db.close()

     return True

if __name__ == "__main__":
    filenms = sys.argv[1:]
    status = upload(filenms)
    if status: print("All files uploaded successfully")
    
