#!/usr/bin/env python
import os, sys, numpy, datetime
import database
from astropy.io import fits
from astropy import coordinates, units

# The following constants and functions were taken from psr_constants
# and psr_utils

RADTODEG    = float('57.295779513082320876798154814105170332405472466564')
DEGTORAD    = float('1.7453292519943295769236907684886127134428718885417e-2')
SECTORAD    = float('7.2722052166430399038487115353692196393452995355905e-5')
ARCSECTORAD = float('4.8481368110953599358991410235794797595635330237270e-6')
SECPERDAY   = float('86400.0')

def hms_to_rad(hour, min, sec):
    """
    hms_to_rad(hour, min, sec):
       Convert hours, minutes, and seconds of arc to radians
    """
    if (hour < 0.0): sign = -1
    else: sign = 1
    return sign * SECTORAD * \
           (60.0 * (60.0 * numpy.fabs(hour) +
                    numpy.fabs(min)) + numpy.fabs(sec))

def dms_to_rad(deg, min, sec):
    """
    dms_to_rad(deg, min, sec):
       Convert degrees, minutes, and seconds of arc to radians.
    """
    if (deg < 0.0):
        sign = -1
    elif (deg==0.0 and (min < 0.0 or sec < 0.0)):
        sign = -1
        
    else:
        sign = 1
    return sign * ARCSECTORAD * \
           (60.0 * (60.0 * numpy.fabs(deg) +
                    numpy.fabs(min)) + numpy.fabs(sec))

def ra_to_rad(ra_string):
    """
    ra_to_rad(ar_string):
       Given a string containing RA information as
       'hh:mm:ss.ssss', return the equivalent decimal
       radians.
    """
    h, m, s = ra_string.split(":")
    return hms_to_rad(int(h), int(m), float(s))

def dec_to_rad(dec_string):
    """
    dec_to_rad(dec_string):
       Given a string containing DEC information as
       'dd:mm:ss.ssss', return the equivalent decimal
       radians.
    """
    d, m, s = dec_string.split(":")
    if "-" in d and int(d)==0:
        m, s = '-'+m, '-'+s
    return dms_to_rad(int(d), int(m), float(s))


def parse_files(filenms):
    info = []
    for filenm in filenms:
        d = {}
        try:
            hdulist = fits.open(filenm, ignore_missing_end=True)
            hdr0 = hdulist[0].header
            hdr1 = hdulist[1].header
            hdulist.close()
            d["BeamID"] = int(hdr0["SRC_NAME"].strip("GBNCC"))
            d["RightAscension"] = ra_to_rad(hdr0["RA"])*RADTODEG
            d["Declination"] = dec_to_rad(hdr0["DEC"])*RADTODEG
            c = coordinates.ICRSCoordinates(ra=d["RightAscension"],
                                            dec=d["Declination"],
                                            unit=(units.degree,units.degree))
            d["GalacticLon"] = c.galactic.l.degrees
            d["GalacticLat"] = c.galactic.b.degrees
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

        except:
            print("Failed to read %s"%filenm)
    
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
    if sys.argv[1] == "--from-file":
        filenms = numpy.loadtxt(sys.argv[2], dtype=str)
    else:
        filenms = sys.argv[1:]
    db = database.Database("observations")
    query = "SELECT FileName FROM GBNCC"
    db.execute(query)
    ret = numpy.array(db.fetchall())
    filenms = [filenm for filenm in filenms if os.path.basename(filenm) not in ret]
    status = upload(filenms)
    if status: print("All files uploaded successfully")
    
