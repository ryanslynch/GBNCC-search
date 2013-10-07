import os, glob, msub

def getqueue(machine):
    if machine == "guillimin":
        jobs    = [job for job in nsub.get_all_jobs().itervalues() \
                       if job["User"] == config.user]
        nqueued = len([job for job in myjobs if job["State"] == "Idle" \
                       or job["State"] == "Blocked"])
    
    return None

def subjob(machine,subfilenm):
    if machine == "guillimin":
        jobid,out,err = nsub.submit_job(subfilenm)
    
    if   len(out) != 0: return jobid, out.strip()
    elif len(err) != 0: return jobid, err.strip()
    else: return jobid, None


def results_status(outdir):
    npfdplots   = len(glob.glob(os.path.join(outdir, "*pfd.png")))
    nspplots    = len(glob.glob(os.path.join(outdir, "*.singlepulse.png")))
    nratings    = len(glob.glob(os.path.join(outdir, "*.pfd.ratings")))
    nrfifinds   = len(glob.glob(os.path.join(outdir, "*rfifind*")))
    ntgzs       = len(glob.glob(os.path.join(outdir, "*.tgz")))
    nreport     = len(glob.glob(os.path.join(outdir, "*.report")))
    ndiagnostic = len(glob.glob(os.path.join(outdir, "*.diagnostics")))
    naccels     = len(glob.glob(os.path.join(outdir, "*accelcands*")))

    if (naccels != 2) or (ntgz !=8) or (nrfifinds < 8) or (nreport != 1) or \
       (npfdplots == 0) or (nspplots != 6):
        return "f"
    elif (npfdplots != nratings) or (ndiagnostic != 1):
        return "w"
    else:
        return "s"
        
