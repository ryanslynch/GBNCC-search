import os, glob, msub, config

def getqueue(machine):
    if machine == "guillimin":
        alljobs = msub.get_all_jobs()
        if alljobs is not None:
            myjobs    = [job for job in alljobs.itervalues() \
                         if job.has_key("User") and job["User"] == config.user]
        else: myjobs = []
        nqueued = len([job for job in myjobs if job["State"] == "Idle" \
                       or job["State"] == "Blocked"])
    
    return nqueued

def subjob(machine, subfilenm, options=""):
    if machine == "guillimin":
        jobid,out,err = msub.submit_job(subfilenm, options)
    
    if   len(out) != 0: return jobid, out.strip()
    elif len(err) != 0: return jobid, err.strip()
    else: return jobid, None


def results_status(outdir, basenm):
    npfdplots   = len(glob.glob(os.path.join(outdir, "%s*.pfd.png"%basenm)))
    nspplots    = len(glob.glob(os.path.join(outdir, "%s*singlepulse.png"%basenm)))
    nratings    = len(glob.glob(os.path.join(outdir, "%s*.pfd.ratings"%basenm)))
    nrfifinds   = len(glob.glob(os.path.join(outdir, "%s*rfifind*"%basenm)))
    ntgzs       = len(glob.glob(os.path.join(outdir, "%s*.tgz"%basenm)))
    nreport     = len(glob.glob(os.path.join(outdir, "%s*.report"%basenm)))
    ndiagnostic = len(glob.glob(os.path.join(outdir, "%s*.diagnostics"%basenm)))
    naccels     = len(glob.glob(os.path.join(outdir, "%s*accelcands*"%basenm)))
    #ngroup      = len(glob.glob(os.path.join(outdir, "groups.txt"%basenm)))
    ngroup = 1
    #ngroupplots = len(glob.glob(os.path.join(outdir, "grouped*png")))
    ngroupplots = 1

    if (naccels != 2) or (ntgzs !=8) or (nrfifinds < 8) or (nreport != 1) or \
       (npfdplots == 0) or (nspplots < 6):
        return "f"
    elif (npfdplots != nratings) or (ndiagnostic != 1) or (ngroup != 1) or \
         (ngroupplots == 0):
        return "w"
    else:
        return "s"

def get_size(path):
    size = 0
    for dirpath,dirnms,filenms in os.walk(path):
        for f in filenms:
            size += os.path.getsize(os.path.join(dirpath,f))

    return size


