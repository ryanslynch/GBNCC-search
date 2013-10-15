#!/usr/bin/env python

# Import necessary packages
import sys, glob
import numpy as N

LO_Z              = 0
HI_Z              = 50
FOLDING_THRESHOLD = 2.0

# Function definitions
def get_rfi_diagnostics(basename):
    """
    Get diagnostic information from the results of rfifind.

    Parameters
    ----------
    basename : str
        The base file name of this beam
    
    Returns
    -------
    mask_percentage : float
        The percentage of the data that has been masked
    """

    f     = open(basename + "_rfifind.out", "r")
    lines = f.readlines()
    f.close()

    for line in lines:
        if line.startswith("  Number of  bad   intervals:"):
            mask_prct = float(line.strip().split()[-1].strip("(").strip("%)"))

    return mask_prct

def get_fold_diagnostics(basename, lo_z, hi_z, folding_threshold):
    """
    Get diagnostic information from the results of rfifind.

    Parameters
    ----------
    basename : str
        The base file name of this beam
    lo_z : int
        The number of Fourier bins searched over in the low acceleration
        searches
    hi_z : int
        The number of Fourier bins searched over in the high acceleration
        searches
    folding_threshold : float
        The minimum presto sigma for which a candidate is passed to folding

    Returns
    -------
    num_cands_lo : int
        The total number of candidates from the low z accelsearch
    num_cands_hi : int
        The total number of candidates from the high z accelsearch
    num_folded_lo : int
        The total number of folded candidates from the low z accelsearch
    num_folded_hi : int
        The total number of folded candidates from the high z accelsearch
    num_above_threshold : int
        The number of candidates above the minimum presto sigma threshold
    min_sigma_folded : float
        The lowest presto sigma of the folded candidates
    """

    lo_accel_fn = basename + ".accelcands_Z%i"%lo_z
    hi_accel_fn = basename + ".accelcands_Z%i"%hi_z

    lo_accel_cands  = []
    lo_accel_sigmas = []
    f               = open(lo_accel_fn, "r")
    lines           = f.readlines()
    f.close()
    for line in lines:
        if line.startswith("guppi"):
            sline      = line.strip().split()
            candnum    = int(sline[0].split(":")[1])
            DM         = sline[1]
            sigma      = float(sline[3])
            accel_cand = "%s_DM%s_Z%i_ACCEL_Cand_%i" % \
                         (basename, DM, lo_z, candnum)
            lo_accel_cands.append(accel_cand)
            lo_accel_sigmas.append(sigma)
        else: pass

    hi_accel_cands  = []
    hi_accel_sigmas = []
    f               = open(hi_accel_fn, "r")
    lines           = f.readlines()
    f.close()
    for line in lines:
        if line.startswith("guppi"):
            sline      = line.strip().split()
            candnum    = int(sline[0].split(":")[1])
            DM         = sline[1]
            sigma      = float(sline[3])
            accel_cand = "%s_DM%s_Z%i_ACCEL_Cand_%i" % \
                         (basename, DM, hi_z, candnum)
            hi_accel_cands.append(accel_cand)
            hi_accel_sigmas.append(sigma)
        else: pass

    all_cands  = N.concatenate((lo_accel_cands, hi_accel_cands))
    all_sigmas = N.concatenate((lo_accel_sigmas, hi_accel_sigmas))

    num_cands_lo        = len(lo_accel_cands)
    num_cands_hi        = len(hi_accel_cands)
    num_above_threshold = sum(all_sigmas >= folding_threshold)

    lo_accel_folds = glob.glob("%s*Z%i*png"%(basename, lo_z))
    hi_accel_folds = glob.glob("%s*Z%i*png"%(basename, hi_z))
    num_folded_lo  = len(lo_accel_folds)
    num_folded_hi  = len(hi_accel_folds)

    all_folds     = N.concatenate((lo_accel_folds, hi_accel_folds))
    folded_sigmas = []
    for cand,sigma in zip(all_cands, all_sigmas):
        fold_name = cand + ".pfd.png"
        if fold_name in all_folds:
            folded_sigmas.append(sigma)
        else: pass
    min_sigma_folded = min(folded_sigmas)

    return (num_cands_lo, num_cands_hi, num_folded_lo, num_folded_hi,
            num_above_threshold, min_sigma_folded)


def write_diagnostics(basename):
    """
    Write diagnostic values to an output file

    Parameters
    ----------
    basename : str
        The base file name of this beam

    Returns
    -------
    None
    """
    mask_percentage = get_rfi_diagnostics(basename)
    num_cands_lo,num_cands_hi,num_folded_lo,num_folded_hi,\
        num_above_threshold,min_sigma_folded = \
        get_fold_diagnostics(basename, LO_Z, HI_Z, FOLDING_THRESHOLD)

    outfile = open(basename + ".diagnostics", "w")
    line    = "RFI mask percentage                 = %.3f\n"%mask_percentage
    outfile.write(line)
    line    = "Total low accel cands               = %i\n"%num_cands_lo
    outfile.write(line)
    line    = "Total high accel cands              = %i\n"%num_cands_hi
    outfile.write(line)
    line    = "Low accel cands folded              = %i\n"%num_folded_lo
    outfile.write(line)
    line    = "High accel cands folded             = %i\n"%num_folded_hi
    outfile.write(line)
    line    = "Total cands above folding threshold = %i\n"%num_above_threshold
    outfile.write(line)
    line    = "Min presto sigma folded             = %.2f\n"%min_sigma_folded
    outfile.write(line)
    outfile.close()

    return None

if __name__ == "__main__":
    beams = sys.argv[1:]
    for beam in beams:
        print "Working on %s"%beam
        report = glob.glob("%s/*.report"%beam)
        if len(report) == 1:
            basename = report[0].strip(".report")
            try:
                write_diagnostics(basename)
            except:
                print "WARNING: Diagnostics failed for %s"%beam
        else: pass

    
