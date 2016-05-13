#!/usr/bin/env python 

# Import the necessary packages
import sys, os
import numpy as N
import scipy as S
import pylab as PL

import prepfold as PF
import psr_utils as PU
import pypsrcat as PSRCAT
import profile_tools as PT
import cPickle
import ubc_AI.data

try:
    from pyslalib.slalib import sla_dsep
except ImportError:
    from slalib import sla_dsep

# Set some global variables

# The threshold probability (as determined by an F-test) for accepting an
# additional Gaussian component in a profile fit
F_STAT_THRESHOLD         = 0.05 
# The maximum number of Gaussian profile components to allow for
# multi-Gaussian fitting
MAX_GAUSSIANS            = 5
# The maximum phase drift to count a subint/subband as "good" in phase wiggle
# ratings
PHASE_DRIFT_TOLERANCE    = 0.02
# The peak SNR threshold for a "good" signal in a subint/subband for
# intermittence/broadbandedness ratings
PEAK_SNR_THRESHOLD       = 3.0
# The integrated SNR threshold for a "good" signal in a subint/subband for
# intermittence/broadbandedness ratings
INTEGRATED_SNR_THRESHOLD = 5.0
# The FWHM of the telescope beam (radians)
BEAM_FWHM                = 30.0/60.0*N.pi/180.0
# The UBC AI classsifier data
ubc_AI_classifier_name = "/gs/project/bgf-180-ad/PALFA3/software/lib/python2.7/site-packages/ubc_AI/trained_AI/clfl2_BD.pkl"
ubc_AI_classifier = cPickle.load(open(ubc_AI_classifier_name, "rb"))

# Utility function definitions
def calc_persistence_metrics(subints, on_pulse_bins, off_pulse_bins):
    """
    Calculate the fraction of sub-intervals above a SNR threshold and the
    larger of the peak or integrated SNR standard deviations.

    Parameters
    ----------
    subints : array_like
        The folded sub-intervals
    on_pulse_bins : array_like
        A list of on-pulse profile bins
    off_pulse_bins : array_like
        A list of off-pulse profile bins

    Returns
    -------
    frac_good_subints : float
        The fraction of \"good\" sub-intervals
    subints_snr_std : float
        The standard deviation of sub-intervals SNRs
    """
    # Initiate some variables
    good_subints            = 0
    zapped_subints          = 0
    subints_peak_snrs       = []
    subints_integrated_snrs = []
    # Loop over each subint and calculate the appropriate metrics
    for subint in subints:
        # Get the RMS
        off_pulse_rms = N.std(subint[off_pulse_bins])
        # If there is signal in this subint...
        if off_pulse_rms != 0.0:
            # Normalize
            subint        -= N.median(subint[off_pulse_bins])
            subint        /= N.std(subint[off_pulse_bins])
            # Calculate the SNRs and store them
            peak_snr       = max(subint[on_pulse_bins])
            integrated_snr = sum(subint[on_pulse_bins])
            subints_peak_snrs.append(peak_snr)
            subints_integrated_snrs.append(integrated_snr)
            # If the SNRs are high enough, count this as a good subint
            if peak_snr > PEAK_SNR_THRESHOLD or \
               integrated_snr > INTEGRATED_SNR_THRESHOLD:
                good_subints += 1
        # If the off_pulse_rms is 0, then this subint was zapped
        else:
            zapped_subints += 1

    # It is unlikely that we would zap all subints, but just to be safe...
    if len(subints) != zapped_subints:
        frac_good_subints = good_subints/float(len(subints) - zapped_subints)
        subints_snr_std   = N.max((N.std(subints_peak_snrs),
                                   N.std(subints_integrated_snrs)))
    else:
        frac_good_subints = 0.0
        subints_snr_std   = 0.0

    return frac_good_subints,subints_snr_std


def calc_broadbandedness_metrics(subbands, on_pulse_bins, off_pulse_bins):
    """
    Calculate the fraction of sub-bands above a SNR threshold and the
    larger of the peak or integrated SNR standard deviations.

    Parameters
    ----------
    subbands : array_like
        The folded sub-bands
    on_pulse_bins : array_like
        A list of on-pulse profile bins
    off_pulse_bins : array_like
        A list of off-pulse profile bins

    Returns
    -------
    frac_good_subbands : float
        The fraction of \"good\" sub-bands
    subbands_snr_std : float
        The standard deviation of sub-bands SNRs
    """
    # Initiate some variables
    good_subbands            = 0
    zapped_subbands          = 0
    subbands_peak_snrs       = []
    subbands_integrated_snrs = []
    # Loop over each subband and calculate the appropriate metrics
    for subband in subbands:
        # Get the RMS
        off_pulse_rms = N.std(subband[off_pulse_bins])
        # If there is signal in this subband...
        if off_pulse_rms != 0.0:
            # Normalize
            subband        -= N.median(subband[off_pulse_bins])
            subband        /= N.std(subband[off_pulse_bins])
            # Calculate the SNRs and store them
            peak_snr        = max(subband[on_pulse_bins])
            integrated_snr  = sum(subband[on_pulse_bins])
            subbands_peak_snrs.append(peak_snr)
            subbands_integrated_snrs.append(integrated_snr)
            # If the SNRs are high enough, count this as a good subband
            if peak_snr > PEAK_SNR_THRESHOLD or \
               integrated_snr > INTEGRATED_SNR_THRESHOLD:
                good_subbands += 1
        # If the off_pulse_rms is 0, then this subband was zapped
        else:
            zapped_subbands += 1

    # It is unlikely that we would zap all subbands, but just to be safe...
    if len(subbands) != zapped_subbands:
        frac_good_subbands = good_subbands/float(len(subbands)-zapped_subbands)
        subbands_snr_std   = N.max((N.std(subbands_peak_snrs),
                                    N.std(subbands_integrated_snrs)))
    else:
        frac_good_subbands = 0.0
        subbands_snr_std   = 0.0

    return frac_good_subbands,subbands_snr_std


def Prepfold_Sigma(pfd):
    """
    Calculate the equivalent Gaussian significance (sigma) from the
    chi-squared as determined by prepfold.  A rescaled value is also
    returned.

    Parameters
    ----------
        pfd : class
        An instance of the prepfold.pfd class
    Returns
    -------
    names : list
        A list of the ratings names
    ratings : list
        A list of the ratings values
    """
    # The rating names
    name1 = "Prepfold_Sigma"
    name2 = "Rescaled_Prepfold_Sigma"

    # De-disperse at the best DM
    pfd.dedisperse(DM=pfd.bestdm, doppler=1)
    # Internally rotate the data cube so that is aligned with best fold values
    pfd.adjust_period()
    # Get the reduced chi^2 from the .bestprof file if it exists, otherwise
    # calculate from the best profile
    if hasattr(pfd.bestprof, "chi_sqr"):
        redchi2 = pfd.bestprof.chi_sqr
    else:
        prof    = N.sum(N.sum(pfd.profs, axis=1), axis=0)
        redchi2 = pfd.calc_redchi2(prof=prof)

    # Calculate the equivalent Gaussian significance
    dof   = pfd.proflen - 1
    chi2  = redchi2*dof
    prob  = S.special.chdtrc(dof,chi2)
    sigma = -S.special.ndtri(prob)
    if sigma > 500: sigma = 999

    # Rescale the chi^2 by the off-pulse chi^2 value (so that off-pulse chi^2
    # is equal to 1) and re-compute the equivalent Gaussian significance
    offpulse_redchi2 = pfd.estimate_offsignal_redchi2()
    chi2_scale       = 1.0/offpulse_redchi2
    rescaled_redchi2 = chi2_scale*redchi2
    rescaled_chi2    = rescaled_redchi2*dof
    rescaled_prob    = S.special.chdtrc(dof, rescaled_chi2)
    rescaled_sigma   = -S.special.ndtri(rescaled_prob)
    if rescaled_sigma > 500: rescaled_sigma = 999

    # Store the ratings values
    rating1 = sigma
    rating2 = rescaled_sigma

    return [name1,name2],[rating1,rating2]
        

# Rating function definitions
def DM0_Comparison_Ratings(pfd):
    """
    Calculate the ratio of the profile RMS, profile peak, and profile reduced
    chi-squared at DM=0 and at the best DM.

    Parameters
    ----------
    pfd : class
        An instance of the prepfold.pfd class

    Returns
    -------
    names : list
        A list of the ratings names
    ratings : list
        A list of the ratings values
    """
    # The rating names
    name1 = "DM0_RMS_Comparison"
    name2 = "DM0_Peak_Comparison"
    name3 = "DM0_ChiSqr_Comparison"
    
    # De-disperse at the best DM and calculate the RMS, peak, and reduced
    # chi-squared
    pfd.dedisperse(DM=pfd.bestdm, doppler=1)
    # Internally rotate the data cube so that is aligned with best fold values
    pfd.adjust_period()
    best_dm_prof = N.sum(N.sum(pfd.profs, axis=1), axis=0)
    best_dm_rms  = N.std(best_dm_prof)
    best_dm_peak = N.max(best_dm_prof) - N.median(best_dm_prof)
    if hasattr(pfd.bestprof, "chi_sqr"):
        best_dm_chi2 = pfd.bestprof.chi_sqr
    else:
        best_dm_chi2 = pfd.calc_redchi2(prof=best_dm_prof)

    # De-dipsers at zero DM and calculate the RMS, peak, and reduced
    # chi-squared
    pfd.dedisperse(DM=0.0, doppler=1)
    # Internally rotate the data cube so that is aligned with best fold values
    pfd.adjust_period()
    zero_dm_prof = N.sum(N.sum(pfd.profs, axis=1), axis=0)
    zero_dm_rms  = N.std(zero_dm_prof)
    zero_dm_peak = N.max(zero_dm_prof) - N.median(zero_dm_prof)
    zero_dm_chi2 = pfd.calc_redchi2(prof=zero_dm_prof)

    # Store the ratings values
    rating1 = zero_dm_rms/best_dm_rms
    rating2 = zero_dm_peak/best_dm_peak
    rating3 = zero_dm_chi2/best_dm_chi2

    return [name1,name2,name3],[rating1,rating2,rating3]


def Single_Gaussian_Ratings(pfd, debug=False):
    """
    Calculate ratings that depend on the fit of a single Gaussian profile
    component.

    Parameters
    ----------
    pfd : class
        An instance of the prepfold.pfd class
    debug : boolean
        If True, plot the profile and best-fit Gaussian

    Returns
    -------
    names : list
        A list of ratings names
    ratings : list
        A list of ratings values
    """
    # The rating names
    name1 = "Gaussian_Amplitude"
    name2 = "Gaussian_Width"
    name3 = "Gaussian_Phase"
    name4 = "Gaussian_GoF"
    name5 = "Fraction_Good_Subints"
    name6 = "Subints_SNR_StdDev"
    name7 = "Fraction_Good_Subbands"
    name8 = "Subbands_SNR_StdDev"

    # De-disperse the profile at the best DM
    pfd.dedisperse(DM=pfd.bestdm, doppler=1)
    # Internally rotate the data cube so that is aligned with best fold values
    pfd.adjust_period()
    # Get the subints, subbands, and profile
    subints  = N.sum(pfd.profs, axis=1)
    subbands = N.sum(pfd.profs, axis=0)
    profile  = N.sum(subints, axis=0)
    # Fit a gaussian to the profile
    gaussian_params, red_chi2, ngaussians = \
        PT.fit_gaussians_presto(profile, pfd.avgprof, N.sqrt(pfd.varprof),
                                1, F_STAT_THRESHOLD)

    # If we were able to fit a Gaussian...
    if ngaussians != 0:
        # Parse the Gaussian parameters
        amplitude = gaussian_params[1]
        fwhm      = gaussian_params[2]
        phase     = gaussian_params[3]
		
        # Get the on-pulse and off-pulse bins
        on_pulse_bins,off_pulse_bins = \
            PT.calc_on_pulse_region(profile,gaussian_params)
        # Get the persistence metrics
        frac_good_subints,subints_snr_std = \
            calc_persistence_metrics(subints,on_pulse_bins,off_pulse_bins)
        # Get the broadbandedness metrics
        frac_good_subbands,subbands_snr_std = \
            calc_broadbandedness_metrics(subbands,on_pulse_bins,off_pulse_bins)

        if debug:
            phases    = N.arange(pfd.proflen, dtype=N.float)/pfd.proflen
            gaussians = PT.make_gaussians_presto(gaussian_params, pfd.proflen)
            print len(on_pulse_bins)/float(pfd.proflen)
            PL.plot(phases, profile, "k")
            PL.plot(phases, gaussians, "b")
            PL.axvline(on_pulse_bins[0]/float(pfd.proflen), c="g")
            PL.axvline(on_pulse_bins[-1]/float(pfd.proflen), c="r")
            PL.show()
                    

    # If no Gaussians could be fit to the profile, store bogus metrics
    else:
        amplitude          = 0.0
        fwhm               = 0.0
        phase              = -1.0
        frac_good_subints  = -1.0
        frac_good_subbands = -1.0
        subints_snr_std    = -1.0
        subbands_snr_std   = -1.0

    # Store the ratings values
    rating1 = amplitude
    rating2 = fwhm
    rating3 = phase
    rating4 = red_chi2
    rating5 = frac_good_subints
    rating6 = subints_snr_std
    rating7 = frac_good_subbands
    rating8 = subbands_snr_std
    
    return [name1,name2,name3,name4,name5,name6,name7,name8],\
           [rating1,rating2,rating3,rating4,rating5,rating6,rating7,rating8]


def Multi_Gaussian_Ratings(pfd, debug=False):
    """
    Calculate ratings that depend on a fit of multiple Gaussian profile
    components.

    Parameters
    ----------
    pfd : class
        An instance of the prepfold.pfd class
    debug : boolean
        If True, plot the profile and best-fit Gaussians

    Returns
    -------
    names : list
        A list of ratings names
    ratings : list
        A list of ratings values
    """
    # The ratings names
    name1 = "Number_of_Gaussians"
    name2 = "Pulse_Width_Rating"

    # De-disperse the profile at the best DM
    pfd.dedisperse(DM=pfd.bestdm, doppler=1)
    # Internally rotate the data cube so that is aligned with best fold values
    pfd.adjust_period()
    # Get the best fold profile
    profile = N.sum(N.sum(pfd.profs, axis=1), axis=0)
    # Get the period from bestprof, if it exists
    if hasattr(pfd.bestprof, "p0"):
        P0  = pfd.bestprof.p0
    # Otherwise use the barycentric or topocentric folding periods
    else:
        if   pfd.bary_p1 != 0.0: P0 = pfd.bary_p1
        elif pfd.topo_p1 != 0.0: P0 = pfd.topo_p1
    # Get the center observing frequency
    f_ctr   = (pfd.hifreq + pfd.lofreq)/2
    # Fit up to MAX_GAUSSIANS Gaussians to the profile
    gaussian_params, red_chi2, ngaussians = \
        PT.fit_gaussians_presto(profile, pfd.avgprof, N.sqrt(pfd.varprof),
                             MAX_GAUSSIANS, F_STAT_THRESHOLD)

    if debug:
        phases    = N.arange(pfd.proflen, dtype=N.float)/pfd.proflen
        gaussians = PT.make_gaussians_presto(gaussian_params, pfd.proflen)
        PL.plot(phases, profile)
        PL.plot(phases, gaussians)
        PL.show()

    # Create a list for storing the FWHMs of the Gaussians
    fwhms = []
    # If we were able to fit at least one Gaussian...
    if ngaussians != 0:
        # Get and store the FWHM for each Gaussian
        for nn in xrange(ngaussians):
            fwhms.append(gaussian_params[1+3*nn+1])
        # Calculate the DM smearing time
        dm_smear_sec = PU.dm_smear(pfd.bestdm, pfd.chan_wid, f_ctr)
        # Calculate the expected minimum pulse width
        width_phase  = N.sqrt(dm_smear_sec**2 + pfd.dt**2)/P0
        # Compare the expected minimum and actual pulse widths
        width_ratio  = width_phase/min(fwhms)
    # If no Gaussians could be fit to the profile, store bogus metrics
    else:
        width_ratio = -1.0

    # Store the ratings
    rating1 = ngaussians
    rating2 = width_ratio

    return [name1,name2],[rating1,rating2]


def Phase_Wiggle_Ratings(pfd, method="GOODFRAC"):
    """
    Calculate a metric for the phase wiggle in time and frequency.

    Parameters
    ----------
    pfd : class
        An instance of the prepfold.pfd class
    method : string
        The method to use to calculate the phase wiggle metrics
        \"GOODFRAC\" - The fraction of sub-intervals/sub-bands with a wiggle
            less than some threshold
        \"WANDER\" - The total phase wander normalized by the number of
            profile bins and the number of sub-intervals/sub-bands

    Returns
    -------
    names : list
        A list of ratings names
    ratings : list
        A list of ratings values
    """
    # The ratings names
    name1 = "Phase_Wiggle_Time"
    name2 = "Phase_Wiggle_Freq"
    
    # De-disperse the profile at the best DM
    pfd.dedisperse(DM=pfd.bestdm, doppler=1)
    # Internally rotate the data cube so that is aligned with best fold values
    pfd.adjust_period()
    # Get the subints.
    subints      = N.sum(pfd.profs, axis=1)
    # Get the subbands.
    subbands     = N.sum(pfd.profs, axis=0)
    # Get the template ("best") profile
    template     = N.sum(subints, axis=0)
    # Make an array for storing a new template
    new_template = N.zeros_like(template)
    # Make an arrary for storing the offsets (in phase bins)
    bin_offsets_time  = N.empty(pfd.npart)
    bin_offsets_freq  = N.empty(pfd.npart)
    # The following loop creats a better template by removing wiggle, but
    # it does not change the actual subints
    for ii,subint in enumerate(subints):
        # Measure the phase offset
        phase_offset = PU.measure_phase_corr(subint, template)
        # The following is needed to put phase offsets on the interval
        # (-0.5,0.5]
        if phase_offset > 0.5: phase_offset -= 1.0
        # Caclulate the offset in bins
        bin_offset    = int(round(pfd.proflen*phase_offset))
        # Update the new template
        new_template += PU.rotate(subint, -bin_offset)

    # Now calculate the wiggle using the updated template
    for ii,(subint,subband) in enumerate(zip(subints,subbands)):
        phase_offset_time = PU.measure_phase_corr(subint, new_template)
        if phase_offset_time > 0.5: phase_offset_time -= 1.0
        bin_offsets_time[ii] = int(round(pfd.proflen*phase_offset_time))

        phase_offset_freq = PU.measure_phase_corr(subband, new_template)
        if phase_offset_freq > 0.5: phase_offset_freq -= 1.0
        bin_offsets_freq[ii] = int(round(pfd.proflen*phase_offset_freq))

    # Calcultae the various metrics
    good_fraction_time = sum(abs(bin_offsets_time) < \
                             PHASE_DRIFT_TOLERANCE*pfd.proflen)/ \
                             float(pfd.npart)
    total_wander_time  = sum(abs(bin_offsets_time))/(pfd.proflen*pfd.npart)
    good_fraction_freq = sum(abs(bin_offsets_freq) < \
                             PHASE_DRIFT_TOLERANCE*pfd.proflen)/ \
                         float(pfd.npart)
    total_wander_freq  = sum(abs(bin_offsets_freq))/(pfd.proflen*pfd.nsub)

    # Make the appropriate metric the rating
    if method == "GOODFRAC":
        rating1 = good_fraction_time
        rating2 = good_fraction_freq
    
    if method == "WANDER"  :
        rating1 = total_wander_time
        rating2 = total_wander_freq

    return [name1,name2],[rating1,rating2]
    

def Known_Pulsar_Rating(pfd):
    """
    Calculate the probability that a candidate is a known pulsar or a harmonic
    of a known pulsar.

    Parameters
    ----------
    pfd : class
        An instance of the prepfold.pfd class

    Returns
    -------
    names : list
        A list of ratings names
    ratings : list
        A list of ratings values
    """
    # The rating name
    name1 = "Known_Pulsar_Rating"
    # A fudge factor for error calculations
    factor = 0.3*pfd.proflen

    # Get the candidate RA and DEC (in radians)
    cand_ra  = PU.ra_to_rad(pfd.rastr)
    cand_dec = PU.dec_to_rad(pfd.decstr)

    # Get the period and folding epoch from bestprof, if it exists
    if hasattr(pfd.bestprof, "p0"):
        cand_p     = pfd.bestprof.p0
        cand_epoch = pfd.bestprof.epochi + pfd.bestprof.epochf
        # Get the candidate period error from bestprof, if it exists
        if hasattr(pfd.bestprof, "p0err"):
            cand_p_err = factor*pfd.bestprof.p0err
        # Otherwise, try to estimate it
        else:
            cand_p_err = factor*cand_p**2/(pfd.proflen*pfd.T)
    # Otherwise, use the barycentric or topocentric values
    elif pfd.bary_p1 != 0.0:
        cand_p     = pfd.bary_p1
        cand_epoch = pfd.bepoch
        cand_p_err = factor*cand_p**2/(pfd.proflen*pfd.T)
    elif pfd.topo_p1 != 0.0:
        cand_p     = pfd.topo_p1
        cand_epoch = pfd.tepoch
        cand_p_err = factor*cand_p**2/(pfd.proflen*pfd.T)
    # Get the min and max frequencies for this observation
    f_min       = pfd.subfreqs.min()
    f_max       = pfd.subfreqs.max()
    # Get the candidate's best DM
    cand_dm     = pfd.bestdm
    # Try to estimate the DM error
    cand_dm_err = factor*cand_p*f_max**2*f_min**2/ \
                  (4.15e3*pfd.proflen*(f_max**2 - f_min**2))

    # Now loop through the catalog of known pulsars and find any that are
    # close to the candidate
    nearby_psrs = [psr for psr in PSRCAT.psrs \
                   if (hasattr(psr, "ra") and hasattr(psr, "dec") and \
                   sla_dsep(cand_ra,cand_dec,psr.ra,psr.dec)<1.3*BEAM_FWHM)]

    # If there were no pulsars nearby, return a rating of 0.0 and exit
    if len(nearby_psrs) == 0:
        rating1 = 0.0
        return [name1],[rating1]
    # Otherwise, try to estimate the probability that the candidate is a
    # known pulsar or its harmonic.  This is done by calculting the difference
    # between the candidates period and DM and those any nearby known pulsars,
    # normalized by the error in the candidate values (i.e., difference in
    # "sigmas").  Calculate the probability assuming Gaussian statistics.
    else:
        # Start off assuming zero probability that the candidate is known
        max_prob     = 0.0
        for psr in nearby_psrs:
            # If possible, calculate the known pulsar period at the observing
            # epoch
            if hasattr(psr,"pepoch"):
                delta_t = (cand_epoch - psr.pepoch)*86400.0 # seconds
                if hasattr(psr,"pd") and hasattr(psr,"pdd"):
                    psr_p = psr.p + psr.pd*delta_t + 0.5*psr.pdd*delta_t**2
                elif hasattr(psr,"pd"):
                    psr_p = psr.p + psr.pd*delta_t
            else:
                psr_p = psr.p
                
            p_prob  = 0.0 # The probability that the periods are the same
            # Try all harmonic ratios with up to 16 harmonics
            for a in xrange(16):
                a += 1.0 # Since xrange starts with 0
                for b in xrange(16):
                    b += 1.0 # Since xrange starts with 0
                    # Difference between candidate and known pulsar harmonic
                    delta_p = abs(a/b*cand_p - psr_p)
                    # Calculate the equivalent Gaussian probability
                    tmp_prob  = \
                         S.special.erfc(delta_p/((a/b)**2*cand_p_err)/ \
                                        N.sqrt(2))
                    # If this is the highest probability so far, store it
                    if tmp_prob > p_prob: p_prob       = tmp_prob

            # Update the total probability
            prob      = p_prob
            # Get the difference in DM
            delta_dm  = abs(cand_dm - psr.dm)
            # Multiply the total probability by the DM probability
            prob     *= S.special.erfc(delta_dm/cand_dm_err/N.sqrt(2))
            # Update max_prob if necessary
            if prob > max_prob: max_prob = prob

        # Store the rating
        rating1 = max_prob

        return [name1],[rating1]


def Known_RFI_Freq_Diff(pfd):
    """
    The minimum absolute difference between a candidate's rotational
    frequency and those of known sources of RFI (mostly the 60 Hz power mains
    signal and its harmonics).

    Parameters
    ----------
    pfd : class
        An instance of the prepfold.pfd class

    Returns
    -------
    names : list
        A list of ratings names
    ratings : list
        A list of ratings values
    """
    # The rating name
    name1 = "Known_RFI_Freq_Diff"

    # An array of known RFI frequencies
    rfi_fs = N.array([.13,.26,.5,.9,1,1.8,2.06,2.1,2.7,3,3.6,4.15,4.7,5.7,6,6.2,6.8,7.3,8.4,8.9,9,9.4,10,11,12,15,17,18,20,23,29,30,40,60,80,90,100,120,180,200,240,300,360,400,420,480,540,600,700,800,900,1000,1100,1220,1300,1350,1400,1500,1700,1800,1900],
                     dtype=N.float)

    # Get the topocentric period if available, otherwise use the barycentric
    # period
    if pfd.topo_p1 != 0.0: p = pfd.topo_p1
    else: p = pfd.bary_p1
    # Convert to rotational frequency
    f = 1.0/p

    # Calculate the minimum absolute difference between the candidate
    # frequency and known RFI frequencies
    rating1 = min(abs(f - rfi_fs))

    return [name1],[rating1]

def UBC_AI_Rating(pfd, classifier=ubc_AI_classifier):
    """
    Run the UBC pfd artificial intelligence rater.

    Parameters
    ----------
    pfd : class
        An instance of the prepfold.pfd class

    Returns
    -------
    names : list
        A list of ratings names
    ratings : list
        A list of ratings values
    """
    # The rating name
    name1 = "UBC_AI"

    rating1 = classifier.report_score([ubc_AI.data.pfdreader(pfd.pfd_filename)])[0]
    
    return [name1],[rating1]


# A list of all the ratings
all_ratings = [Prepfold_Sigma, DM0_Comparison_Ratings,
               Single_Gaussian_Ratings, Multi_Gaussian_Ratings,
               Phase_Wiggle_Ratings, Known_Pulsar_Rating,
               Known_RFI_Freq_Diff, UBC_AI_Rating]


# Driver function definition
def rate_candidate(pfdfile, ratings=all_ratings):
    """
    Calculate all the ratings for a candidate.

    Parameters
    ----------
    pfdfile : string
        The name of a prepfold .pfd file

    Returns
    -------
    status : int
        The exit stats.  status = 0 for success, 1 for failure.
    """
    try:
        # Create an instance of the pfd class
        pfd = PF.pfd(pfdfile)
        # Open a file for storing the ratings
        outfile = open("%s.ratings"%pfdfile, "w")
        # Calculate all the ratings
        for rating in all_ratings:
            rating_names,rating_values = rating(pfd)
            # Write each rating
            for n,v in zip(rating_names,rating_values):
                outfile.write("%-25s  %9.6f\n"%(n,v))
        outfile.close()
        status = 0

    # If anything went wrong, return a non-zero exit status
    except:
        status = 1

    return status


# If running interactively
if __name__ == "__main__":
    pfdfiles = sys.argv[1:]
    for pfdfile in pfdfiles:
        status = rate_candidate(pfdfile)
        print status, pfdfile
