import numpy as N
from scipy import stats
from mpfit import mpfit
import psr_utils as PU


"""
Some useful functions for fitting Gaussian components to pulse profiles.

By Ryan Lynch (Last modified 09/30/11)
"""


#==============================================================================

def rotate_profile(profile, ctr_phs):
    """
    Rotate a pulse profile so that the peak is at phase = 0.5.  This helps
    to avoid problems when trying to fit a Gaussian to a pulse that is
    split across phases 0 or 1.

    Parameters
    ----------
    profile : array_like
        The pulse profile to rotate
    ctr_phs : float
        The phase at which to center the peak profile bin
        
    Returns
    -------
    new_profile : array_like
    """
    # Get the number of profile bins
    nbins    = len(profile)
    # Find the bin containing the max profile value
    bin_max  = N.argmax(profile)
    # Determine the amount to rotate
    rotate   = int(ctr_phs*nbins - bin_max)
    # Now make a new list of bins
    new_bins = []
    for bin in xrange(nbins):
        new_bin = bin - rotate
        if new_bin >= nbins:
            new_bin -= nbins
        if new_bin < 0:
            new_bin += nbins
        new_bins.append(new_bin)

    # Return the rotated profile
    return profile[new_bins]


#==============================================================================

def make_gaussians(params, nbins):
    """
    Produce one or more Gaussian profiles with the specified parameters.

    Parameters
    ----------
    params : array_like
        The parameters of the desried Gaussian profiles.  The first entry is
        the offset from a zero baseline and the remaining entries are
        amplitudes, standard deviations, and phases, i.e.
        params = [offset, amp_1, std_1, phase_1, amp_2, std_2, phase_2, etc.]
    nbins : int
        The number of bins in the profile

    Returns
    -------
    gaussians : ndarray, shape(nbins,)
        
    """
    # Determine the number of Gaussian profiles to make
    ngaussians = (len(params) - 1)/3
    # Set up some lists to store the individual parameters
    amplitudes = []
    std_devs   = []
    phases     = []
    # Get each parameter and store it
    offset = params[0]
    for nn in xrange(ngaussians):
        amplitudes.append(params[1+3*nn])
        std_devs.append(params[1+3*nn+1])
        phases.append(params[1+3*nn+2])
        
    # Create an array of phase bins going from 0 --> 1
    bins     = N.arange(nbins, dtype=N.float)/nbins
    # Create an array for the Gaussian profile
    gaussians = N.zeros_like(bins) + offset
    # Add each individual Gaussian to the full profile
    for a,s,p in zip(amplitudes, std_devs, phases):
        gaussians += a*N.sqrt(2*N.pi*s*s)*stats.norm.pdf(bins, loc=p, scale=s)

    return gaussians

#==============================================================================

def make_gaussians_presto(params, nbins):
    """
    Produce one or more Gaussian profiles with the specified parameters.

    Parameters
    ----------
    params : array_like
        The parameters of the desried Gaussian profiles.  The first entry is
        the offset from a zero baseline and the remaining entries are
        amplitudes, full widths at half max, and phases, i.e.
        params = [offset, amp_1, fwhm_1, phase_1, amp_2, fwhm_2, phase_2, etc.]
    nbins : int
        The number of bins in the profile

    Returns
    -------
    gaussians : ndarray, shape(nbins,)
        
    """
    # Determine the number of Gaussian profiles to make
    ngaussians = (len(params) - 1)/3
    # Set up some lists to store the individual parameters
    amplitudes = []
    fwhms      = []
    phases     = []
    # Get each parameter and store it
    offset = params[0]
    for nn in xrange(ngaussians):
        amplitudes.append(params[1+3*nn])
        fwhms.append(params[1+3*nn+1])
        phases.append(params[1+3*nn+2])
        
    # Create an array for the Gaussian profile
    gaussians = N.zeros(nbins) + offset
    # Add each individual Gaussian to the full profile
    for a,w,p in zip(amplitudes, fwhms, phases):
        gaussians += a*w/2*N.sqrt(N.pi/N.log(2))*PU.gaussian_profile(nbins,p,w)

    return gaussians


#==============================================================================

def make_vonmises(params, nbins):
    """
    Produce one or more von Mises profiles with the specified parameters.

    Parameters
    ----------
    params : array_like
        The parameters of the desried von Mises profiles.  The first entry is
        the offset from a zero baseline and the remaining entries are
        amplitudes, concentrations, and phases, i.e.
        params = [offset, amp_1, conc_1, phase_1, amp_2, conc_2, phase_2, etc.]
    nbins : int
        The number of bins in the profile

    Returns
    -------
    vonmises : ndarray, shape(nbins,)
        
    """
    # Determine the number of von Mises profiles to make
    nvonmises = (len(params) - 1)/3
    # Set up some lists to store the individual parameters
    amplitudes     = []
    concentrations = []
    locations      = []
    # Get each parameter and store it
    offset = params[0]
    for nn in xrange(nvonmises):
        amplitudes.append(params[1+3*nn])
        concentrations.append(params[1+3*nn+1])
        # Need to convert location from a phase (0 --> 1) to an angle
        # (-pi --> pi)
        location = 2*N.pi*(params[1+3*nn+2] - 0.5)
        locations.append(location)
        
    # Create an array of phase bins going from -pi --> pi
    bins     = N.arange(-N.pi, N.pi, 2*N.pi/nbins)
    # Create an array for the von Mises profile
    vonmises = N.zeros_like(bins) + offset
    # Add each individual von Mises to the full profile
    for a,c,l in zip(amplitudes, concentrations, locations):
        vonmises += a*stats.vonmises.pdf(bins, c, loc=l)/ \
                    stats.vonmises.pdf(l, c, loc=l)

    return vonmises


#==============================================================================

def fit_gaussians(data, data_avg, data_std, max_gaussians, F_stat_threshold):
    """
    Fit one or more Gaussian profiles to data using mpfit.

    Parameters
    ----------
    data : ndarray
        The data to fit
    data_avg : float
        The average of the data (used for normalization)
    data_std : float
        The standard deviation of the data (used for normalization)
    max_gaussians : int
        The maximum number of Gaussian profile components to fit
    F_stat_threshold : float
        The threshold probability that a profile component is not required for
        rejecting more components

    Returns
    -------
    params : list
        The best fit Gaussian parameters.  The first entry is the offset from
        zero baseline and the remaining entries are amplitudes, standard
        deviations, and phases, i.e.
        params = [offset, amp_1, std_1, phase_1, amp_2, std_2, phase_2, etc.]
    ngaussians : int
        The number of Gaussian profile components used to fit the data
    red_chi2 : float
        The reduced chi^2 of the fit
    """
    # Normalize the data to have average = 0 and std_dev = 1 (It would be best
    # to do so based on values for the off-pulse region, but we don't know what
    # that is...)
    data -= data_avg
    data /= data_std

    # Initialize some starting values
    nbins      = len(data)
    ngaussians = 0
    # After normalization the first parameter (offset) should be close to zero
    prev_params = [0.0]
    # Nothing fit yet, so residuals are just the data values
    prev_residuals = data - N.zeros_like(data) 
    # No need to normalize chi^2 by variance since we already did that to the
    # data
    prev_chi2  = sum(prev_residuals*prev_residuals)
    prev_dof   = nbins
    fit        = True

    # We will now start fitting Gaussian profile components until the
    # additional components are no longer statistically needed to improve the
    # fit.  The starting parameter guesses for each new component will come
    # from the highest remaining residual and from the previous best-fit values
    # for previous components
    while fit:
        ngaussians  += 1
        # Update values based on results of previous run
        trial_params = list(prev_params)

        # Guess the parameters for the next profile component
        amplitude = max(prev_residuals)
        # Base std_dev on stats.norm normalization
        std_dev   = 1/(N.sqrt(2*N.pi)*amplitude)
        phase     = N.argmax(prev_residuals)/float(nbins)
        trial_params.append(amplitude)
        trial_params.append(std_dev)
        trial_params.append(phase)
        # params_dict is used by mpfit to get initial values and constraints on
        # parameters
        params_dict = []
        for ii,param in enumerate(trial_params):
            if ii == 0:
                # The first parameter is the offset, which can be negative and
                # should be allowed to vary more
                params_dict.append({"value"  : param,
                                    "fixed"  : False,
                                    "limited": [False,False],
                                    "limits" : [0.0,0.0]})
            else:
                # Limits are set assuming that our initial guesses were correct
                # to within 25%...
                params_dict.append({"value"  : param,
                                    "fixed"  : False,
                                    "limited": [True,True],
                                    "limits" : [0.25*param,1.75*param]})

        # Define the fitting function for mpfit
        def func(params, fjac=None, errs=None):
            # Return values are [status, residuals]
            return [0, data - make_gaussians(params, nbins)]

        # Now fit
        mpfit_out     = mpfit(func, parinfo=params_dict, quiet=True)
        # Store the new best-fit parameters
        new_params    = mpfit_out.params
        # Calculate the new residuals and statistics
        new_residuals = data - make_gaussians(new_params, nbins)
        new_chi2      = mpfit_out.fnorm
        new_dof       = nbins - len(new_params) # Degrees-of-freedom
        # Calculate the F-statistic for the fit, i.e. the probability that the
        # additional profile component is /not/ required by the data
        F_stat        = PU.Ftest(prev_chi2, prev_dof, new_chi2, new_dof)

        # If the F-test probability is greater than some threshold, then the
        # additional Gaussian did not significantly improve the fit and we
        # should stop.  The nan test is needed because if the fit is /worse/
        # then Ftest doesn't return a valid number.  Also stop if we reach
        # the maximum number of Gaussian profile components
        if F_stat > F_stat_threshold or N.isnan(F_stat) \
               or ngaussians > max_gaussians:
            fit    = False
        # Otherwise, keep fitting and update the parameters for the next pass
        else:
            fit            = True
            prev_params    = new_params
            prev_residuals = new_residuals
            prev_chi2      = new_chi2
            prev_dof       = new_dof

    # We stop when a fit is no longer needed, so we have to return the values
    # from the /previous/ run (otherwise we return the unneeded fit)
    return prev_params, prev_chi2/prev_dof, ngaussians-1

#==============================================================================

def fit_gaussians_presto(data, data_avg, data_std, max_gaussians,
                         F_stat_threshold):
    """
    Fit one or more Gaussian profiles to data using mpfit.

    Parameters
    ----------
    data : ndarray
        The data to fit
    data_avg : float
        The average of the data (used for normalization)
    data_std : float
        The standard deviation of the data (used for normalization)
    max_gaussians : int
        The maximum number of Gaussian profile components to fit
    F_stat_threshold : float
        The threshold probability that a profile component is not required for
        rejecting more components

    Returns
    -------
    params : list
        The best fit Gaussian parameters.  The first entry is the offset from
        zero baseline and the remaining entries are amplitudes, full widths at
        half max, and phases, i.e.
        params = [offset, amp_1, fwhm_1, phase_1, amp_2, fwhm_2, phase_2, etc.]
    ngaussians : int
        The number of Gaussian profile components used to fit the data
    red_chi2 : float
        The reduced chi^2 of the fit
    """
    # Normalize the data to have average = 0 and std_dev = 1 (It would be best
    # to do so based on values for the off-pulse region, but we don't know what
    # that is...)
    data -= data_avg
    data /= data_std

    # Initialize some starting values
    nbins      = len(data)
    ngaussians = 0
    # After normalization the first parameter (offset) should be close to zero
    prev_params = [0.0]
    # Nothing fit yet, so residuals are just the data values
    prev_residuals = data - N.zeros_like(data) 
    # No need to normalize chi^2 by variance since we already did that to the
    # data
    prev_chi2  = sum(prev_residuals*prev_residuals)
    prev_dof   = nbins
    fit        = True

    # We will now start fitting Gaussian profile components until the
    # additional components are no longer statistically needed to improve the
    # fit.  The starting parameter guesses for each new component will come
    # from the highest remaining residual and from the previous best-fit values
    # for previous components
    while fit:
        ngaussians  += 1
        # Update values based on results of previous run
        trial_params = list(prev_params)

        # Guess the parameters for the next profile component
        amplitude = max(prev_residuals)
        # Stat with a FWHM guess of 0.1
        fwhm      = 2*N.sqrt(2*N.log(2))/(N.sqrt(2*N.pi)*amplitude)
        phase     = N.argmax(prev_residuals)/float(nbins)
        trial_params.append(amplitude)
        trial_params.append(fwhm)
        trial_params.append(phase)
        # params_dict is used by mpfit to get initial values and constraints on
        # parameters
        params_dict = []
        for ii,param in enumerate(trial_params):
            if ii == 0:
                # The first parameter is the offset, which can be negative and
                # should be allowed to vary more
                params_dict.append({"value"  : param,
                                    "fixed"  : False,
                                    "limited": [False,False],
                                    "limits" : [0.0,0.0]})
            elif (ii - 1)%3 == 1:
                # This is the FWHM, and is allowed to vary between
                # 0.01 and 1.0
                params_dict.append({"value"  : param,
                                    "fixed"  : False,
                                    "limited": [True,True],
                                    "limits" : [0.01,1.0]})
            else:
                # Limits are set assuming that our initial guesses were correct
                # to within 25%...
                params_dict.append({"value"  : param,
                                    "fixed"  : False,
                                    "limited": [True,True],
                                    "limits" : [0.25*param,1.75*param]})

        # Define the fitting function for mpfit
        def func(params, fjac=None, errs=None):
            # Return values are [status, residuals]
            return [0, data - make_gaussians_presto(params, nbins)]

        # Now fit
        mpfit_out     = mpfit(func, parinfo=params_dict, quiet=True)
        # Store the new best-fit parameters
        new_params    = mpfit_out.params
        # Calculate the new residuals and statistics
        new_residuals = data - make_gaussians_presto(new_params, nbins)
        new_chi2      = mpfit_out.fnorm
        new_dof       = nbins - len(new_params) # Degrees-of-freedom
        # Calculate the F-statistic for the fit, i.e. the probability that the
        # additional profile component is /not/ required by the data
        F_stat        = PU.Ftest(prev_chi2, prev_dof, new_chi2, new_dof)

        # If the F-test probability is greater than some threshold, then the
        # additional Gaussian did not significantly improve the fit and we
        # should stop.  The nan test is needed because if the fit is /worse/
        # then Ftest doesn't return a valid number.  Also stop if we reach
        # the maximum number of Gaussian profile components
        if F_stat > F_stat_threshold or N.isnan(F_stat) \
               or ngaussians > max_gaussians:
            fit    = False
        # Otherwise, keep fitting and update the parameters for the next pass
        else:
            fit            = True
            prev_params    = new_params
            prev_residuals = new_residuals
            prev_chi2      = new_chi2
            prev_dof       = new_dof

    # We stop when a fit is no longer needed, so we have to return the values
    # from the /previous/ run (otherwise we return the unneeded fit)
    return prev_params, prev_chi2/prev_dof, ngaussians-1


#==============================================================================

def fit_vonmises(data, data_avg, data_std, max_vonmises, F_stat_threshold):
    """
    Fit one or more von Mises profiles to data using mpfit.

    Parameters
    ----------
    data : ndarray
        The data to fit
    data_avg : float
        The average of the data (used for normalization)
    data_std : float
        The standard deviation of the data (used for normalization)
    max_vonmises : int
        The maximum number of von Mises profile components to fit
    F_stat_threshold : float
        The threshold probability that a profile component is not required for
        rejecting more components

    Returns
    -------
    params : list
        The best fit von Mises parameters.  The first entry is the offset from
        zero baseline and the remaining entries are amplitudes, concentrations,
        and phases, i.e.
        params = [offset, amp_1, conc_1, phase_1, amp_2, conc_2, phase_2, etc.]
    nvonmises : int
        The number of von Mises profile components used to fit the data
    red_chi2 : float
        The reduced chi^2 of the fit
    """
    # Normalize the data to have average = 0 and std_dev = 1 (It would be best
    # to do so based on values for the off-pulse region, but we don't know what
    # that is...)
    data -= data_avg
    data /= data_std

    # Initialize some starting values
    nbins      = len(data)
    nvonmises = 0
    # After normalization the first parameter (offset) should be close to zero
    prev_params = [0.0]
    # Nothing fit yet, so residuals are just the data values
    prev_residuals = data - N.zeros_like(data) 
    # No need to normalize chi^2 by variance since we already did that to the
    # data
    prev_chi2  = sum(prev_residuals*prev_residuals)
    prev_dof   = nbins
    fit        = True

    # We will now start fitting von Mises profile components until the
    # additional components are no longer statistically needed to improve the
    # fit.  The starting parameter guesses for each new component will come
    # from the highest remaining residual and from the previous best-fit values
    # for previous components
    while fit:
        nvonmises  += 1
        # Update values based on results of previous run
        trial_params = list(prev_params)

        # Guess the parameters for the next profile component
        amplitude = max(prev_residuals)
        # Assume a concentration appropriate for a FWHM of 0.075 (in phase)
        concentration = 25.0
        location      = N.argmax(prev_residuals)/float(nbins)
        trial_params.append(amplitude)
        trial_params.append(concentration)
        trial_params.append(location)
        # params_dict is used by mpfit to get initial values and constraints on
        # parameters
        params_dict = []
        for ii,param in enumerate(trial_params):
            if ii == 0:
                # The first parameter is the offset, which can be negative and
                # should be allowed to vary more
                params_dict.append({"value"  : param,
                                    "fixed"  : False,
                                    "limited": [False,False],
                                    "limits" : [0.0,0.0]})
            elif (ii - 1)%3 == 1:
                # This is the concentration, and is allowed to vary between
                # values appropriate for FWHMs of ~0.015 to ~0.5
                params_dict.append({"value"  : param,
                                    "fixed"  : False,
                                    "limited": [True,True],
                                    "limits" : [0.1,600.0]})
            else:
                # Limits are set assuming that our initial guesses were correct
                # to within 25%...
                params_dict.append({"value"  : param,
                                    "fixed"  : False,
                                    "limited": [True,True],
                                    "limits" : [0.25*param,1.75*param]})

        # Define the fitting function for mpfit
        def func(params, fjac=None, errs=None):
            # Return values are [status, residuals]
            return [0, data - make_vonmises(params, nbins)]

        # Now fit
        mpfit_out     = mpfit(func, parinfo=params_dict, quiet=True)
        # Store the new best-fit parameters
        new_params    = mpfit_out.params
        # Calculate the new residuals and statistics
        new_residuals = data - make_vonmises(new_params, nbins)
        new_chi2      = mpfit_out.fnorm
        new_dof       = nbins - len(new_params) # Degrees-of-freedom
        # Calculate the F-statistic for the fit, i.e. the probability that the
        # additional profile component is /not/ required by the data
        F_stat        = PU.Ftest(prev_chi2, prev_dof, new_chi2, new_dof)

        # If the F-test probability is greater than some threshold, then the
        # additional Gaussian did not significantly improve the fit and we
        # should stop.  The nan test is needed because if the fit is /worse/
        # then Ftest doesn't return a valid number.  Also stop if we reach
        # the maximum number of Gaussian profile components
        if F_stat > F_stat_threshold or N.isnan(F_stat) \
               or nvonmises > max_vonmises:
            fit    = False
        # Otherwise, keep fitting and update the parameters for the next pass
        else:
            fit            = True
            prev_params    = new_params
            prev_residuals = new_residuals
            prev_chi2      = new_chi2
            prev_dof       = new_dof

    # We stop when a fit is no longer needed, so we have to return the values
    # from the /previous/ run (otherwise we return the unneeded fit)
    return prev_params, prev_chi2/prev_dof, nvonmises-1


#==============================================================================
def calc_on_pulse_region(profile, params):
    """
    Calculate the on-pulse and off-pulse regions.  The on-pulse region is
    defined as beginning and ending where the pulse amplitude equals the
    off-pulse RMS noise level.

    Parameters
    ----------
    profile : array_like
        The best-fold pulse profile
    gaussian_params : array_like
        The parameters of the best-fit Gaussian profile.  The first entry is
        the offset from a zero baseline and the remaining entries are
        amplitude, FWHM, and phase, i.e. params = [offset, amp, FWHM, phase].

    Returns
    -------
    on_pulse_bins : list
        The profile bins for the on-pulse region
    off_pulse_bins : list
        The profile bins for the off-pulse region
    """
    # Parse the Gaussian parameters and get the number of profile bins
    amplitude = params[1]
    fwhm      = params[2]
    phase     = params[3]
    nbins     = len(profile)
    # Make a first guess at the on-pulse region
    on_pulse_start = phase - 0.5*fwhm
    on_pulse_end   = phase + 0.5*fwhm
    # Adjust start/end phases if they are outside the interval [0,1)
    if on_pulse_start <  0.0: on_pulse_start += 1.0
    if on_pulse_end   >= 1.0: on_pulse_end   -= 1.0
    # Make the list of on-pulse and off-pulse bins
    start_bin = int(round(nbins*on_pulse_start))
    end_bin   = int(round(nbins*on_pulse_end))
    if start_bin > end_bin:
        on_pulse_bins = range(start_bin, nbins) + range(0, end_bin+1)
    else:
        on_pulse_bins = range(start_bin, end_bin+1)
    off_pulse_bins = list(set(range(nbins)) - set(on_pulse_bins))
    """
    off_pulse_rms  = N.std(profile[off_pulse_bins])
    # Now that we know the approximate off-pulse RMS, re-calculate the 
    # on-pulse region so that it begins and ends when the pulse amplitude
    # is equal to the noise level.  Re-calculate the off-pulse region, too.
    on_pulse_start = phase - 0.5*fwhm*N.sqrt(N.log(amplitude/off_pulse_rms)/\
                                             N.log(2.0))
    on_pulse_end   = phase + 0.5*fwhm*N.sqrt(N.log(amplitude/off_pulse_rms)/\
                                             N.log(2.0))

    print on_pulse_start,on_pulse_end

    if on_pulse_start <  0.0: on_pulse_start += 1.0
    if on_pulse_end   >= 1.0: on_pulse_end   -= 1.0

    start_bin = int(round(nbins*on_pulse_start))
    end_bin   = int(round(nbins*on_pulse_end))
    if start_bin > end_bin:
        on_pulse_bins = range(start_bin, nbins) + range(0, end_bin+1)
    else:
        on_pulse_bins = range(start_bin, end_bin+1)
    off_pulse_bins = list(set(range(nbins)) - set(on_pulse_bins))
    """
    return on_pulse_bins, off_pulse_bins
