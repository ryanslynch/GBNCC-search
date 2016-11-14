import warnings

import utils

registered_raters = ["peak_over_rms",\
                     "wiggle",\
                     "gaussian_amplitude",\
                     "gaussian_goodness",\
                     "gaussian_fwhm",\
                     "frac_good_subbands",\
                     "subbands_snrstd",\
                     "known_pulsar",\
                     "max_dm_ratio",\
                     "ubc_spd_ai"]

__all__ = registered_raters

for ii in reversed(range(len(registered_raters))):
    rater_name = registered_raters[ii]
    try:
        __import__(rater_name, globals())
    except:
        warnings.warn("The rater '%s' could not be loaded!" % rater_name, \
                utils.RaterLoadWarning)
        registered_raters.pop(ii)
