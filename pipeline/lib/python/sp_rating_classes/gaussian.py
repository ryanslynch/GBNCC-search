import numpy as np

import utils
from sp_rating_classes import profile

class GaussianProfileClass(profile.ProfileClass):
    data_key = "gaussfit"
    
    def _compute_data(self, cand):
        """Fit the candidate's profile with multiple gaussian
            components and return the fit's parameters.

            Input:
                cand: A ratings2.0 SPCandidate object.

            Output:
                multigaussfit: The corresponding fit. A MultiGaussFit object.
        """
        data = utils.get_scaled_profile(cand.profile, cand.spd.varprof)

        # Initialize some starting values
        nbins = len(data)
 
        trial_params = [0.0]
 
        amplitude = max(data[(0.1*nbins):(0.4*nbins)])
        fwhm = 0.02 # full window should be 50 times estimated pulse width
        phase = 0.25 # this is where the single pulse should be placed
        trial_params.append(amplitude)
        trial_params.append(fwhm)
        trial_params.append(phase)
                
        from scipy.optimize import leastsq
        def func(params):
            #print "DEBUG: params", params
            # since this is single gaussian, params is just [offset, amp, std, phs]
            fit = utils.multigaussfit_from_paramlist(params)
            return fit.get_resids(data)

        new_params, status = leastsq(func, trial_params)
        if status not in (1,2,3,4):
            raise utils.RatingError("Status returned by " \
                                "scipy.optimize.leastsq (%d) " \
                                "indicates the fit failed!" % status)

        new_fit = utils.multigaussfit_from_paramlist(new_params)
        
        return new_fit
        
