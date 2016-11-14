import base
from sp_rating_classes import gaussian

class GaussianFWHMRater(base.BaseRater):
    short_name = "gaussfwhm"
    long_name = "Gaussian FWHM Rating"
    description = "The full-width at half-maximum of a single-Gaussian " \
                  "fit to the profile."
    version = 1

    rat_cls = gaussian.GaussianProfileClass()

    def _compute_rating(self, cand):
        """Return a rating for the candidate. The rating value is the
            FWHM of the gaussian function fit to the candidate's 
            profile.
        
            Input:
                cand: An SPCandidate object to rate.

            Output:
                value: The rating value.
        """
        sgauss = cand.gaussfit
        ncomp = len(sgauss.components)
        if ncomp == 1:
            fwhm = sgauss.components[0].fwhm
        elif ncomp  == 0:
            fwhm = 0.0
        else:
            raise utils.RatingError("Bad number of components for single " \
                                    "gaussian fit (%d)" % ncomp)
        return fwhm


Rater = GaussianFWHMRater
