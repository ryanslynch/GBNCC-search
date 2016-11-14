import numpy as np
import psr_utils

import base
from sp_rating_classes import waterfall_dd

TOL = 0.02
method = "GOODFRAC"

class WiggleRater(base.BaseRater):
    short_name = "wiggle"
    long_name = "Wiggle Rating"
    description = "The fraction of subbands (0.0 - 1.0) that 'wiggle' " \
                  "in the pulse window by <= 0.02. The offset is calculated " \
                  "by cross-correlating each interval with the summed " \
                  "profile."
    version = 1

    rat_cls = waterfall_dd.WaterfallDDClass()

    def _compute_rating(self, cand):
        """Return a rating for the candidate. The rating value is the
            the fraction of subbands that deviate from the position of
            the pulse.

            Input:
                cand: A SPCandidate object to rate.

            Output:
                value: The rating value.
        """
        wdd = cand.waterfall_dd
        spd = cand.spd

        prof = wdd.get_profile()
        new_template = np.zeros_like(prof)
        bin_offsets = np.empty(spd.waterfall_nsubs)

        # The following loop creates a better template by removing wiggle, but
        # it does not change the actual subbands
        for ii, subband in enumerate(wdd.data):
            # Measure the pulse offset
            pulse_offset = psr_utils.measure_phase_corr(subband, prof, zoom=1)
           # Caclulate the offset in bins
            bin_offset = int(round(spd.waterfall_nbins*pulse_offset))
            # Update the new template
            new_template += psr_utils.rotate(subband, -bin_offset)

        # Now calculate the wiggle using the updated template
        for ii, subband in enumerate(wdd.data):
            pulse_offset = psr_utils.measure_phase_corr(subband, new_template, zoom=1)
            bin_offsets[ii] = int(round(spd.waterfall_nbins*pulse_offset))

        # Calculate the various metrics
        if method == "GOODFRAC":
            # good fraction 
            wigglescore = sum(abs(bin_offsets) < TOL*spd.waterfall_nbins)/ \
                        float(spd.waterfall_nsubs)
        elif method == "WANDER":
            # total wander
            wigglescore = sum(abs(bin_offsets))/(spd.waterfall_nbins*spd.waterfall_nsubs)
        elif method == "OFFSTD":
            # offset std
            wigglescore = bin_offsets.std()
        elif method == "OFFMAX":
            # offset max
            wigglescore = bin_offsets.max() 
        else:
            raise utils.RatingError("Unrecognized method for wiggle " \
                                    "rating (%s)" % method)

        return wigglescore
    

Rater = WiggleRater

