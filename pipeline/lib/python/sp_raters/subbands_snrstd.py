import base
from sp_rating_classes import waterfall_dd_stats

class SubbandsSNRStdev(base.BaseRater):
    short_name = "subbands_snrstd"
    long_name = "Subband SNR Std Dev"
    description = "Determines the standard deviation of the subband " \
                  "signal-to-noise ratios."
    version = 1
    rat_cls = waterfall_dd_stats.WaterfallDDStats()

    def _compute_rating(self, cand):
        """Return a rating for the candidate. The rating value is the
            standard deviation of the sub-band SNRs.

            Input:
                cand: An SPCandidate object to rate.

            Output:
                value: The rating value.
        """
        chanstats = cand.waterfall_dd_stats
        return max((chanstats.get_snr_stddev(), chanstats.get_peak_snr_stddev()))


Rater = SubbandsSNRStdev 
