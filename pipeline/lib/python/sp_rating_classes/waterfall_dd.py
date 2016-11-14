import copy
import psr_utils
import spd
import dataproducts

class WaterfallDDClass(spd.SpdRatingClass):
    data_key = "waterfall_dd"

    def _compute_data(self, cand):
        """Create a WaterfallDD object for the candidate.

            Input:
                cand: A ratings2.0 SPCandidate object.

            Output:
                wdd: The corresponding WaterfallDD object.
        """
        spd = cand.spd
        data = spd.data_zerodm_dedisp
        times = spd.waterfall_time_axis()
        freqs = spd.waterfall_freq_axis()
        wdd = dataproducts.WaterfallDD(data, spd.best_dm, times, freqs)
                                          
        return wdd
