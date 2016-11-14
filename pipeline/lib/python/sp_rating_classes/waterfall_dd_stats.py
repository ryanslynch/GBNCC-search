import numpy as np

import dataproducts
import gaussian
import waterfall_dd

class WaterfallDDStats(gaussian.GaussianProfileClass, waterfall_dd.WaterfallDDClass):
    data_key = "waterfall_dd_stats"

    def _compute_data(self, cand):
        """Compute statistics for waterfall plot and return an object storing them.

            Input:
                cand: A Ratings2.0 SPCandidate object.

            Output:
                subband_stats: The resulting PulseWindowStats object.
        """
        gauss = cand.gaussfit
        wdd = cand.waterfall_dd

        onpulse_region = gauss.get_onpulse_region(wdd.nbin)
        offpulse_region = np.bitwise_not(onpulse_region)

        zapped_profs = np.zeros(wdd.nchan, dtype=bool)
        snrs = np.empty(wdd.nchan)
        peak_snrs = np.empty(wdd.nchan)
        corr_coefs = np.empty(wdd.nchan)
        gaussprof = gauss.make_gaussians(wdd.nbin)
        for ichan in np.arange(wdd.nchan):
            profile = wdd.data[ichan,:].copy()
            offpulse = profile[offpulse_region]
            
            stddev = np.std(offpulse)
            if stddev == 0:
                zapped_profs[ichan] = True
                continue

            # Scale profile
            profile -= np.mean(offpulse)
            profile /= stddev
            
            # Calculate snrs
            onpulse = profile[onpulse_region] # Get on-pulse now so it is scaled
            snrs[ichan] = np.sum(onpulse)
            peak_snrs[ichan] = np.mean(onpulse)
            
            # Determine correlation coeff
            corr_coefs[ichan] = np.corrcoef(gaussprof, profile)[0][1]
        
        snrs = np.ma.masked_array(snrs, zapped_profs)
        peak_snrs = np.ma.masked_array(peak_snrs, zapped_profs)
        corr_coefs = np.ma.masked_array(corr_coefs, zapped_profs)
        pwstats = dataproducts.PulseWindowStats(snrs, peak_snrs, corr_coefs)
        return pwstats
