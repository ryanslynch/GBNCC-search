import os.path
import numpy as np
from astropy.coordinates import SkyCoord
from scipy.spatial import cKDTree

import base
from sp_rating_classes import spd

class MaxDMRatioRater(base.BaseRater):
    short_name = "max_dm_ratio"
    long_name = "Maximum DM Ratio Rating"
    description = "Return the ratio of the candidate DM to the maximum "\
                  "DM in the candidate direction according to the NE2001 "\
                  "model for the Milky Way. A value greater than 1 indicates "\
                  "a DM higher than the nominal maximum."
    version = 1

    rat_cls = spd.SpdRatingClass()

    def _setup(self):
        """A setup method to be called when the Rater is initialised
        
            Inputs:
                None

            Outputs:
                None
        """
          
        self.cl = 0.  
        self.cb = 0. 
        self.max_DM = 1595.88 #pc/cc
        
    def _compute_rating(self, cand):
        """Return a rating for the candidate. The rating value is a ratio
            of the candidate DM to the max DM along the line of sight in the
            Galaxy.

            Input:
                cand: An SPCandidate object to rate.

            Output:
                value: The rating value.
        """

        dm = cand.info['dm']
        ra = cand.info['raj_deg']
        dec = cand.info['decj_deg']

        pos = SkyCoord(ra=ra, dec=dec, unit='degree').galactic
        cl_cand = pos.l.radian
        cb_cand = pos.b.radian

        #Call to NE2001 only if max DM has not been computed for that position previously
        if self.cl != cl_cand or self.cb != cb_cand: 
            import NE2001
            dist = 55. #Dist > Dist to LMC 
            DM = 0. #Dummy Value
            self.max_DM = NE2001.dmdsm(cl_cand,cb_cand,-1,DM,dist)[0]
            self.cl = cl_cand
            self.cb = cb_cand

        return dm / self.max_DM

Rater = MaxDMRatioRater
