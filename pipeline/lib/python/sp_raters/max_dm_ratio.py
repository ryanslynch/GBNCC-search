import os.path
import numpy as np
from astropy.coordinates import SkyCoord
from scipy.spatial import cKDTree

import base
from sp_rating_classes import spd

# Initialise constants
NE2001_FILENM = os.path.join(os.path.split(__file__)[0], "../NE2001_grid.npz")

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

        dat = np.load(NE2001_FILENM)
        self.lb_pairs = dat['lb_pairs']
        self.max_DM = dat['max_DM']
        self.search_tree = cKDTree(self.lb_pairs)

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
        cl = pos.l.degree
        cb = pos.b.degree

        grid_distance, index = \
          self.search_tree.query([cl, cb], distance_upper_bound=1)

        return dm / self.max_DM[index]

Rater = MaxDMRatioRater
