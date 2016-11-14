import os.path
import numpy as np

import base
from sp_rating_classes import spd
#from sp_utils import ddm_response
from singlepulse.read_spd import ddm_response

# Initialise constants
KNOWNPSR_FILENM = os.path.join(os.path.split(__file__)[0], "../knownpulsars.csv")
BEAM_FWHM_ARCMIN = 3.35 
MOCK_BAND = [1214., 1537.]

def gaussian_response(sep, sigma):
    return np.exp(-0.5*pow(sep/sigma, 2))

class KnownPulsarRater(base.BaseRater):
    short_name = "knownpsr"
    long_name = "Known Pulsar Rating"
    description = "Evaluate how similar the position and DM are to a " \
                    "known pulsar.  The value is between 0 and 1, with " \
                    "values closer to 1 indicating similarity to a known" \
                    "pulsar."
    version = 2

    rat_cls = spd.SpdRatingClass()

    def _setup(self):
        """A setup method to be called when the Rater is initialised
        
            Inputs:
                None

            Outputs:
                None
        """

        known_pulsars = np.recfromcsv(KNOWNPSR_FILENM, delimiter=";",\
          comments='#', usecols=(1,2,3,4,5))
          #np.dtype=(np.str, np.float, np.float, np.float, np.float))
          #dtype=(str, float, float, float, float))

        self.known_names = known_pulsars['name']
        self.known_ras = known_pulsars['rajd']
        self.known_decs = known_pulsars['decjd']
        self.known_dms = known_pulsars['dm']

    def _compute_rating(self, cand):
        """Return a rating for the candidate. The rating value encodes 
            how close the candidate's position and DM are to that of a
            known pulsar.

            Input:
                cand: An SPCandidate object to rate.

            Output:
                value: The rating value.
        """
        dm = cand.info['dm']
        ra = cand.info['raj_deg']
        dec = cand.info['decj_deg']
        pulsewidth = cand.spd.pulsewidth_seconds

        offsets_x = (self.known_ras - ra)*np.cos(0.5*(self.known_decs + dec))
        offsets_y = self.known_decs - dec
        offsets_deg = np.sqrt(offsets_x**2 + offsets_y**2)
        offsets_arcmin = offsets_deg * 60.

        dm_diffs = np.abs(self.known_dms - dm)

        # If the signal is at exactly the position and DM of a known pulsar,
        # this will return 1.0.  The beam response drops off as a gaussian,
        # the dDM response drops off according to Cordes & McLaughlin 2003
        # using the observed pulse width
        gresp = gaussian_response(offsets_arcmin, \
          BEAM_FWHM_ARCMIN / (2. * np.sqrt(2.*np.log(2))))
        dresp = ddm_response(dm_diffs, pulsewidth*1000., MOCK_BAND)

        # Picking the value that presumably corresponds to the best pulsar match
        return np.max(np.sqrt(gresp*dresp))

Rater = KnownPulsarRater
