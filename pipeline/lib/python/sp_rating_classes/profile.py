import waterfall_dd

class ProfileClass(waterfall_dd.WaterfallDDClass):
    data_key = "profile"

    def _compute_data(self, cand):
        """Create the dedispersed optimal profile for the candidate.

            Input:
                cand: A ratings2.0 SPCandidate object.

            Output:
                prof: The corresponding Profile object.
        """
        spd = cand.spd
        wdd = cand.waterfall_dd
        return wdd.get_profile()
