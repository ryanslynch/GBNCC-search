import base

class CandInfoRatingClass(base.BaseRatingClass):
    data_key = "info"

    def _compute_data(self, cand):
        info = {}
        info["dm"] = cand.dm
        info["raj_deg"] = cand.raj_deg
        info["decj_deg"] = cand.decj_deg
        info["spdfn"] = cand.spdfn
        return info
