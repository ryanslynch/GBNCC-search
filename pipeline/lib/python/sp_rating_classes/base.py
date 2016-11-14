class BaseRatingClass(object):
    data_key = "dummy"

    def add_data(self, cand):
        if not self.has_data(cand):
            parents = [cl for cl in self.__class__.__bases__ if cl != self.__class__]
            #print self.__class__ # Debug
            #print parents # Debug
            for parent in parents:
                if issubclass(parent, BaseRatingClass):
                    
                    parent().add_data(cand)
            data = self._compute_data(cand)
            cand.add_to_cache(self.data_key, data)

    def has_data(self, cand):
        return cand.is_in_cache(self.data_key)

    def get_data(self, cand):
        self.add_data(cand)
        return cand.get_from_cache(self.data_key)

    def _compute_data(self, cand):
        return None
        
