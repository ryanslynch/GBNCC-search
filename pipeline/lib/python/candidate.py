"""
A module part of the ratings2.0 package.
The module defines a candidate object that is rated by other
modules of this package.
Patrick Lazarus, Dec 15. 2011 - mid-flight
"""
import prepfold
import psr_utils
#import sp_utils
import singlepulse.read_spd as read_spd

class Candidate(object):
    """
    A candidate object wrapping a PRESTO prepfold.pfd object
    adding extra information and methods useful for automatic
    candidate rating as used in surveys.
    """
    def __init__(self, topo_period, bary_period, dm, raj_deg, decj_deg, pfdfn):
        """Constructor for Candidate objects
            
            Inputs:
                topo_period: Topocentric period (in ms)
                bary_period: Barycentric period (in ms)
                dm: Dispersion measure (in pc cm^-3)
                raj_deg: Right Ascension (J2000 in degrees)
                decj_deg: Declination (J2000 in degrees)
                pfdfn: The *.pfd file for this candidate.
 
            Output:
                cand: Candidate object
        """
        self.topo_period = topo_period
        self.bary_period = bary_period
        self.dm = dm
        self.raj_deg = raj_deg
        self.decj_deg = decj_deg
        self.pfdfn = pfdfn
        self.rating_values = []
        self.cache = {}

    def add_rating(self, ratval):
        self.rating_values.append(ratval)

    def get_ratings_string(self, sep=('-'*45+'\n')):
        return sep.join([str(rv)+'\n' for rv in self.rating_values])
    
    def get_ratings_overview(self):
        return '\n'.join([rv.get_short_string() for rv in self.rating_values])

    def write_ratings_to_file(self, fn=None, *args, **kwargs):
        """Write candidate's ratings to file.
            
            Inputs:
                fn: File name to write ratings to.
                    (Default: generate the file name based on
                            the candidate's PFD file's name.)
                **Additional arguments are passed to 'get_ratings_string(...)'
            Outputs:
                fn: The output file name.
        """
        if fn is None:
            fn = self.pfdfn+".rat"
        f = open(fn, 'w')
        f.write(self.get_ratings_string(*args, **kwargs))
        f.close()
        return fn

    def add_to_cache(self, key, val):
        self.cache[key] = val
        
    def get_from_cache(self, key):
        return self.cache[key]
        
    def is_in_cache(self, key):
        return (key in self.cache)

    def clear_cache(self):
        self.cache.clear()
        self.cache = {}


class SPCandidate(object):
    """
    Analagous to Candidate class, but for single pulse candidates.
    """
    def __init__(self, dm, raj_deg, decj_deg, spdfn):
        """Constructor for SPCandidate objects
            
            Inputs:
                dm: Dispersion measure (in pc cm^-3)
                raj_deg: Right Ascension (J2000 in degrees)
                decj_deg: Declination (J2000 in degrees)
                spdfn: The *.spd file for this candidate.
 
            Output:
                cand: SPCandidate object
        """
        self.dm = dm
        self.raj_deg = raj_deg
        self.decj_deg = decj_deg
        self.spdfn = spdfn
        self.rating_values = []

    def add_rating(self, ratval):
        self.rating_values.append(ratval)

    def get_ratings_string(self, sep=('-'*45+'\n')):
        return sep.join([str(rv)+'\n' for rv in self.rating_values])
    
    def get_ratings_overview(self):
        return '\n'.join([rv.get_short_string() for rv in self.rating_values])

    def write_ratings_to_file(self, fn=None, *args, **kwargs):
        """Write candidate's ratings to file.
            
            Inputs:
                fn: File name to write ratings to.
                    (Default: generate the file name based on
                            the candidate's SPD file's name.)
                **Additional arguments are passed to 'get_ratings_string(...)'
            Outputs:
                fn: The output file name.
        """
        if fn is None:
            fn = self.spdfn+".rat"
        f = open(fn, 'w')
        f.write(self.get_ratings_string(*args, **kwargs))
        f.close()
        return fn

    def add_to_cache(self, key, val):
        setattr(self, key, val)
        
    def get_from_cache(self, key):
        return getattr(self, key)
        
    def is_in_cache(self, key):
        return hasattr(self, key)


def read_pfd_file(pfdfn):
    """Return a Candidate object for the pfd given.
        Input:
            pfdfn: PRESTO *.pfd file name.
        Output:
            cand: Candidate object constructed from the given pfd
                file name.
    """
    pfd = prepfold.pfd(pfdfn)
    cand = Candidate(pfd.topo_p1, pfd.bary_p1, pfd.bestdm, \
                    psr_utils.ra_to_rad(pfd.rastr)*psr_utils.RADTODEG, \
                    psr_utils.dec_to_rad(pfd.decstr)*psr_utils.RADTODEG, \
                    pfdfn)
    return cand

def read_spd_file(spdfn):
    """Return a SPCandidate object for the spd given.
        Input:
            spdfn: *.spd file name.
        Output:
            cand: Candidate object constructed from the given pfd
                file name.
    """
    #spd = sp_utils.spd(spdfn)
    spd = read_spd.spd(spdfn)
    cand = SPCandidate(spd.best_dm, spd.ra_deg, spd.dec_deg, spdfn)
    return cand
