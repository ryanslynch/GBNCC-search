"""
Create Single Pulse Group object and methods to group Single Pulse Groups together.
Chen Karako
Oct. 28, 2010
"""
import numpy as np

# Defining some constants
# These thresholds are for determining whether two groups are close
# enough to combine
TIME_THRESH = 0.1 # s
DM_THRESH = 0.5 # pc cm-3

class SinglePulseGroup:
    """Define single pulse group
    """
    def __init__(self,dm,sigma,time,sample,downfact):
        """SinglePulseGroup constructor.
            Takes as input one single pulse (creates a group of one)
            inputs DM,Sigma,Time,Sample,Downfact.
        """
        self.min_dm = dm
        self.max_dm = dm
        self.max_sigma = sigma
        self.center_time = time
        if sample == 0:
            dt = 0 # this will ignore events with sample=0. better would be to use the inf files to determine the dt for these events
        else:
            dt = time/sample
#        dt = time/sample
        self.min_time = time-downfact/2.0*dt
        self.max_time = time+downfact/2.0*dt
        self.singlepulses = [(dm,sigma,time,sample,downfact)]
        self.numpulses = 1
        self.rank = 0
#        self.sp_params = np.array([(dm,sigma,time,sample,downfact)])
#        self.sp_params_trans = self.sp_params.transpose()

    def isclose(self,other):
        """Checks whether other is close enough to self to group the two
            SinglePulseGroup objects.
            Takes as input other, a SinglePulseGroup object.
        """
        timeisclose = (other.min_time <= self.center_time and \
                        other.max_time >= self.center_time) or \
                        (other.max_time <= self.center_time and \
                        other.max_time >= self.center_time-TIME_THRESH) or \
                        (other.min_time >= self.center_time and \
                        other.min_time <= self.center_time+TIME_THRESH)
        
        dmisclose = (other.min_dm >= self.min_dm and \
                    other.max_dm <= self.max_dm) or \
                    (other.min_dm <= self.min_dm and \
                    other.max_dm >= self.max_dm) or \
                    (other.min_dm <= self.max_dm+DM_THRESH and \
                    other.max_dm >= self.max_dm) or \
                    (other.max_dm >= self.min_dm-DM_THRESH and \
                    other.min_dm <= self.min_dm)
        
        return (timeisclose and dmisclose)

    def timeisclose(self,other,time_thresh=TIME_THRESH):
        """Checks whether the overlap in time of self and other is within
            time_thresh. Takes as input other, a SinglePulseGroup object,
            as well as the optional input time_thresh (in s).
        """
        timeisclose = (other.min_time <= self.center_time and \
                        other.max_time >= self.center_time) or \
                        (other.max_time <= self.center_time and \
                        other.max_time >= self.center_time-time_thresh) or \
                        (other.min_time >= self.center_time and \
                        other.min_time <= self.center_time+time_thresh)
        
        return timeisclose

    def dmisclose(self,other,dm_thresh=DM_THRESH):
        """Checks whether the DM of self and other is within dm_thresh of one
            another. Takes as input other, a SinglePulseGroup object, as well as the optional input dm_thresh (in pc cm-3).
        """
        dmisclose = (other.min_dm >= self.min_dm and \
                    other.max_dm <= self.max_dm) or \
                    (other.min_dm <= self.min_dm and \
                    other.max_dm >= self.max_dm) or \
                    (other.min_dm <= self.max_dm+dm_thresh and \
                    other.max_dm >= self.max_dm) or \
                    (other.max_dm >= self.min_dm-dm_thresh and \
                    other.min_dm <= self.min_dm)
        
        return dmisclose

    def combine(self,other):
        """combines self and other SinglePulseGroup objects.
            takes as input other, a SinglePulseGroup object.
            combines in place; nothing returned.
        """
        self.min_dm = min(self.min_dm, other.min_dm)
        self.max_dm = max(self.max_dm, other.max_dm)
        self.min_time = min(self.min_time, other.min_time)
        self.max_time = max(self.max_time, other.max_time)
        self.max_sigma = max(self.max_sigma, other.max_sigma)
        self.center_time = (self.min_time + self.max_time)/2.0
        self.numpulses = self.numpulses + other.numpulses
        self.singlepulses.extend(other.singlepulses)

    def __str__(self):
        s = ["Group of %d single pulses: " % len(self.singlepulses), \
             "\tMin DM (cm-3 pc): %f" % self.min_dm, \
             "\tMax DM (cm-3 pc): %f" % self.max_dm, \
             "\tCenter time (s):  %f" % self.center_time, \
             "\tDuration (s):     %f" % (self.max_time-self.min_time), \
             "\tMax sigma:        %f" % self.max_sigma, \
             "\tRank:             %f" % self.rank]
        return '\n'.join(s)
