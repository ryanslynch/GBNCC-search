import os, glob, config, psr_utils
from subprocess import Popen, PIPE
import numpy as np
import matplotlib.pyplot as plt

def getqueue(machine,queue):
    if machine == "guillimin" or machine == "nimrod":
        alljobs = queue.getjobs()
        if alljobs is not None:
            myjobs    = [job for job in alljobs.itervalues() \
                         if job.has_key("euser") and \
                         job["euser"][0] == config.user]
        else: myjobs = []
        nqueued = len([job for job in myjobs if job["job_state"][0] == "Q" \
                       or job["job_state"] == "B"])
    
    return nqueued

def subjob(machine, subfilenm, options=""):
    if machine == "guillimin" or machine == "nimrod":
        cmd = "qsub " + options + " " + subfilenm
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        jobid,err = process.communicate()
            
    if len(err) != 0: return None, err.strip()
    else: return jobid, None


def results_status(outdir, basenm):
    npfdplots   = len(glob.glob(os.path.join(outdir, "%s*.pfd.png"%basenm)))
    nspplots    = len(glob.glob(os.path.join(outdir, "%s*singlepulse.png"%basenm)))
    nratings    = len(glob.glob(os.path.join(outdir, "%s*.pfd.ratings"%basenm)))
    nrfifinds   = len(glob.glob(os.path.join(outdir, "%s*rfifind*"%basenm)))
    ntgzs       = len(glob.glob(os.path.join(outdir, "%s*.tgz"%basenm)))
    nreport     = len(glob.glob(os.path.join(outdir, "%s*.report"%basenm)))
    ndiagnostic = len(glob.glob(os.path.join(outdir, "%s*.diagnostics"%basenm)))
    naccels     = len(glob.glob(os.path.join(outdir, "%s*accelcands*"%basenm)))
    ngroup      = len(glob.glob(os.path.join(outdir, "groups.txt")))
    #ngroup = 1
    ngroupplots = len(glob.glob(os.path.join(outdir, "grouped*png")))
    #ngroupplots = 1

    if (naccels != 2) or (ntgzs !=9) or (nrfifinds < 8) or (nreport != 1) or \
       (npfdplots == 0) or (nspplots < 6):
        return "f"
    elif (npfdplots != nratings) or (ndiagnostic != 1) or (ngroup != 1) or \
         (ngroupplots == 0):
        return "w"
    else:
        return "s"

def get_size(path):
    size = 0
    for dirpath,dirnms,filenms in os.walk(path):
        for f in filenms:
            size += os.path.getsize(os.path.join(dirpath,f))

    return size


def print_sp_raters_list(verbosity=0):
    """Print the list of imported single pulse raters to stdout.
        
        Input:
            verbosity: If True, print description of raters.
                (Default: Don't be verbose.)
        Outputs:
            None
    """
    import textwrap
    import sp_raters
    print "Number of single pulse raters registered: %d" % len(sp_raters.registered_raters)
    for rater_name in sp_raters.registered_raters:
        rater_module = getattr(sp_raters, rater_name)
        rater = rater_module.Rater()
        print "'%s': %s (v%d)" % (rater_name, rater.long_name, rater.version)
        if verbosity:
            print ""
            for line in rater.description.split('\n'):
                print textwrap.fill(line, width=70)
            print "-"*25

def get_scaled_profile(profile, varprof):
    scaled = profile.copy()
    scaled /= np.sqrt(varprof)
    scaled -= scaled.mean()
    return scaled

def multigaussfit_from_paramlist(params):
    comps = []
    for ii in range(1, len(params), 3):
        amp = params[ii]
        std = np.abs(params[ii+1])
        phs = params[ii+2]
        comps.append(MultiGaussComponent(amp, std, phs))
    fit = MultiGaussFit(offset=params[0], components=comps)
    return fit

class MultiGaussComponent(object):
    def __init__(self, amp, fwhm, phs):
        """Constructor for MultiGaussComponent, an object to represent
            a single gaussian component of a multiple-gaussian fit to
            a profile.

            Inputs:
                amp: The amplitude of the gaussian component.
                fwhm: The full-width at half-maximum of the gaussian component.
                phs: The phase of the gaussian component.

            Output:
                component: The MultiGaussComponent object.
        """
        self.amp = amp
        self.fwhm = fwhm
        self.phs = phs

    def __str__(self):
        s = "Amplitude: %g, FWHM: %g, Phase: %g" % \
                    (self.amp, self.fwhm, self.phs)
        return s

    def make_gaussian(self, nbins):
        """Return an aray of length 'nbins' containing the gaussian component.

           Inputs:
               nbins: The number of bins in the profile

           Output:
               gaussian: Array of data

        """
        # Create an array for the Gaussian profile
        gaussian = self.amp*self.fwhm/2*np.sqrt(np.pi/np.log(2)) * \
                        psr_utils.gaussian_profile(nbins,self.phs,self.fwhm)
        return gaussian

    def get_onpulse_region(self, nbins):
        """Return a tuple of phases that represent the on-pulse window.

            Inputs:
                nbins: Number of phase bins.

            Output:
                onpulse: A tuple of phases, between which are the
                    on-pulse region.
        """
        # Determine fudge factor depending on width
        if self.fwhm < 0.1:
            fudge_factor = 4.0
        elif self.fwhm < 0.2:
            fudge_factor = 2.0
        elif self.fwhm < 0.4:
            fudge_factor = 1.5
        else:
            fudge_factor = 1.5

        if self.fwhm*fudge_factor > 1.0:
            raise RatingError("Fudge factored FWHM is larger than 1.0 in phase")

        start_phase = self.phs - (self.fwhm*fudge_factor)/2.0
        end_phase = self.phs + (self.fwhm*fudge_factor)/2.0

        start_phase %= 1
        end_phase %= 1

        start_bin = int(start_phase*nbins+0.5) # Round to nearest integer
        end_bin = int(end_phase*nbins+0.5) # Round to nearest integer
        onpulse_length = (end_bin - start_bin) % nbins
        onpulse_indices = np.arange(start_bin, start_bin+onpulse_length) % nbins
        onpulse_region = np.zeros(nbins, dtype=bool)
        onpulse_region[onpulse_indices] = True
        return onpulse_region



class MultiGaussFit(object):
    def __init__(self, offset=0.0, components=[]):
        """Constructor for MultiGaussFit, a multiple-gaussian fit to
            a profile.

            Inputs:
                offset: The DC offset of the fit. (Default: 0.0)
                components: A list of MultiGaussComponents making up the
                    fit to the profile. (Default: No components)

            Output:
                fit: The MultiGaussFit object.
        """
        self.offset = offset
        self.components = components

    def __str__(self):
        lines = ["Multi-Gaussian fit with %d components" % \
                    len(self.components)]
        for ii, comp in enumerate(self.components):
            lines.append("    Component %d: %s" % (ii+1, str(comp)))
        return '\n'.join(lines)

    def add_component(self, comp):
        """Add a component to the MultiGaussianFit.

            Input:
                comp: A MultiGaussComponent to add.

            Outputs:
                None
        """
        self.components.append(comp)

    def make_gaussians(self, nbins):
        """Return an array of length 'nbins' containing the gaussian fit.

           Inputs:
               nbins: The number of bins in the profile

           Output:
               gaussian: Array of data

        """
        # Determine the number of Gaussian profiles to make
        ngaussians = len(self.components)

        # Create an array for the Gaussian profile
        gaussians = np.zeros(nbins) + self.offset

        # Add each individual Gaussian to the full profile
        for comp in self.components:
            #print "DEBUG: comp.amp, comp.std, comp.phs", comp.amp, comp.std, comp.phs
            gaussians += comp.make_gaussian(nbins)
        return gaussians

    def get_resids(self, data):
        model = self.make_gaussians(len(data))
        resids = data - model
        return resids

    def get_chisqr(self, data):
        resids = self.get_resids(data)
        return np.sum(resids**2)

    def get_dof(self, nbins):
        return nbins - self.get_num_params()

    def get_num_params(self):
        return 1 + 3*len(self.components)

    def plot_comparison(self, data, individual=False):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        phases = np.linspace(0.0, 1.0, len(data), endpoint=False)
        phases_10x = np.linspace(0.0, 1.0, len(data)*10, endpoint=False)

        ax = plt.axes([0.1, 0.1, 0.85, 0.65])
        plt.plot(phases, data, c='k', label="Profile")

        plt.plot(phases_10x, self.make_gaussians(len(data)*10), c='r', label="Fit")
        if individual:
            for comp in self.components:
                plt.plot(phases_10x, self.offset+comp.make_gaussian(len(data)*10), ls='--')
        plt.xlabel("Phase")
        plt.ylabel("Intensity (arbitrary units)")
        plt.legend(loc='best')

        ax = plt.axes([0.1, 0.75, 0.85, 0.15], sharex=ax)
        plt.plot(phases, self.get_resids(data), c='k')
        plt.axhline(0.0, color='k', ls='--')
        plt.ylabel("Residuals")
        plt.setp(ax.xaxis.get_ticklabels(), visible=False)
        plt.show()

    def get_onpulse_region(self, nbins):
        """Return a tuple of phases that represent the on-pulse window.

            Inputs:
                nbins: Number of phase bins.

            Output:
                onpulse: A tuple of phases, between which are the
                    on-pulse region.
        """
        if not self.components:
            raise RatingError("Multi-Gauss fit has no components " \
                                    "(i.e. no on-pulse region)")
        onpulse_region = np.zeros(nbins, dtype=bool)
        for comp in self.components:
            onpulse_region |= comp.get_onpulse_region(nbins)
        return onpulse_region



class RatingError(Exception):
    pass

class RatingWarning(Warning):
    pass

class RaterLoadWarning(RatingWarning):
    pass

