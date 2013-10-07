"""
Single Pulse utilities: reading in single pulse files in current directory,
plotting DM vs. t for a list of groups (plotting the first five in colour,
the others in gray)
"""

import numpy as np
import glob
import pickle
import os.path
import matplotlib.pyplot as plt

def read_sp_files():
#def read_sp_files(filelist):
    """Read all *.singlepulse files in the current directory.
	Return 5 arrays (properties of all single pulses):
		DM, sigma, time, sample, downfact.
    """
    tmp_sp_params = np.empty((0,5))

    for file in glob.glob('*.singlepulse'):
#    for file in filelist:
	if os.path.getsize(file):
	    curr = np.atleast_2d(np.loadtxt(file))
	    tmp_sp_params = np.concatenate([tmp_sp_params, curr], axis=0)

# storing columns of dm, sigma, time, sample, downfact:   
    dm, sigma, time, sample, downfact = tmp_sp_params.transpose()

#    return tmp_sp_params.transpose()
    return tmp_sp_params


def save_groups(groups, filename="saved_groups.pkl"):
    """Given a list of groups save them to a file.
        
        Inputs:
            groups: a list of groups
            filename: filename to save as (Default: saved_groups.pkl)
    """
    file = open(filename, 'wb')
    pickle.dump(groups, file)
    file.close()


def load_groups(filename="saved_groups.pkl"):
    """Unpickle and return a list of groups from the given file.

        Input:
            filename: filename to save as (Default: saved_groups.pkl)
    """
    
    file = open(filename, 'r')
    groups = pickle.load(file)
    file.close()
    return groups


#def read_groups()
#    """Read in groups from text files containing the parameters of each sp in          a given group.
#    """
#    tmp_sp_params = np.empty((0,5))

#    for file in glob.glob('group*.txt'):
#	if os.path.getsize(file):
#	    curr = np.atleast_2d(np.loadtxt(file))
#	    tmp_sp_params = np.concatenate([tmp_sp_params, curr], axis=0)

## storing columns of dm, sigma, time, sample, downfact:   
#    dm, sigma, time, sample, downfact = tmp_sp_params.transpose()

##    return tmp_sp_params.transpose()
#    return tmp_sp_params

def cmp_numpulses(grp1,grp2):
    """Compare the number of pulses in two Single Pulse Group objects.
        Return the usual -1,0,1 depending on relative sizes of two arguments.
    """
    return cmp(grp1.numpulses,grp2.numpulses)


def plot_sp_largest(groups, nlargest=0):
    """Take in list of Single Pulse Groups, sort them by decreasing size,
        ie. numpulses, and plot the DM vs. t for all, distinguishing the
        nlargest largest groups by plotting them in different colours.
    """
    groups.sort(cmp=cmp_numpulses, reverse=True)
    
    dm = np.empty(0)
    time = np.empty(0)
    sigma = np.empty(0)

    for grp in groups[nlargest:]:
        sp_array = np.array(grp.singlepulses).transpose()
        dm = np.concatenate([dm, sp_array[0]], axis=0)
        time = np.concatenate([time, sp_array[2]], axis=0)
        sigma = np.concatenate([sigma,sp_array[1]], axis=0)
    
    plotsize = np.clip(3*sigma-14,0,50)
    plt.scatter(time, dm, c='k', marker='o', s=plotsize, edgecolor='none')    
    plt.xlabel('Time (s)')
    plt.ylabel('DM (pc cm-3)')

    from matplotlib import cm
    
    for i in range(0,nlargest):
        plt.scatter(np.array(groups[i].singlepulses).transpose()[2], \
                    np.array(groups[i].singlepulses).transpose()[0], marker='o', \
                    edgecolor='none',color=cm.jet(1.*i/nlargest), \
                    s=np.clip(3*np.array(groups[i].singlepulses).transpose()[1]\
                    -14,0,50))
    
    plt.ylim((0,100))
    plt.xlim((0,155))

    from time import localtime
    datetime_stamp = '%4d-%02d-%02dT%02d-%02d-%02d' % localtime()[:6]
    title = 'grouped_sps'
    plotted = 'largest'
    ext = 'eps'
#    plt.savefig('%s-%s-%s.%s' % (title, plotted, datetime_stamp, ext))
#    plt.show()

    
def cmp_maxsigma(grp1,grp2):
    """Compare the maximum sigma in two Single Pulse Group objects.
        Return the usual -1,0,1 depending on relative sizes of two arguments.
    """
    return cmp(grp1.max_sigma,grp2.max_sigma)


def plot_sp_mostsig(groups, nmostsig=0):
    """Take in list of Single Pulse Groups, sort them by decreasing sigma_max,
        and plot the DM vs. t for all, distinguishing the
        nmostsig most significant groups by plotting them in different colours.
    """
    groups.sort(cmp=cmp_maxsigma, reverse=True)
    
    dm = np.empty(0)
    time = np.empty(0)
    sigma = np.empty(0)

    for grp in groups[nmostsig:]:
        sp_array = np.array(grp.singlepulses).transpose()
        dm = np.concatenate([dm, sp_array[0]], axis=0)
        time = np.concatenate([time, sp_array[2]], axis=0)
        sigma = np.concatenate([sigma,sp_array[1]], axis=0)
    
    plotsize = np.clip(3*sigma-14,0,50)
    plt.scatter(time, dm, c='k', marker='o', s=plotsize)    
    plt.xlabel('Time (s)')
    plt.ylabel('DM (pc cm-3)')

    from matplotlib import cm
    
    for i in range(0,nmostsig):
        plt.scatter(np.array(groups[i].singlepulses).transpose()[2], \
                    np.array(groups[i].singlepulses).transpose()[0], marker='o', \
                    edgecolor='none',color=cm.jet(1.*i/nmostsig), \
                    s=np.clip(3*np.array(groups[i].singlepulses).transpose()[1]\
                    -14,0,50))

    plt.ylim((0,100))
    plt.xlim((0,155))

    from time import localtime
    datetime_stamp = '%4d-%02d-%02dT%02d-%02d-%02d' % localtime()[:6]
    title = 'grouped_sps'
    plotted = 'mostsig'
    ext = 'eps'
    plt.savefig('%s-%s-%s.%s' % (title, plotted, datetime_stamp, ext))
    plt.show()

def plot_dmvsig_largest(groups, nlargest=0):
    """Take in list of Single Pulse Groups, sort them by decreasing size,
        ie. numpulses, and plot the DM vs. sigma for all, distinguishing the
        nlargest largest groups by plotting them in different colours.
    """
    groups.sort(cmp=cmp_numpulses, reverse=True)
    
    dm = np.empty(0)
    time = np.empty(0)
    sigma = np.empty(0)

    for grp in groups[nlargest:]:
        sp_array = np.array(grp.singlepulses).transpose()
        dm = np.concatenate([dm, sp_array[0]], axis=0)
        time = np.concatenate([time, sp_array[2]], axis=0)
        sigma = np.concatenate([sigma,sp_array[1]], axis=0)
    
#    plotsize = np.clip(3*sigma-14,0,50)
#    plt.scatter(dm, sigma, c='k', marker='o', s=plotsize, edgecolor='none')    
    plt.ylabel('Sigma')
    plt.xlabel('DM (pc cm-3)')

    from matplotlib import cm
    
    for i in range(0,nlargest):
        plt.scatter(np.array(groups[i].singlepulses).transpose()[0], \
                    np.array(groups[i].singlepulses).transpose()[1], marker='o', \
                    edgecolor='none',color=cm.jet(1.*i/nlargest),s=4)
#                    s=np.clip(3*np.array(groups[i].singlepulses).transpose()[1]\
#                    -14,0,50))
    
    plt.xlim((0,100))
    plt.ylim((4.8,8))

    from time import localtime
    datetime_stamp = '%4d-%02d-%02dT%02d-%02d-%02d' % localtime()[:6]
    title = 'dmvsig_grouped_sps'
    plotted = 'largestonly'
    ext = 'eps'
    plt.savefig('%s-%s-%s.%s' % (title, plotted, datetime_stamp, ext))
    plt.show()

def plot_single_sp_largest(groups, nlargest=1):
    """Take in list of Single Pulse Groups, sort them by decreasing size,
        ie. numpulses, and plot the DM vs. t for the nlargest largest
        groups only, in different colours.
    """
    groups.sort(cmp=cmp_numpulses, reverse=True)
    
    dm = np.empty(0)
    time = np.empty(0)
    sigma = np.empty(0)

#    for grp in groups[nlargest:]:
#        sp_array = np.array(grp.singlepulses).transpose()
#        dm = np.concatenate([dm, sp_array[0]], axis=0)
#        time = np.concatenate([time, sp_array[2]], axis=0)
#        sigma = np.concatenate([sigma,sp_array[1]], axis=0)
    
#    plotsize = np.clip(3*sigma-14,0,50)
#    plt.scatter(time, dm, c='k', marker='o', s=plotsize, edgecolor='none')    
    plt.xlabel('Time (s)')
    plt.ylabel('DM (pc cm-3)')

    from matplotlib import cm
 
    for i in range(0,nlargest):
        sp_array = np.array(groups[i].singlepulses).transpose()
        dm = np.array(groups[i].singlepulses).transpose()[0]
        time = np.array(groups[i].singlepulses).transpose()[2]
        sigma = np.array(groups[i].singlepulses).transpose()[1]
        plotsize = np.clip(3*np.array(groups[i].singlepulses).transpose()[1]-14,0,50)

        plt.scatter(time, dm, marker='o', edgecolor='none', \
                    color=cm.jet(1.*i/nlargest), s=plotsize)

#    for i in range(0,nlargest):
#        plt.scatter(np.array(groups[i].singlepulses).transpose()[2], \
#                    np.array(groups[i].singlepulses).transpose()[0], marker='o', \
#                    edgecolor='none',color=cm.jet(1.*i/nlargest), \
#                    s=np.clip(3*np.array(groups[i].singlepulses).transpose()[1]\
#                    -14,0,50))
    
    plt.ylim((0,100))
    plt.xlim((0,155))

#    from time import localtime
#    datetime_stamp = '%4d-%02d-%02dT%02d-%02d-%02d' % localtime()[:6]
#    title = 'grouped_sps'
#    plotted = 'largestonly'
#    ext = 'eps'
#    plt.savefig('%s-%s-%s.%s' % (title, plotted, datetime_stamp, ext))
#    plt.show()

    
def plot_sp_rated(groups, color='k', ylow=0, yhigh=100, xlow=0, xhigh=120):
    """Take in list of Single Pulse Groups
        and plot the DM vs. t for all, with the plotted colour automatically black
        or another, when specified.
    """
    groups.sort(cmp=cmp_numpulses, reverse=True)
    
    dm = np.empty(0)
    time = np.empty(0)
    sigma = np.empty(0)

    for grp in groups:
        sp_array = np.array(grp.singlepulses).transpose()
        dm = np.concatenate([dm, sp_array[0]], axis=0)
        time = np.concatenate([time, sp_array[2]], axis=0)
        sigma = np.concatenate([sigma,sp_array[1]], axis=0)
    
    plotsize = np.clip(3*sigma-14,0,50)
    plt.scatter(time, dm, c=color, marker='o', s=plotsize, edgecolor='none')    
    plt.xlabel('Time (s)')
    plt.ylabel('DM (pc cm-3)')

    from matplotlib import cm
    
    plt.ylim((ylow, yhigh))
    plt.xlim((xlow, xhigh))

#    from time import localtime
#    datetime_stamp = '%4d-%02d-%02dT%02d-%02d-%02d' % localtime()[:6]
#    title = 'grouped_sps'
#    plotted = 'rated'
#    ext = 'eps'
#    plt.savefig('%s-%s-%s.%s' % (title, plotted, datetime_stamp, ext))
#    plt.show()


def plot_sp_rated_onesize(groups, color='k', ylow=0, yhigh=100, xlow=0, xhigh=140, plotsize=10):
    """Take in list of Single Pulse Groups
        and plot the DM vs. t for all, with the plotted colour automatically black or another, when specified.
    """
    groups.sort(cmp=cmp_numpulses, reverse=True)
    
    dm = np.empty(0)
    time = np.empty(0)
    sigma = np.empty(0)

    for grp in groups:
        sp_array = np.array(grp.singlepulses).transpose()
        dm = np.concatenate([dm, sp_array[0]], axis=0)
        time = np.concatenate([time, sp_array[2]], axis=0)
#        sigma = np.concatenate([sigma,sp_array[1]], axis=0)
    
#    plotsize = np.clip(3*sigma-14,0,50)
    plt.scatter(time, dm, c=color, marker='o', s=plotsize, edgecolor='none')    
    plt.xlabel('Time (s)')
    plt.ylabel(r'$DM (pc cm^{-3})$')

    from matplotlib import cm
    
    plt.ylim((ylow, yhigh))
    plt.xlim((xlow, xhigh))

#    from time import localtime
#    datetime_stamp = '%4d-%02d-%02dT%02d-%02d-%02d' % localtime()[:6]
#    title = 'grouped_sps'
#    plotted = 'rated'
#    ext = 'eps'
#    plt.savefig('%s-%s-%s.%s' % (title, plotted, datetime_stamp, ext))
#    plt.show()


