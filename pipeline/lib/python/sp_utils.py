"""
Single Pulse utilities: reading in single pulse files in current directory,
plotting DM vs. t for a list of groups (plotting the first five in colour,
the others in gray)
"""

import numpy as np
import glob
import os.path
#import matplotlib and set for non-interactive plots
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

    dm2 = np.empty(0)
    time2 = np.empty(0)
    sigma2 = np.empty(0)
   
# roundabout way to plot 8 (inflexible) largest groups in different colours  
    plt.scatter(np.array(groups[0].singlepulses).transpose()[2], \
                np.array(groups[0].singlepulses).transpose()[0], marker='o', \
                edgecolor='none', c='r',\
                s=np.clip(3*np.array(groups[0].singlepulses).transpose()[1]\
                -14,0,50))

    plt.scatter(np.array(groups[1].singlepulses).transpose()[2], \
                np.array(groups[1].singlepulses).transpose()[0], marker='o', \
                edgecolor='none', c='b',\
                s=np.clip(3*np.array(groups[1].singlepulses).transpose()[1]\
                -14,0,50))

    plt.scatter(np.array(groups[2].singlepulses).transpose()[2], \
                np.array(groups[2].singlepulses).transpose()[0], marker='o', \
                edgecolor='none', c='g',\
                s=np.clip(3*np.array(groups[2].singlepulses).transpose()[1]\
                -14,0,50))
    
    plt.scatter(np.array(groups[3].singlepulses).transpose()[2], \
                np.array(groups[3].singlepulses).transpose()[0], marker='o', \
                edgecolor='none', c='m',\
                s=np.clip(3*np.array(groups[3].singlepulses).transpose()[1]\
                -14,0,50))
    
    plt.scatter(np.array(groups[4].singlepulses).transpose()[2], \
                np.array(groups[4].singlepulses).transpose()[0], marker='o', \
                edgecolor='none', c='c',\
                s=np.clip(3*np.array(groups[4].singlepulses).transpose()[1]\
                -14,0,50))

    plt.scatter(np.array(groups[5].singlepulses).transpose()[2], \
                np.array(groups[5].singlepulses).transpose()[0], marker='o', \
                edgecolor='none', c='y',\
                s=np.clip(3*np.array(groups[5].singlepulses).transpose()[1]\
                -14,0,50))
    
    plt.scatter(np.array(groups[6].singlepulses).transpose()[2], \
                np.array(groups[6].singlepulses).transpose()[0], marker='o', \
                edgecolor='none', c='#00cc00',\
                s=np.clip(3*np.array(groups[6].singlepulses).transpose()[1]\
                -14,0,50))

    plt.scatter(np.array(groups[7].singlepulses).transpose()[2], \
                np.array(groups[7].singlepulses).transpose()[0], marker='o', \
                edgecolor='none', c='#ff6600',\
                s=np.clip(3*np.array(groups[7].singlepulses).transpose()[1]\
                -14,0,50))

    plt.ylim((0,100))
    plt.xlim((0,140))
#    plt.show()

    from time import localtime
    datetime_stamp = '%4d-%02d-%02dT%02d-%02d-%02d' % localtime()[:6]
    title = 'grouped_sps'
    ext = 'eps'
    plt.savefig('%s-%s.%s' % (title, datetime_stamp,ext))
    plt.show()

"""
    for grp in groups[0:nlargest]:
        sp_array2 = np.array(grp.singlepulses).transpose()
        dm2 = np.concatenate([dm2, sp_array2[0]], axis=0)
        time2 = np.concatenate([time2, sp_array2[2]], axis=0)
        sigma2 = np.concatenate([sigma2, sp_array2[1]], axis=0)
        plotsize = np.clip(3*sp_array2[1]-14,0,50)
        plt.scatter(sp_array2[2], sp_array2[0], marker='o', edgecolor='none', \
        s=plotsize)
        plt.ylim((0,100))
        plt.xlim((0,140))

#save obtained plot with time generated in filename
    from time import localtime
    datetime_stamp = '%4d-%02d-%02dT%02d-%02d-%02d' % localtime()[:6]
    title = 'grouped_sps'
    ext = 'eps'
    plt.savefig('%s-%s.%s' % (title, datetime_stamp,ext))
    plt.show()

    for grp in groups[nlargest:nlargest+4]:
        sp_array3 = np.array(grp.singlepulses).transpose()
        plotsize = np.clip(3*sp_array3[1]-14,0,50)
        plt.scatter(sp_array3[2], sp_array3[0], marker='o', edgecolor='none', \
        s=plotsize)
        plt.ylim((0,100))
    
    plt.show()
"""
"""
#        plt.scatter(sp_array2[2], sp_array2[0], marker='o', edgecolor='none')
    plotsize = np.clip(3*sigma2-14,0,50)
    plt.scatter(time2, dm2, c='r', marker='o', edgecolor='none',s=plotsize)
    plt.ylim((0,100))
    plt.show()
""" 
