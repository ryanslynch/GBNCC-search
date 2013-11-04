"""
Single Pulse utilities: reading in single pulse files in current directory,
plotting DM vs. t for a list of groups (plotting the first five in colour,
the others in gray)
"""

import numpy as np
import glob
import pickle
import os.path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import infodata
import group_sp

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

def read_groups(grpfn='groups.txt'):
    grps = []
    rank0groups = []
    rfigroups = []
    rank3groups = []
    rank4groups = []
    rank5groups = []
    rank6groups = []
    size = os.path.getsize(grpfn)
    grpfile = open(grpfn, 'r')

    while True:
        if grpfile.tell() == size: # check if we've reached the end of the file
            break
        line = grpfile.readline()
        if line.startswith("Group of"):
            spgrp = []
            Nsp = int(line.split()[2])
            for j in range(6): # skip ahead to 6 lines down
                line = grpfile.readline()
            rank = float(line.split()[1])
#            for j in range(3):
#                grpfile.readline() # skip 3 rows down: works for GBNCC output
            for j in range(2):
                grpfile.readline() # skip 2 rows down: for PALFA output
            spdata = np.fromfile(grpfile, count=Nsp*5, sep=' ')
            spdata.shape = (Nsp, 5)
            spgrp = group_sp.SinglePulseGroup(*spdata[0])
            for sp in spdata[1:]:
                spgrp.combine(group_sp.SinglePulseGroup(*sp))
            # now assign this group the saved rank
            spgrp.rank = rank
            if spgrp.rank == 2:
                rfigroups.append(spgrp) # add to rfigroups list but not rest of groups
            else:
                grps.append(spgrp)
                if spgrp.rank == 0:
                    rank0groups.append(spgrp)
                if spgrp.rank == 3:
                    rank3groups.append(spgrp)
                if spgrp.rank == 4:
                    rank4groups.append(spgrp)
                if spgrp.rank == 5:
                    rank5groups.append(spgrp)
                if spgrp.rank == 6:
                    rank6groups.append(spgrp)
    grpfile.close()

#    outfile = open("groups-new.txt", 'w')
#    for grp in grps:
#        outfile.write(str(grp) + '\n')
#        outfile.write('\n')
#        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
#        for sp in grp.singlepulses:
#            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
#        outfile.write('\n')
#    outfile.close()

    return {'groups': grps, 'rfigroups': rfigroups, 'rank0groups': rank0groups, 'rank3groups': rank3groups, 'rank4groups': rank4groups, 'rank5groups': rank5groups, 'rank6groups': rank6groups}

def read_groups_ignoreobsend(grpfn='groups.txt'):
    """Read in groups from groups.txt file. Save groups as usual, but add any group with center time > (obs_end_time - 10s) to rfigrps. This is a patch to handle bad signals in ends of observations with zero-dm output, for beams that have already been processed.
    """
    grps = []
    rank0groups = []
    rfigroups = []
    rank3groups = []
    rank4groups = []
    rank5groups = []
    rank6groups = []
    size = os.path.getsize(grpfn)
    grpfile = open(grpfn, 'r')
    obsend = get_obs_info()['T'] 
    nearend = obsend - 10.0 #obs end time - 10 seconds 

    while True:
        if grpfile.tell() == size: # check if we've reached the end of the file
            break
        line = grpfile.readline()
        if line.startswith("Group of"):
            spgrp = []
            Nsp = int(line.split()[2])
            for j in range(6): # skip ahead to 6 lines down
                line = grpfile.readline()
            rank = float(line.split()[1])
#            for j in range(3):
#                grpfile.readline() # skip 3 rows down: works for GBNCC output
            for j in range(2):
                grpfile.readline() # skip 2 rows down: for PALFA output
            spdata = np.fromfile(grpfile, count=Nsp*5, sep=' ')
            spdata.shape = (Nsp, 5)
            spgrp = group_sp.SinglePulseGroup(*spdata[0])
            for sp in spdata[1:]:
                spgrp.combine(group_sp.SinglePulseGroup(*sp))
            # now assign this group the saved rank
            spgrp.rank = rank
            if spgrp.rank == 2 or spgrp.center_time >= nearend:
                spgrp.rank = 2 #assign RFI rank also to groups near obs end
                rfigroups.append(spgrp) # add to rfigroups list but not rest of groups
            else:
                grps.append(spgrp)
                if spgrp.rank == 0:
                    rank0groups.append(spgrp)
                if spgrp.rank == 3:
                    rank3groups.append(spgrp)
                if spgrp.rank == 4:
                    rank4groups.append(spgrp)
                if spgrp.rank == 5:
                    rank5groups.append(spgrp)
                if spgrp.rank == 6:
                    rank6groups.append(spgrp)
    grpfile.close()

    outfile = open("groups-new.txt", 'w')
    
    # beam summary
    outfile.write("Number of rank 0 groups: %d \n" %len(rank0groups))
    outfile.write("Number of rank 2 groups: %d \n" %len(rfigroups))
    outfile.write("Number of rank 3 groups: %d \n" %len(rank3groups))
    outfile.write("Number of rank 4 groups: %d \n" %len(rank4groups))
    outfile.write("Number of rank 5 groups: %d \n" %len(rank5groups))
    outfile.write("Number of rank 6 groups: %d \n" %len(rank6groups))
    outfile.write("\n")

    # print list of events in each group
    for grp in grps:
        outfile.write(str(grp) + '\n')
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')
    for grp in rfigroups:
        outfile.write(str(grp) + '\n')
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')
    outfile.close()

    return {'groups': grps, 'rfigroups': rfigroups, 'rank0groups': rank0groups, 'rank3groups': rank3groups, 'rank4groups': rank4groups, 'rank5groups': rank5groups, 'rank6groups': rank6groups}

def recover_plots(grpfn='groups.txt'):
    """Function to use when the groups.txt file is available, but the plots need to be remade. Loads groups using read_groups and the default groups.txt file, and make plots as in the original analyse_sp code. Input is groups.txt file (by default), and outputs several DM vs. t single pulse plots. Note an inf file needs to be available in the working dir."""

    groups = read_groups(grpfn)['groups']
    rfigroups = read_groups(grpfn)['rfigroups']
    rank3groups = read_groups(grpfn)['rank3groups']
    rank4groups = read_groups(grpfn)['rank4groups']
    rank5groups = read_groups(grpfn)['rank5groups']
    rank6groups = read_groups(grpfn)['rank6groups']

# Plot groups, using the same plotting scheme as in analyse_sp
#    DM 0-30
    plt.figure(figsize=(10,8))
    if len(rfigroups) != 0:
        plot_sp_rated(rfigroups, 'r', 0, 30)
    if len(groups) != 0:
        plot_sp_rated(groups, 'k', 0, 30)
    if len(rank3groups) != 0:
        plot_sp_rated(rank3groups, 'g', 0, 30)
    if len(rank4groups) != 0:
        plot_sp_rated(rank4groups, 'b', 0, 30)
    if len(rank5groups) != 0:
        plot_sp_rated(rank5groups, 'm', 0, 30)
    if len(rank6groups) != 0:
        plot_sp_rated(rank6groups, 'c', 0, 30)
    plt.savefig('grouped_sps_DMs0-30.png', dpi = 300)

#   DM 20-110
    plt.figure(figsize=(10,8))
    if len(rfigroups) != 0:
        plot_sp_rated(rfigroups, 'r', 20, 110)
    if len(groups) != 0:
        plot_sp_rated(groups, 'k', 20, 110)
    if len(rank3groups) != 0:
        plot_sp_rated(rank3groups, 'g', 20, 110)
    if len(rank4groups) != 0:
        plot_sp_rated(rank4groups, 'b', 20, 110)
    if len(rank5groups) != 0:
        plot_sp_rated(rank5groups, 'm', 20, 110)
    if len(rank6groups) != 0:
        plot_sp_rated(rank6groups, 'c', 20, 110)
    plt.savefig('grouped_sps_DMs20-110.png', dpi = 300)

#    DM 100-310
    plt.figure(figsize=(10,8))
    if len(rfigroups) != 0:
        plot_sp_rated(rfigroups, 'r', 100, 310)
    if len(groups) != 0:
        plot_sp_rated(groups, 'k', 100, 310)
    if len(rank3groups) != 0:
        plot_sp_rated(rank3groups, 'g', 100, 310)
    if len(rank4groups) != 0:
        plot_sp_rated(rank4groups, 'b', 100, 310)
    if len(rank5groups) != 0:
        plot_sp_rated(rank5groups, 'm', 100, 310)
    if len(rank6groups) != 0:
        plot_sp_rated(rank6groups, 'c', 100, 310)
    plt.savefig('grouped_sps_DMs100-310.png', dpi = 300)

#    DM 300-1100
    plt.figure(figsize=(10,8))
    if len(rfigroups) != 0:
        plot_sp_rated(rfigroups, 'r', 300, 1100)
    if len(groups) != 0:
        plot_sp_rated(groups, 'k', 300, 1100)
    if len(rank3groups) != 0:
        plot_sp_rated(rank3groups, 'g', 300, 1100)
    if len(rank4groups) != 0:
        plot_sp_rated(rank4groups, 'b', 300, 1100)
    if len(rank5groups) != 0:
        plot_sp_rated(rank5groups, 'm', 300, 1100)
    if len(rank6groups) != 0:
        plot_sp_rated(rank6groups, 'c', 300, 1100)
    plt.savefig('grouped_sps_DMs300-1100.png', dpi = 300)

#    DM 1000-5000
    plt.figure(figsize=(10,8))
    if len(rfigroups) != 0:
        plot_sp_rated(rfigroups, 'r', 1000, 5000)
    if len(groups) != 0:
        plot_sp_rated(groups, 'k', 1000, 5000)
    if len(rank3groups) != 0:
        plot_sp_rated(rank3groups, 'g', 1000, 5000)
    if len(rank4groups) != 0:
        plot_sp_rated(rank4groups, 'b', 1000, 5000)
    if len(rank5groups) != 0:
        plot_sp_rated(rank5groups, 'm', 1000, 5000)
    if len(rank6groups) != 0:
        plot_sp_rated(rank6groups, 'c', 1000, 5000)
    plt.savefig('grouped_sps_DMs1000-5000.png', dpi = 300)

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

def get_obs_info():
    """Read in an .inf file to extract observation information.
        Return observation RA, Dec, duration, and source name.
    """
    inffiles = glob.glob('*.inf')
    if len(inffiles) == 0: # no inf files exist
        print "No inf files available!"
        return None
    else:
        inffile = inffiles[0] # take the first inf file in current dir
        inf = infodata.infodata(inffile)
        T = inf.dt * inf.N # total observation time (s)
        RA = inf.RA
        dec = inf.DEC
        src = inf.object
        return {'T': T, 'RA': RA, 'dec': dec, 'src': src}

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

    
def plot_sp_rated(groups, color='k', ylow=0, yhigh=100, xlow=0, xhigh=140):
    """Take in list of Single Pulse Groups
        and plot the DM vs. t for all, with the plotted colour automatically black
        or another, when specified.
    """
#    groups.sort(cmp=cmp_numpulses, reverse=True)
    
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
    plt.ylabel('DM (pc cm$^{-3}$)')
    if get_obs_info() is not None: # if inf files exist, can get obs info
        plt.title('Single Pulse Results for %s\nRA: %s Dec: %s' % (get_obs_info()['src'], get_obs_info()['RA'], get_obs_info()['dec']))
        xhigh = get_obs_info()['T'] # set xhigh to observation duration
    else:
        xhigh = 120
    xlow = 0
    plt.xlim((xlow, xhigh)) 
    plt.ylim((ylow, yhigh))

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


