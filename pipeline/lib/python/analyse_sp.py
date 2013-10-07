#!/usr/bin/env python

"""
Analyse all single pulse files in a given beam by grouping them, ranking the groups, and plotting them in colours corresponding to their rating.
Modifications to analyse_sp3.py script which aim to improve the grouping and rating algorithms.
Modifications to analyse_sp_v2.py which aim to better classify RFI at high DMs and avoid it being rated highly.
Chen Karako
Nov. 4, 2011

Making changes to analyse_sp_v3.py in order to fix bugs and eventually integrate code into GBNCC pipelines.
Oct 2, 2012
"""
import numpy as np
import sp_utils
import sp_utils2_gbncc as sp_utils2 #note sp_utils2_palfa differs from _gbncc in that does not have predefined xlimits (ie. decided automatically, not nec 0 to 120s)
import group_sp
import sys

CLOSE_DM = 2 # pc cm-3 #consider changing
FRACTIONAL_SIGMA = 0.9
MIN_GROUP = 20 #minimum group size that is not considered noise

def grouping_sp_dmt(singlepulses):
    """Groups SinglePulse objects based on proximity in time, DM. Takes as input list of singlepulses to group.
        Outputs list of Single Pulse Groups.
    """
    groups = []
    
    for dm, sigma, time, sample, downfact in singlepulses:
#        if dm > 100: #remove this if statement for general case!
#            continue
        grp = group_sp.SinglePulseGroup(dm, sigma, time, sample, downfact)
        groups.append(grp)

    didcombine = True
    while didcombine:
        didcombine = False
        for i, grp1 in enumerate(groups):
            for j in range(len(groups)-1,i,-1):
                if grp1.timeisclose(groups[j]):
                    if grp1.dmisclose(groups[j]):
                        grp1.combine(groups.pop(j))
                        didcombine = True
    return groups

def grouping_sp_t(groups):
    """Groups SinglePulse objects based on proximity in time. Takes as input list of Single Pulse Groups.
        Outputs new list of Single Pulse Groups.
    """
    didcombine = True
    while didcombine:
        didcombine = False
        for i, grp1 in enumerate(groups):
            for j in range(len(groups)-1,i,-1):
                if grp1.timeisclose(groups[j]):
                    grp1.combine(groups.pop(j))
                    didcombine = True
    return groups

def remove_rfi(groups):
    """
    Remove RFI groups based on sigma behavior. Takes as input list of Single Pulse Groups. Outputs new list of Single PU    Pulse Groups, from which RFI groups have been removed. Also outputs list of RFI groups.
    """

#   flag potential RFI groups:
    rfigroups = []
    for j in reversed(range(len(groups))):
        if groups[j].min_dm <= CLOSE_DM:
            for sp in groups[j].singlepulses:
                if sp[0] <= CLOSE_DM:
                    if sp[1] >= FRACTIONAL_SIGMA * groups[j].max_sigma:
                        groups[j].rank = 2 #group marked as "bad" group, RFI
#                        rfigrp = groups.pop(j) #remove rfi group from groups list
                        if groups[j].numpulses >= MIN_GROUP:
                            rfigroups.append(groups[j]) #add rfi group to rfigroups list, unless it's noise (should move this up to before assigning rank 2?)
                        break
        if groups[j].numpulses < MIN_GROUP:
            groups[j].rank = 1 #group marked as noise
        if groups[j].rank == 1 or groups[j].rank == 2:
            groups.pop(j) #remove noise group or rfi group from groups list

    return {'groups':groups, 'rfigroups':rfigroups}

def main():
    import matplotlib.pyplot as plt

#   singlepulses = sp_utils.read_sp_files(filelist)
#    CLOSE_DM = 2 # pc cm-3 #consider changing
#    FRACTIONAL_SIGMA = 0.9
#    MIN_GROUP = 20 #minimum group size that is not considered noise
    MIN_SIGMA = 8
    MAX_DMRANGE = 40 # pc cm-3

    singlepulses = list(sp_utils.read_sp_files())
    singlepulses.sort(cmp = (lambda x,y: cmp(x[1],y[1])),reverse = True) #sort list by time
    
    groups = grouping_sp_dmt(singlepulses)
#    print 'original length of "groups" list:'
#    print len(groups)    

#   flag potential RFI groups:
    rfigroups = remove_rfi(groups)['rfigroups']
    groups = remove_rfi(groups)['groups']

#   regroup good groups based on proximity in time only (compensate for missing middles):
    groups = grouping_sp_t(groups)

#   remove groups that are likely RFI, based on their large span in DM
# to consider: not just DM range, but number of events in extremities.. ie. the distribution of events over DM range. If just a few events grouped in accidentally, shouldn't be discounted.
    for j in reversed(range(len(groups))):
        if groups[j].max_dm-groups[j].min_dm > MAX_DMRANGE:
            groups[j].rank = 2 #group marked as RFI
            rfigrp = groups.pop(j) #remove RFI group from groups list
            rfigroups.append(rfigrp)
# should the above DM span check be done BEFORE grouping_sp_t ?? test this!

#   divide groups into 5 parts (based on number events) to examine sigma behaviour
#    optDMrange = []
#    optDMrangelist = []
    for j in range(len(groups)):
        sigmalist1 = []
        sigmalist2 = []
        sigmalist3 = []
        sigmalist4 = []
        sigmalist5 = []
        groups[j].singlepulses.sort(cmp = (lambda x,y: cmp(x[0],y[0]))) #sort list by incr DM
        L = len(groups[j].singlepulses)
        d = L/5
        for k in range(d):
            sigmalist1.append(groups[j].singlepulses[k][1])
        for k in range(d,2*d):
            sigmalist2.append(groups[j].singlepulses[k][1])
        for k in range(2*d,3*d):
            sigmalist3.append(groups[j].singlepulses[k][1])
        for k in range(3*d,4*d):
            sigmalist4.append(groups[j].singlepulses[k][1])
        for k in range(4*d,5*d):
            sigmalist5.append(groups[j].singlepulses[k][1])
        maxsigma1 = max(sigmalist1) #max sigma value for list 1
        maxsigma2 = max(sigmalist2)
        maxsigma3 = max(sigmalist3)
        maxsigma4 = max(sigmalist4)
        maxsigma5 = max(sigmalist5)
        maxsigma3index = np.argmax(sigmalist3) #
        if maxsigma3 > maxsigma2:
            if maxsigma3 > maxsigma4: #ie. nearest neighbour subgroups both have smaller sigma
                groups[j].rank = 3
                if maxsigma4 > maxsigma5:
                    if maxsigma2 > maxsigma1: #ie. next-nearest subgps have sigma < nearest neighbours
                        groups[j].rank = 4
                        if maxsigma3 > MIN_SIGMA:
                            groups[j].rank = 5
            else: #ie. maxsigma3 <= maxsigma4, allowing for asymmetry:
                if maxsigma2 > maxsigma1:
                    groups[j].rank = 3
                    if maxsigma4 > maxsigma5:
                        groups[j].rank = 4
                        if maxsigma4 > MIN_SIGMA:
                            groups[j].rank = 5
        else: #ie. maxsigma2 >= maxsigma3, allowing for asymmetry:
            if maxsigma2 > maxsigma1:
                if maxsigma3 > maxsigma4:
                    groups[j].rank = 3
                    if maxsigma4 > maxsigma5:
                        groups[j].rank = 4
                        if maxsigma2 > MIN_SIGMA:
                            groups[j].rank = 5


    rank0groups = []
    rank3groups = []
    rank4groups = []
    rank5groups = []
    rank5low_dm_groups = []
    for grp in groups:
        if grp.rank == 0:
            rank0groups.append(grp)
        if grp.rank == 3:
            rank3groups.append(grp)
#            print grp
        if grp.rank == 4:
            rank4groups.append(grp)
#            print grp
        if grp.rank == 5:
            rank5groups.append(grp)
            if grp.min_dm <= 250:
                rank5low_dm_groups.append(grp)

    outfile = open('groups.txt', 'w')

    outfile.write("Number of rank 0 groups: %d \n" %len(rank0groups))
    outfile.write("Number of rank 2 groups: %d \n" %len(rfigroups))
    outfile.write("Number of rank 3 groups: %d \n" %len(rank3groups))
    outfile.write("Number of rank 4 groups: %d \n" %len(rank4groups))
    outfile.write("Number of rank 5 groups: %d \n" %len(rank5groups))
    outfile.write("Number of rank 5 groups with DM < 250: %d \n" %len(rank5low_dm_groups))
    outfile.write("\n")

# print list of events in each group
    for grp in rank5groups:
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')
    
    for grp in rank4groups:
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')

    for grp in rank3groups:
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')

    for grp in rank0groups:
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')

    for grp in rfigroups:
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')

    outfile.close()

    from time import localtime
    datetime_stamp = '%4d-%02d-%02dT%02d-%02d-%02d' %localtime()[:6]
    title = 'grouped_sps'
    plotted = 'rated'
#    plotted = 'largest'
#    ext = 'eps'
    ext = 'png'

#   create four DM vs t plots, splitting up DM in overlapping intervals

#    DM 0-30
    if len(rfigroups) != 0:
        sp_utils2.plot_sp_rated(rfigroups, 'r', 0, 30)
    if len(groups) != 0:
        sp_utils2.plot_sp_rated(groups, 'k', 0, 30)
    if len(rank3groups) != 0:
        sp_utils2.plot_sp_rated(rank3groups, 'g', 0, 30)
    if len(rank4groups) != 0:
        sp_utils2.plot_sp_rated(rank4groups, 'b', 0, 30)
    if len(rank5groups) != 0:
        sp_utils2.plot_sp_rated(rank5groups, 'm', 0, 30)
    plt.savefig('%s_DMs0-30.%s' % (title, ext), dpi = 300)

#   DM 20-110
    if len(rfigroups) != 0:
        sp_utils2.plot_sp_rated(rfigroups, 'r', 20, 110)
    if len(groups) != 0:
        sp_utils2.plot_sp_rated(groups, 'k', 20, 110)
    if len(rank3groups) != 0:
        sp_utils2.plot_sp_rated(rank3groups, 'g', 20, 110)
    if len(rank4groups) != 0:
        sp_utils2.plot_sp_rated(rank4groups, 'b', 20, 110)
    if len(rank5groups) != 0:
        sp_utils2.plot_sp_rated(rank5groups, 'm', 20, 110)
    plt.savefig('%s_DMs20-110.%s' % (title, ext), dpi = 300)

#    DM 100-310
    if len(rfigroups) != 0:
        sp_utils2.plot_sp_rated(rfigroups, 'r', 100, 310)
    if len(groups) != 0:
        sp_utils2.plot_sp_rated(groups, 'k', 100, 310)
    if len(rank3groups) != 0:
        sp_utils2.plot_sp_rated(rank3groups, 'g', 100, 310)
    if len(rank4groups) != 0:
        sp_utils2.plot_sp_rated(rank4groups, 'b', 100, 310)
    if len(rank5groups) != 0:
        sp_utils2.plot_sp_rated(rank5groups, 'm', 100, 310)
    plt.savefig('%s_DMs100-310.%s' % (title, ext), dpi = 300)

#    DM 300-1000
    if len(rfigroups) != 0:
        sp_utils2.plot_sp_rated(rfigroups, 'r', 300, 1000)
    if len(groups) != 0:
        sp_utils2.plot_sp_rated(groups, 'k', 300, 1000)
    if len(rank3groups) != 0:
        sp_utils2.plot_sp_rated(rank3groups, 'g', 300, 1000)
    if len(rank4groups) != 0:
        sp_utils2.plot_sp_rated(rank4groups, 'b', 300, 1000)
    if len(rank5groups) != 0:
        sp_utils2.plot_sp_rated(rank5groups, 'm', 300, 1000)
    plt.savefig('%s_DMs300-1000.%s' % (title, ext), dpi = 300)

if __name__ == '__main__':
    main()
