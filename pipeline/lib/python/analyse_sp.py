#!/usr/bin/env python

"""
Analyse all single pulse files in a given beam by grouping them, ranking the groups, and plotting them in colours corresponding to their rating.
Modifications to analyse_sp3.py script which aim to improve the grouping and rating algorithms.
Modifications to analyse_sp_v2.py which aim to better classify RFI at high DMs and avoid it being rated highly.
Chen Karako
Nov. 4, 2011

Making changes to analyse_sp_v3.py in order to fix bugs and eventually integrate code into GBNCC pipelines.
Oct 2, 2012

Combining Felix's modifications (speed ups, rank 6, memory minimization) with my own. Uses inf files to get obs info (source position, observation length).
Aug 26, 2013

Increasing min group size to get rid of noise.
Oct 16, 2013
"""
import numpy as np
import sp_utils
import group_sp
import sys
import time
#import memory
#from operator import itemgetter

CLOSE_DM = 2 # pc cm-3 #consider changing
FRACTIONAL_SIGMA = 0.9
MIN_GROUP = 100 #minimum group size that is not considered noise
TIME_THRESH = 0.2 # What value would be optimal?
DM_THRESH = 2
MIN_SIGMA = 8


def grouping_sp_dmt(singlepulses):
    """Groups SinglePulse objects based on proximity in time, DM. Takes as input list of singlepulses to group.
        Outputs list of Single Pulse Groups.
    """
    groups = []
    for dm, sigma, time, sample, downfact in singlepulses:
#        if dm > 100: #remove this if statement for general case!
#            continue
        grp = group_sp.SinglePulseGroup(dm, sigma, time+0.004*dm, sample, downfact) #The 0.004*dm is to compensate for the slope
        groups.append(grp)
    groups.sort(key=lambda group: group.center_time) #Sort by time
    
    didcombine = True
    while didcombine:
        didcombine = False
        for i, grp1 in enumerate(groups):
            j=i+1
            while (j<len(groups) and groups[i].center_time+TIME_THRESH > groups[j].center_time): #Only look at groups that are close in time
               if grp1.dmisclose(groups[j]): # Consider switching this if statement with the next one if TIME_THRESH is bigger than 0.2
                    if grp1.timeisclose(groups[j]):
                        grp1.combine(groups.pop(j))
                        didcombine = True
               j=j+1

            j=i-1
            while ( j>-1 and i<len(groups) and groups[i].center_time-TIME_THRESH < groups[j].center_time): # We must check the length of i, since some elements might have popped up in the last loop
               if grp1.dmisclose(groups[j]): # Consider switching this if statement with the next one if TIME_THRESH is bigger than 0.2
                     if grp1.timeisclose(groups[j]): 
                        grp1.combine(groups.pop(j))
                        didcombine = True
               j=j-1

    return groups
    """
   # This method is in general less efficient and precise than sorting by time, but might be useful for particular cases
    groups=sorted(groups, key=lambda group: group.max_dm) #Sort by DM
  
    didcombine = True
    while didcombine:
        didcombine = False
        for i, grp1 in enumerate(groups):
            j = i + 1
            while (j < len(groups) and groups[i].max_dm+DM_THRESH>groups[j].max_dm):
                if grp1.timeisclose(groups[j]):
                    if grp1.dmisclose(groups[j]): 
                        grp1.combine(groups.pop(j))
                        didcombine = True
                j = j + 1

            j = i - 1
            while ( j > -1 and i<len(groups) and groups[i].max_dm-DM_THRESH<groups[j].max_dm): # We must check the length of i, since some elements might have poped up in the last loop
                if grp1.timeisclose(groups[j]):
                    if grp1.dmisclose(groups[j]):  
                        grp1.combine(groups.pop(j))
                        didcombine = True
                j = j - 1
#    printTime(start_time)
    return groups   
    """
    """
def is_rank6(groups,j):
        
        #This determines if a a group is a potential rank 6 candidate (by looking at the max sigma only) and returns true or false.
       
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
    #    avgsigma1 = sum(sigmalist1)/len(sigmalist1) #Not using for the moment
     #   avgsigma2 = sum(sigmalist2)/len(sigmalist2)
    #    avgsigma3 = sum(sigmalist3)/len(sigmalist3)
     #   avgsigma4 = sum(sigmalist4)/len(sigmalist4)
     #   avgsigma5 = sum(sigmalist5)/len(sigmalist5)
        if maxsigma3 > maxsigma2 :
            if maxsigma3 > maxsigma4: #ie. nearest neighbour subgroups both have smaller sigma
                if maxsigma4 > maxsigma5 :
                    if maxsigma2 > maxsigma1: #ie. next-nearest subgps have sigma < nearest neighbours
                        if maxsigma3 > MIN_SIGMA:
                            return False
            if maxsigma3 <= maxsigma4 : #ie. maxsigma3 <= maxsigma4, allowing for asymmetry:
                if maxsigma2 > maxsigma1 :
                    if maxsigma4 > maxsigma5 :
                        if maxsigma4 > MIN_SIGMA:
                            return False
        if maxsigma2 >= maxsigma3  : #ie. maxsigma2 >= maxsigma3, allowing for asymmetry:
            if maxsigma2 > maxsigma1 :
                if maxsigma3 > maxsigma4 :
                    if maxsigma4 > maxsigma5  :
                        if maxsigma2 > MIN_SIGMA:
                            return False
        return True
  """
def grouping_rfi(groups, rfigroups):
  """
  Groups together close groups of rfi and consider as rfi other groups that may be close to groups of rfi.
  Takes as input the rfigroup and the group, and return the new groups.
  """
  # First we want to group the RFI with a similar time together
  for i, grp1 in enumerate(rfigroups):
     for j in range(len(rfigroups)-1,i,-1):
         if grp1.dmisclose(rfigroups[j],10) and grp1.timeisclose(rfigroups[j],0.2):
             rfigroups[j].combine(rfigroups.pop(i))
             break
  k=0
  # If a group is very close to a group of rfi, set it as rfi  
  while k<5: # We loop five times, in case bigger rfi groups are formed
    for i, grp1 in enumerate(groups):
        for j, grp2 in enumerate(rfigroups):
            if grp1.dmisclose(grp2,10) and grp1.timeisclose(grp2) and grp2.numpulses>50 and grp2.max_dm>2: # We do not want to group with noise or low DM rfi
               groups[i].rank=2 # Set as rfi
               grp2.combine(groups.pop(i))
               break
    k=k+1
  return {'groups':groups, 'rfigroups':rfigroups}

"""  This part could be used to group everything under 15 DM as RFI if there is already a lot of interference under 15 DM. It also removes groups that have their max_sigma under 5 DM and a max DM over 
  num_rfi=0
  for j in reversed(range(len(groups))): # Remove groups that have a max_sigma at a DM under 5 and a max DM over 20, since they are most likely rfi.
     max_sigma=0
     center_dm=0
     for i in range(len(groups[j].singlepulses)):
         if groups[j].singlepulses[i][1]>max_sigma:
               max_sigma=groups[j].singlepulses[i][1]
               center_dm=groups[j].singlepulses[i][0]
     if center_dm < 5 and groups[j].max_dm> 15:
         groups[j].rank == 2
         rfigroups.append(groups.pop(j))
         num_rfi=num_rfi+1
  if num_rfi>10:
     for j in reversed(range(len(groups))): # If we have a lot of rfi under 15 DM, we consider everything under 15 DM as rfi
       if groups[j].max_dm< 15:
          groups[j].rank == 2
          rfigroups.append(groups.pop(j))
  return {'groups':groups, 'rfigroups':rfigroups}
"""
def grouping_sp_t(groups):
    """Groups SinglePulse objects based on proximity in time. Takes as input list of Single Pulse Groups.
        Outputs new list of Single Pulse Groups.
    """
    didcombine = True
    while didcombine:
        didcombine = False
        for i, grp1 in enumerate(groups):
            for j in range(len(groups)-1,i,-1):
                if grp1.timeisclose(groups[j]) and grp1.dmisclose(groups[j],10): # We check if two events have a very similar time and a DM difference under 10
                    grp1.combine(groups.pop(j))
                    didcombine = True
    return groups

def grouping_diagonal(groups):
    """Groups SinglePulse objects based on proximity in time and DM in order to compensate for the slope. Takes as input list of Single Pulse Groups.
        Outputs new list of Single Pulse Groups.
    """
    didcombine = True
    for i, grp1 in enumerate(groups):
        count = 0
        for j in range(len(groups)-1,i,-1):
            if grp1.timeisclose(groups[j],4) and groups[j].min_dm>15: # We check how many events have a similar time
                count = count + 1
        if count > 15: # If there is more than 15 groups within a 4 seconds interval we combine those that are close in DM and time 
            while didcombine:
                didcombine = False
                for i, grp1 in enumerate(groups):
                    for j in range(len(groups)-1,i,-1):
                        if grp1.timeisclose(groups[j],0.4) and grp1.dmisclose(groups[j],groups[j].max_dm/10):
                            grp1.combine(groups.pop(j))
                            didcombine = True
    return groups


def remove_rfi(groups):
    """
    Remove RFI groups based on sigma behavior. Takes as input list of Single Pulse Groups. Outputs new list of Single Pulse Groups, from which RFI groups have been removed. Also outputs list of RFI groups.
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
    start_time = time.time()    
#    print "memory is " + str(memory.memory()/1000000)
#   singlepulses = sp_utils.read_sp_files(filelist)
    MAX_DMRANGE = 80 # pc cm-3
#    print "memory is " + str(memory.memory()/1000000)
    singlepulses = list(sp_utils.read_sp_files())         
    singlepulses.sort(cmp = (lambda x,y: cmp(x[1],y[1])),reverse= True)

    elapsed_time = time.time() - start_time
#    print "after extracting", elapsed_time
    groups = grouping_sp_dmt(singlepulses)
#    print "memory before deleting is " + str(memory.memory()/1000000)
    elapsed_time = time.time() - start_time
#    print "after grouping", elapsed_time

 
#   flag potential RFI groups:
    rfigroups = remove_rfi(groups)['rfigroups']
    groups = remove_rfi(groups)['groups']

  
#   regroup good groups based on proximity in time only (compensate for missing middles):
    groups = grouping_diagonal(groups)
    groups = grouping_sp_t(groups)
#   remove groups that are likely RFI, based on their large span in DM
# to consider: not just DM range, but number of events in extremities.. ie. the distribution of events over DM range. If just a few events grouped in accidentally, shouldn't be discounted.
    for j in reversed(range(len(groups))):
        if groups[j].max_dm-groups[j].min_dm > MAX_DMRANGE:
            groups[j].rank = 2 #group marked as RFI
            rfigrp = groups.pop(j) #remove RFI group from groups list
            rfigroups.append(rfigrp)
# should the above DM span check be done BEFORE grouping_sp_t ?? test this! (in most cases not)
# Group rfi with very close groups
    temp = grouping_rfi(groups, rfigroups)
    groups =temp['groups']
    rfigroups =temp['rfigroups']
    del temp
#    print "memory usage is " + str(memory.memory()/1000000)
#   divide groups into 5 parts (based on number events) to examine sigma behaviour
    sumslope=0
    for j in range(len(groups)):
        sigmalist1 = []
        sigmalist2 = []
        sigmalist3 = []
        sigmalist4 = []
        sigmalist5 = []
        groups[j].singlepulses.sort(cmp = (lambda x,y: cmp(x[0],y[0]))) #sort list by incr DM
        L = len(groups[j].singlepulses)
        d = L/5
        sumdm=0
        sumdm2=0
        sumtime=0
        sumtime2=0
        for k in range(d):
            sigmalist1.append(groups[j].singlepulses[k][1])
            sumdm=sumdm+groups[j].singlepulses[k][0]
            sumtime=sumtime+groups[j].singlepulses[k][2]
        for k in range(d,2*d):
            sigmalist2.append(groups[j].singlepulses[k][1])
        for k in range(2*d,3*d):
            sigmalist3.append(groups[j].singlepulses[k][1])
        for k in range(3*d,4*d):
            sigmalist4.append(groups[j].singlepulses[k][1])
        for k in range(4*d,5*d):
            sumdm2=sumdm2+groups[j].singlepulses[k][0]
            sumtime2=sumtime2+groups[j].singlepulses[k][2]
            sigmalist5.append(groups[j].singlepulses[k][1])

        maxsigma1 = max(sigmalist1) #max sigma value for list 1
        maxsigma2 = max(sigmalist2)
        maxsigma3 = max(sigmalist3)
        maxsigma4 = max(sigmalist4)
        maxsigma5 = max(sigmalist5)
        avgsigma1 = sum(sigmalist1)/len(sigmalist1) #max average sigma value for list 1
        avgsigma2 = sum(sigmalist2)/len(sigmalist2)
        avgsigma3 = sum(sigmalist3)/len(sigmalist3)
        avgsigma4 = sum(sigmalist4)/len(sigmalist4)
        avgsigma5 = sum(sigmalist5)/len(sigmalist5)
        maxsigma3index = np.argmax(sigmalist3) #
        maxsigma=max(maxsigma1,maxsigma2,maxsigma3,maxsigma4,maxsigma5) # The largest maxsigma
        minsigma=min(maxsigma1,maxsigma2,maxsigma3,maxsigma4,maxsigma5) # The smallest maxsigma
        maxavgsigma=max(avgsigma1,avgsigma2,avgsigma3,avgsigma4,avgsigma5)
        minavgsigma=min(avgsigma1,avgsigma2,avgsigma3,avgsigma4,avgsigma5)
        if maxavgsigma<1.05*minavgsigma: # If the sigma is pretty much constant, then it is rfi
           groups[j].rank=2
        minsigma=min(maxsigma1,maxsigma2,maxsigma3,maxsigma4,maxsigma5) # The smallest maxsigma
        if maxsigma3 > maxsigma2 :
            if maxsigma3 > maxsigma4: #ie. nearest neighbour subgroups both have smaller sigma
                groups[j].rank = 3
                if maxsigma4 > maxsigma5:
                    if maxsigma2 > maxsigma1: #ie. next-nearest subgps have sigma < nearest neighbours
                        groups[j].rank = 4
                        if maxsigma3 > MIN_SIGMA:  # We want the largest maxsigma to be at least 1.25 times bigger than the smallest
                            groups[j].rank = 5
                            if avgsigma3 > avgsigma1 and avgsigma3 > avgsigma5 and maxsigma>1.15*minsigma:
                               groups[j].rank = 6 
            if maxsigma3 <= maxsigma4: #ie. maxsigma3 <= maxsigma4, allowing for asymmetry:
                if maxsigma2 > maxsigma1 :
                    groups[j].rank = 3
                    if maxsigma4 > maxsigma5 :
                        groups[j].rank = 4
                        if maxsigma4 > MIN_SIGMA :
                            groups[j].rank = 5
                            if avgsigma4 > avgsigma1  and avgsigma4 > avgsigma5 and maxsigma>1.15*minsigma:
                               groups[j].rank = 6                             
        if maxsigma2 >= maxsigma3 : #ie. maxsigma2 >= maxsigma3, allowing for asymmetry:
            if maxsigma2 > maxsigma1:
                if maxsigma3 > maxsigma4 :
                    groups[j].rank = 3
                    if maxsigma4 > maxsigma5  :
                        groups[j].rank = 4
                        if maxsigma2 > MIN_SIGMA:
                            groups[j].rank = 5
                            if avgsigma2 >= avgsigma1 and avgsigma2 > avgsigma5 and maxsigma>1.15*minsigma:
                               groups[j].rank = 6 
    rank0groups = []
    rank3groups = []
    rank4groups = []
    rank5groups = []
    rank6groups = []
    for grp in reversed(groups):
        if grp.rank == 3:
            rank3groups.append(grp)
        elif grp.rank == 4:
            rank4groups.append(grp)
        elif grp.rank == 5:
            rank5groups.append(grp)
        elif grp.rank == 6:
            rank6groups.append(grp)
        else:
           rank0groups.append(grp)
        groups.pop() # not sure if actually frees memory, try del groups.pop()
#    print "memory after ranks is " + str(memory.memory()/1000000)
    outfile = open('groups.txt', 'w')

    outfile.write("Number of rank 0 groups: %d \n" %len(rank0groups))
    outfile.write("Number of rank 2 groups: %d \n" %len(rfigroups))
    outfile.write("Number of rank 3 groups: %d \n" %len(rank3groups))
    outfile.write("Number of rank 4 groups: %d \n" %len(rank4groups))
    outfile.write("Number of rank 5 groups: %d \n" %len(rank5groups))
    outfile.write("Number of rank 6 groups: %d \n" %len(rank6groups))
    outfile.write("\n")
    rank6groupsA=[]
    rank6groupsB=[]
    rank6groupsC=[]
    rank6groupsD=[]
    rank5groupsA=[]
    rank5groupsB=[]
    rank5groupsC=[]
    rank5groupsD=[]
    rank4groupsA=[]
    rank4groupsB=[]
    rank4groupsC=[]
    rank4groupsD=[]
    rank3groupsA=[]
    rank3groupsB=[]
    rank3groupsC=[]
    rank3groupsD=[]
    rank0groupsA=[]
    rank0groupsB=[]
    rank0groupsC=[]
    rank0groupsD=[]
    rfigroupsA=[]
    rfigroupsB=[]
    rfigroupsC=[]
    rfigroupsD=[]
# print list of events in each group
    for grp in reversed(rank6groups):
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')
        if grp.min_dm<=30:
           rank6groupsA.append(grp)
        if grp.min_dm<=110 and grp.max_dm>=20:
           rank6groupsB.append(grp)
        if grp.min_dm<=310 and grp.max_dm>=100:
           rank6groupsC.append(grp)
        if grp.max_dm>=300:
           rank6groupsD.append(grp)
        rank6groups.pop()
                   

    for grp in reversed(rank5groups):
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')
        if grp.min_dm<=30:
           rank5groupsA.append(grp)
        if grp.min_dm<=110 and grp.max_dm>=20:
           rank5groupsB.append(grp)
        if grp.min_dm<=310 and grp.max_dm>=100:
           rank5groupsC.append(grp)
        if grp.max_dm>=300:
           rank5groupsD.append(grp)
        rank5groups.pop()
    
    for grp in reversed(rank4groups):
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')
        if grp.min_dm<=30:
           rank4groupsA.append(grp)
        if grp.min_dm<=110 and grp.max_dm>=20:
          rank4groupsB.append(grp)
        if grp.min_dm<=310 and grp.max_dm>=100:
           rank4groupsC.append(grp)
        if grp.max_dm>=300:
           rank4groupsD.append(grp)
        rank4groups.pop()

    for grp in reversed(rank3groups):
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')
        if grp.min_dm<=30:
           rank3groupsA.append(grp)
        if grp.min_dm<=110 and grp.max_dm>=20:
           rank3groupsB.append(grp)
        if grp.min_dm<=310 and grp.max_dm>=100:
           rank3groupsC.append(grp)
        if grp.max_dm>=300:
           rank3groupsD.append(grp)
        rank3groups.pop()
         
    for grp in reversed(rank0groups):
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')
        if grp.min_dm<=30:
           rank0groupsA.append(grp)
        if grp.min_dm<=110 and grp.max_dm>=20:
           rank0groupsB.append(grp)
        if grp.min_dm<=310 and grp.max_dm>=100:
           rank0groupsC.append(grp)
        if grp.max_dm>=300:
           rank0groupsD.append(grp)
        rank0groups.pop()

    for grp in reversed(rfigroups):
        outfile.write(str(grp) + '\n') #print group summary
        outfile.write('\n')
        outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
        for sp in grp.singlepulses:
            outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
        outfile.write('\n')
        if grp.min_dm<=30:
           rfigroupsA.append(grp)
        if grp.min_dm<=110 and grp.max_dm>=20:
           rfigroupsB.append(grp)
        if grp.min_dm<=310 and grp.max_dm>=100:
           rfigroupsC.append(grp)
        if grp.max_dm>=300:
           rfigroupsD.append(grp)
        rfigroups.pop()


    outfile.close()
    title = 'grouped_sps'
    plotted = 'rated'
    ext = 'png'
    elapsed_time = time.time() - start_time
#    print "before plots", elapsed_time
#   create four DM vs t plots, splitting up DM in overlapping intervals 
#    print "memory after merge is " + str(memory.memory()/1000000)   
#    DM 0-30
    if len(rfigroupsA) != 0:
        sp_utils.plot_sp_rated(rfigroupsA, 'r', 0, 30)
    if len(rank0groupsA) != 0:
        sp_utils.plot_sp_rated(rank0groupsA, 'k', 0, 30) ### AND RANK 0
    if len(rank3groupsA) != 0:
        sp_utils.plot_sp_rated(rank3groupsA, 'g', 0, 30)
    if len(rank4groupsA) != 0:
        sp_utils.plot_sp_rated(rank4groupsA, 'b', 0, 30)
    if len(rank5groupsA) != 0:
        sp_utils.plot_sp_rated(rank5groupsA, 'm', 0, 30)
    if len(rank6groupsA) != 0:
        sp_utils.plot_sp_rated(rank6groupsA, 'c', 0, 30)
    plt.savefig('%s_DMs0-30.%s' % (title, ext), dpi = 300)
    plt.close()
    del rfigroupsA
    del rank6groupsA
    del rank5groupsA
    del rank4groupsA
    del rank3groupsA
    del rank0groupsA
#   DM 20-110
    if len(rfigroupsB) != 0:
        sp_utils.plot_sp_rated(rfigroupsB, 'r', 20, 110)
    if len(rank0groupsB) != 0:
        sp_utils.plot_sp_rated(rank0groupsB, 'k', 20, 110)
    if len(rank3groupsB) != 0:
        sp_utils.plot_sp_rated(rank3groupsB, 'g', 20, 110)
    if len(rank4groupsB) != 0:
        sp_utils.plot_sp_rated(rank4groupsB, 'b', 20, 110)
    if len(rank5groupsB) != 0:
        sp_utils.plot_sp_rated(rank5groupsB, 'm', 20, 110)
    if len(rank6groupsB) != 0:
        sp_utils.plot_sp_rated(rank6groupsB, 'c', 20, 110)
    plt.savefig('%s_DMs20-110.%s' % (title, ext), dpi = 300)
    plt.close()
    del rfigroupsB
    del rank6groupsB
    del rank5groupsB
    del rank4groupsB
    del rank3groupsB
    del rank0groupsB
#    DM 100-310
    if len(rfigroupsC) != 0:
        sp_utils.plot_sp_rated(rfigroupsC, 'r', 100, 310)
    if len(rank0groupsC) != 0:
        sp_utils.plot_sp_rated(rank0groupsC, 'k', 100, 310)
    if len(rank3groupsC) != 0:
        sp_utils.plot_sp_rated(rank3groupsC, 'g', 100, 310)
    if len(rank4groupsC) != 0:
        sp_utils.plot_sp_rated(rank4groupsC, 'b', 100, 310)
    if len(rank5groupsC) != 0:
        sp_utils.plot_sp_rated(rank5groupsC, 'm', 100, 310)
    if len(rank6groupsC) != 0:
        sp_utils.plot_sp_rated(rank6groupsC, 'c', 100, 310)
    plt.savefig('%s_DMs100-310.%s' % (title, ext), dpi = 300)
    plt.close()
    del rfigroupsC
    del rank6groupsC
    del rank5groupsC
    del rank4groupsC
    del rank3groupsC
    del rank0groupsC

#    DM 300-1000
    if len(rfigroupsD) != 0:
        sp_utils.plot_sp_rated(rfigroupsD, 'r', 300, 1000)
    if len(rank0groupsD) != 0:
        sp_utils.plot_sp_rated(rank0groupsD, 'k', 300, 1000)
    if len(rank3groupsD) != 0:
        sp_utils.plot_sp_rated(rank3groupsD, 'g', 300, 1000)
    if len(rank4groupsD) != 0:
        sp_utils.plot_sp_rated(rank4groupsD, 'b', 300, 1000)
    if len(rank5groupsD) != 0:
        sp_utils.plot_sp_rated(rank5groupsD, 'm', 300, 1000)
    if len(rank6groupsD) != 0:
        sp_utils.plot_sp_rated(rank6groupsD, 'c', 300, 1000)
    plt.savefig('%s_DMs300-1000.%s' % (title, ext), dpi = 300)
    plt.close()
    del rfigroupsD
    del rank6groupsD
    del rank5groupsD
    del rank4groupsD
    del rank3groupsD
    del rank0groupsD

    elapsed_time = time.time() - start_time
#    print "end", elapsed_time
#    print "memory usage is " + str(memory.memory()/1000000)
if __name__ == '__main__':
    main()
