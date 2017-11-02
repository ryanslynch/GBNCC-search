def get_ffa_folding_command(cand, obs, ddplans, maskfilenm):
    """
    get_ffa_folding_command(cand, obs):
        Return a command for prepfold for folding the subbands using
            an obs_info instance, and a candidate instance that 
            describes the observations and searches.
    """
    # Folding rules are same as those for folding accel cands.
    outfilenm = obs.basefilenm+"_DM%s_ffa"%cand.DMstr

    # Note:  the following calculations should probably only be done once,
    #        but in general, these calculation are effectively instantaneous
    #        compared to the folding itself
    hidms = [x.lodm for x in ddplans[1:]] + [2000]
    dfacts = [x.downsamp for x in ddplans]
    numsubs = [x.numsub for x in ddplans]
    for hidm, dfact,numsub in zip(hidms, dfacts, numsubs):
        if cand.DM < hidm:
            downsamp = dfact
            nsub = numsub
            break
    if downsamp==1:
        fitsfile = obs.fits_filenm
    else:
        fitsfile =  obs.dsbasefilenm+"_DS%d%s"%\
                    (downsamp,obs.fits_filenm[obs.fits_filenm.rfind("_"):])
    #p = cand.p
    if cand.p < 0.5:
        N = 100
        npart = 40
        otheropts = "-pstep 1 -pdstep 2 -dmstep 1 -nodmsearch"
    elif cand.p < 2.0:
        N = 200
        npart = 40
        otheropts = "-nosearch -slow" 
    elif cand.p < 5.0:
        N = 200
        npart = 30
        otheropts = "-nosearch -slow" 
    elif cand.p < 10.0:
        N = 200
        npart = 20
        otheropts = "-nosearch -slow" 
    else:
        N = 200
        npart = 10
        otheropts = "-nosearch -slow"

    return "prepfold -noxwin -mask %s -dm %.2f -p %f -o %s " \
                "-nsub %d -npart %d %s -n %d %s" % \
           (maskfilenm,cand.DM, cand.p, outfilenm,
            nsub, npart, otheropts, N, fitsfile)

