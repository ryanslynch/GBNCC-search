#!/usr/bin/env python

import os 
import singlepulse.spio as spio
import rfifind #To implement masking for sub-banded files
import singlepulse.spcand as spcand
import singlepulse.read_spd as read_spd
import subprocess
import psrfits
import numpy as np
import glob
import singlepulse.plot_spd as plot_spd
import argparse
from os.path import isfile #To check if subbanded file exists
from operator import itemgetter #To sort candidates list based on DM
import infodata #To read rfifind file

#Modified to sort candidates in order of increasing DM. Subbands the fits file into 128 subbands using psrfits_subband 
#Checks if 5*smearing time because of subbanding at the DM of the previous candidate is less than the pulse width of 
#the current candidate. If yes, uses file subbanded at the previous candidate's DM. Otherwise, subbands at the DM of 
#the current candidate.

def get_obs_info(fitsfilenm):
    """Read in an .inf file to extract observation information.
        Return observation RA, Dec, duration, and source name.
        Function copied from RRATtrap.py 
    """
    inffile = fitsfilenm[:-5] + "_rfifind.inf"
    inf = infodata.infodata(inffile) 
    lofreq = inf.lofreq
    hifreq = (inf.numchan - 1)*inf.chan_width + inf.lofreq
    return {'lofreq': lofreq, 'hifreq': hifreq}
 
def pick_rawdatafile(dm, fitsfilenm, freq_lim_1, freq_lim_2, pulse_width, prevDM, prevsubfile, init_flag):
    center_freq = 350. #MHz
    if init_flag == 1:
        subband = 1 #If this is the first candidate, sub-band at its DM
    else:
        dm_diff = dm - prevDM #error in DM 
        #calculate dispersive smearing in each channel
        t_smear = 8.3e3 * dm_diff * ( center_freq**-3 ) * ( np.abs(freq_lim_1 - freq_lim_2) / 128.) #Number of channels = 128 
        
        if pulse_width > (10 * t_smear):
            #pulse does not get smeared out by subbanding at DM of prev. cand
            subband = 0
            rawdatafile = prevsubfile
        else:
            #pulse gets smeared out by subbanding at DM of prev. cand
            subband = 1
    if subband: #Subband at DM of candidate
        basenm = fitsfilenm[:-5] #basenm = *_0001, fitsfilenm = *_0001.fits 
        subband_file = "%s_subband_%.1f_0001.fits" %(basenm,dm)
        #check if file subbanded at same DM exists, subband at the given DM if it doesn't
        if isfile(subband_file):
            print "Sub-banded file exists at DM of candidate"
        else:
            mask_subband = rfifind.rfifind("%s_rfifind.mask"%(basenm))  
            mask_subband.set_zap_chans(power=1000,plot=False)  
            mask_subband.set_weights_and_offsets()
            mask_subband.write_weights(filename="%s_weights.txt"%(basenm))
            chan, weights = np.loadtxt("%s_weights.txt"%(basenm), unpack = True, skiprows=1, dtype=int)
            rev_weights = weights[::-1]
            rev_weights[819:1640] = 0
            with open('%s_weights.txt'%(basenm), 'r') as f:
	        header = f.readline()
            data = np.column_stack([chan,rev_weights])
            with open('%s_rev_weights.txt'%(basenm), 'w') as f:
                f.write(header)
                np.savetxt(f, data, fmt="%d", delimiter="\t")
            cmd = "psrfits_subband -dm %.1f -nsub 128 -o %s_subband_%.1f -weights %s_rev_weights.txt %s"%(dm,basenm,dm,basenm,fitsfilenm)
            print "executing %s" %cmd
            subprocess.call(cmd, shell=True)
        rawdatafile = subband_file
        prevsubfile = subband_file
        prevDM = dm
    return rawdatafile, prevsubfile, prevDM

def GBNCC_wrapper(txtfile, maskfile, fitsfilenm, path_sp_files):
    """
	The pipeline should pass job.fits_filenm as argument.
    """
    files = spio.get_textfile(txtfile)
    min_rank = 3
    groups = [i for i in range(7) if(i>=min_rank)][::-1]
    numcands=0 # counter for max number of candidates
    loop_must_break = False # dont break the loop unless num of cands >100.
    
    values = []
    lis = [] 
    ranks = []
    for group in groups:
        rank = group+1
        if files[group] != "Number of rank %i groups: 0 "%rank:
            add_values = spio.split_parameters(rank, txtfile) 
            values.extend(add_values)
            lis = np.append(lis, np.where(files == '\tRank:             %i.000000'%rank)[0])
            add_ranks = np.ones(len(add_values)) * rank
            ranks = np.append(ranks, add_ranks)
     
    if len(values) > 100:
        values = values[0:100]
        lis = lis[0:100]

    #Sort candidates based on DM
    zip_list = zip(values, lis, ranks)
    zip_list = sorted(zip_list, key=itemgetter(0,0))
    values = [x[0] for x in zip_list]
    lis = [x[1] for x in zip_list]
    ranks = [x[2] for x in zip_list]
    
    basename = fitsfilenm[:-5]
    #generate subbanded file at a DM of 0 to extract observation parameters
    cmd = "psrfits_subband -dm 0.0 -nsub 128 -o %s_subband_0.0 %s"%(basename,fitsfilenm)
    print "executing %s" %cmd
    subprocess.call(cmd, shell=True)
    subfilenm = basename + "_subband_0.0_0001.fits"
    subfile = psrfits.PsrfitsFile(subfilenm)
    
    for ii in range(len(values)):
        dm_list, time_list, dm_arr, sigma_arr, width_arr = spio.read_RRATrap_info(txtfile, lis[ii], int(ranks[ii]))
        wrapper_cand = spcand.params() 
        wrapper_cand.read_from_file(values[ii], subfile.tsamp, subfile.specinfo.N, \
                                       get_obs_info(fitsfilenm)['hifreq'], get_obs_info(fitsfilenm)['lofreq'], \
                                       subfile, loc_pulse=0.5, dedisp = True,\
                                       scaleindep = None, zerodm = None, mask = None,\
                                       barytime=True, nsub = None, bandpass_corr = False) 
        temp_filename = basename+"_rank_%i"%int(ranks[ii])
        if ii == 0: #check if index is 0
            correct_rawdatafile, prevsubfile, prevDM  = pick_rawdatafile(wrapper_cand.subdm, fitsfilenm, \
                                                                get_obs_info(fitsfilenm)['hifreq'], get_obs_info(fitsfilenm)['lofreq'], \
                                                                wrapper_cand.pulse_width, prevDM = 0, prevsubfile = '', \
                                                                init_flag = 1)
        else:
            correct_rawdatafile, prevsubfile, prevDM  = pick_rawdatafile(wrapper_cand.subdm, fitsfilenm, \
                                                                get_obs_info(fitsfilenm)['hifreq'], get_obs_info(fitsfilenm)['lofreq'], \
                                                                wrapper_cand.pulse_width, prevDM = prevDM, \
                                                                prevsubfile=prevsubfile, init_flag = 0)
            
        
        path = os.path.dirname(os.path.abspath(__file__))
        cmd = "python " + path + "/make_spd.py --use_manual_params --subdm %f --nsub %d" \
              	      " --dm %f -T %f -t %f --width-bins %d --downsamp %d --show-spec"\
                      " --noplot --notopo -o %s %s " %(wrapper_cand.subdm, \
                                                 wrapper_cand.nsub, wrapper_cand.dm, \
                                                 wrapper_cand.topo_start_time, wrapper_cand.duration, \
                                                 wrapper_cand.width_bins, wrapper_cand.downsamp,\
                                                 temp_filename, correct_rawdatafile)
        print "executing %s" %cmd
        subprocess.call(cmd, shell=True)       
        # Add additional information to the header information array
        text_array = np.array([correct_rawdatafile, subfile.specinfo.telescope, \
                           subfile.specinfo.ra_str, subfile.specinfo.dec_str, \
                           subfile.specinfo.start_MJD[0], int(ranks[ii]), \
                           wrapper_cand.nsub, wrapper_cand.nbins, \
                           wrapper_cand.subdm, wrapper_cand.sigma, wrapper_cand.sample_number, \
                           wrapper_cand.duration, wrapper_cand.width_bins, wrapper_cand.pulse_width, \
                           subfile.tsamp, subfile.specinfo.T, wrapper_cand.topo_start_time])
                
        temp_filename +="_DM%.1f_%.1fs"%(wrapper_cand.subdm, wrapper_cand.topo_start_time)
        spd = read_spd.spd(temp_filename+'.spd')
        spd.man_params = None 
        text_array = np.append(text_array, spd.waterfall_start_time)
        text_array = np.append(text_array, spd.waterfall_tsamp)
        text_array = np.append(text_array, spd.waterfall_prededisp_nbins)
        text_array = np.append(text_array, spd.min_freq)
        text_array = np.append(text_array, spd.max_freq)
        text_array = np.append(text_array, spd.sweep_duration)
        text_array = np.append(text_array, spd.sweep_start_time)
        text_array = np.append(text_array, spd.bary_pulse_peak_time)
        text_array = np.append(text_array, spd.man_params)

        with open(temp_filename+".spd", 'wb') as f:
            np.savez_compressed(f, \
                                        Data_dedisp_nozerodm = spd.data_nozerodm_dedisp,\
                                        Data_dedisp_zerodm = spd.data_zerodm_dedisp,\
                                        Data_nozerodm = spd.data_nozerodm,\
                                        delays_nozerodm = spd.dmsweep_delays, \
                                        freqs_nozerodm = spd.dmsweep_freqs,\
                                        Data_zerodm = spd.data_zerodm, \
                                        dm_arr= map(np.float16, dm_arr),\
                                        sigma_arr = map(np.float16, sigma_arr), \
                                        width_arr =map(np.uint8, width_arr),\
                                        dm_list= map(np.float16, dm_list), \
                                        time_list = map(np.float16, time_list), \
                                        text_array = text_array)
        plot_spd.plot(temp_filename+".spd", glob.glob(path_sp_files+'/*.singlepulse'), maskfile, \
                              outfile=basename, just_waterfall=False, \
                              integrate_spec=True, integrate_ts=True, \
                              disp_pulse=False, bandpass_corr = True, tar = None)
        numcands+= 1
        print 'Finished sp_candidate : %i'%numcands
        if numcands >= 100:    # Max number of candidates to plot 100.
            loop_must_break = True
            break

if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Code to waterfall single pulse candidates grouped"\
                                     "and ranked by RRATtrap")
    parser.add_argument("--groups", dest="groupsfile", default=None,
                        help="groups.txt file containing RRATtrap candidate")
    parser.add_argument("--mask", dest="maskfile",default=None,
                        help="Mask file produced by rfifind")
    parser.add_argument("--fits", dest ="fitsfile",default=None,
                        help="A psrfits file from the GBNCC survey")
    parser.add_argument("--dir", dest ="workdir",default=None,
                        help="Path of the working directory")
    args = parser.parse_args()

    GBNCC_wrapper(args.groupsfile, args.maskfile, args.fitsfile, args.workdir) 
