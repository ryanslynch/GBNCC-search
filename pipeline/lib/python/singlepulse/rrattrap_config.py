# /usr/bin/env python
from scipy.special import erf
import numpy as np

tel = 'GBNCC'
#tel = 'PALFA'

CLOSE_DM = 2 # pc cm-3
# MIN_GROUP, DM_THRESH, TIME_THRESH will change later on depending on the DDplan.
MIN_GROUP = 30 #minimum group size that is not considered noise
TIME_THRESH = 0.1
MIN_SIGMA = 8.0
PLOT = True
PLOTTYPE = 'pgplot' # 'pgplot' or 'matplotlib'
RANKS_TO_WRITE = [2,0,3,4,5,6,7]
RANKS_TO_PLOT = [2,3,4,5,6,7]

#Inputs for varying cluster size based on DM, S/N and pulse width
if tel == 'GBNCC':
    DM_step = np.array([0.01,0.02,0.03,0.05,0.1,0.3,0.5])
    low_DM = np.array([0.0,70.52,125.56,217.36,390.76,778.36,1882.36])
    DM_THRESH = 0.05
    ####factors that mutliply DM_THRESH and TIME_THRESH 
    dm_fac_arr = np.array([1,2,3,5,10,30,50])
    t_fac_arr = np.array([1,1,1,2,4,6,8])
elif tel == 'PALFA':
    DM_step = np.array([0.1,0.3,0.3,0.5,0.5,1.0,2.0,3.0,5.0])
    low_DM = np.array([0.0,212.8,443.2,534.4,876.4,990.4,1826.4,3266.4,5546.4])
    DM_THRESH = 0.5
    ####factors that mutliply DM_THRESH and TIME_THRESH
    ####To search 5 neighbouring DM trials#####
    #dm_fac_arr = np.array([1,3,3,5,5,10,20,30,50]) 
    ####To search 5 neighbouring DM trials at low DMs and 3 DM trials for pulses with DM > 3266####
    dm_fac_arr = np.array([1,3,3,5,5,10,20,18,30])
    t_fac_arr = np.array([1,2,3,5,6,10,10,10,10])          

def use_dmplan(DM):                                                                                              
    """ Sets a factor which multiplies the DMthreshold and time_threshold.                                                                      
        This makes the DM_THRESH and TIME_THRESH depend on the DM instead of having fixed                                                       
        values throughout. This helps at higher DMs where the DM step size is > 0.5 pc cm-3.                                                    
    """                                                                                                                                             
    idx = np.digitize([DM],low_DM)[0] - 1 #choosing the DM step size for the given DM
    dm_fac = dm_fac_arr[idx]
    t_fac = t_fac_arr[idx]    
    return dm_fac, t_fac

def t_chan(DM, c_freq_MHz, chan_width):
    #Get dispersive smearing within each channel in seconds 
    return 8.297616e-6 * (chan_width) * (c_freq_MHz*1e-3)**-3 * (DM)
    

def vary_group_size(S_N_peak,W_B_s,DM,c_freq_MHz, dt, chan_width, BW_MHz):
    """S_N_peak : Peak S/N of the group
       W_B_s : Broadened Pulse Width in s 
       DM : Optimal DM of pulse
    """
    if S_N_peak < MIN_SIGMA:
        min_group = MIN_GROUP
    else: #Vary Group size
        
        #Calculate DM step size
        idx = np.digitize([DM],low_DM)[0] - 1 

        #Calculate W_i**2 + t_scatt**2 by subtracting sampling time and expected dispersive smearing
        sum_quad = W_B_s**2 - dt**2 - t_chan(DM,c_freq_MHz,chan_width)**2  
        if sum_quad > 0:
            W_i_s = np.sqrt(sum_quad)
        else:
            W_i_s = W_B_s 
       	W_i_ms = W_i_s * 1e3 
    
    	S_peak = S_N_peak * np.sqrt(W_B_s) / W_i_s

    	#Getting S for an error in DM
    	cluster_size = 100
    	DM_err = np.arange(0,DM_step[idx]*0.5*cluster_size+1,DM_step[idx])

    	zeta = 6.91e-3 * DM_err * BW_MHz / (W_i_ms * (c_freq_MHz*1e-3)**3)  #Cordes and McLaughlin Eq13
    	S_meas = np.zeros_like(zeta)
    	S_meas[np.where(zeta != 0)] = S_peak * 0.5 * np.sqrt(np.pi) * erf(zeta[np.where(zeta != 0)]) / zeta[np.where(zeta != 0)] 
    	S_meas[np.where(zeta == 0)] = S_peak

    	W_i_smear_s = S_peak * W_i_s / S_meas #Area of pulse = S*W is conserved 
    	W_b_smear_s = np.sqrt(dt**2 + t_chan(DM,c_freq_MHz,chan_width)**2 + W_i_smear_s**2)

    	S_N_smeared = S_meas * W_i_smear_s / np.sqrt(W_b_smear_s)  
    
    	min_group = min(MIN_GROUP,int(0.8*2*len(S_N_smeared[S_N_smeared>5])-1))
    return min_group
