from __future__ import print_function
import os
import glob
import pandas as pd
import numpy  as np
from .misc_tools import *

fileDir = os.path.dirname(os.path.abspath(__file__))
# For IMF convertion
Cha2Sal     = np.log10(1.8)
Kro2Sal     = np.log10(1.5)
DietSal2Sal = -np.log10(0.7)
list_IMF    = ['Salpeter', 'DietSalpeter', 'Kroupa', 'Chabrier']

# For Magnitude convertion
alpha_q      = 1.57     # quasar spectral index in UV
alpha_q_op   = 0.44  # quasar spectral index in optical
_1450_to_912 = (1200./1450)**alpha_q_op*(912/1200.)**alpha_q
_Mi_to_Mg    = -2.5*alpha_q_op*np.log10(4670./7471)
_MB_to_Mg    = -2.5*(1.+alpha_q_op)*np.log10(3)-0.14+0.209
_Mg_to_MB    = -_MB_to_Mg
_MABB_to_MB  = 0.09
_MB_to_MABB  = -_MABB_to_MB
list_MAG     = ['MB2MABB', 'M1450', 'LBsol2M1450', 'MABB', 'Mg2MB', 'M14502MB', 'M1500', 'M1600', 'Bolometric']

# For y Error Format of number density
list_Errors  = ['ULLimits', 'LULimits', 'ULDeltas', 'LUDeltas', 'Delta', 'None', np.NaN]

class number_density:
    
    def __init__(self, feature=None, z_target=None, quiet=False, \
                 h = 1, z_tol=0.25, folder=fileDir+'/data/', IMF_out="Kroupa"):
        
        '''
        feature:  GSMF, GSMF_Blue, GSMF_Red, GSMF_Bulge, GSMF_Disk, GSMF_Quiescent
                  BHM, BHMF, QLF_bolometric, QLF_optical, QLF_UV
                  SFRF, GLF_UV
        z_target: target redshift
        z_tol:    tolerance for the target redshift, default = 0.25
        h:        Hubble constant in units of 100 km/s/Mpc, default = 1.0
        quiet:    if True, no warning or processing is give, default = False
        folder:   where observational data is, default = ./data/
        IMF_out:  the IMF to use to convert data, default = Kroupa
        '''
        self.feature                  = feature
        self.folder                   = folder+'/'
        self.z_target                 = z_target
        self.z_tol                    = z_tol
        self.h                        = h
        self.quiet                    = quiet
        self.info_header              = np.zeros(20, dtype=object)
        self.info_skiprows            = None
        self.info                     = None
        self.n_available_observation  = 0
        self.z_available_observation  = available_redshifts(feature)
        self.available_observation    = None
        self.n_target_observation     = 0
        self.target_observation       = None
        self.IMF_out                  = IMF_out
        if self.feature is not None and self.z_target is not None:
            if not self.quiet:
                print("You are requesting %s at z_target=%.2f with a tolerance of z_tol=%.2f and h=%.3f"%(self.feature, self.z_target, self.z_tol,self.h))
                print("quiet=True to silent")
            with open(self.folder + self.feature + '/info.txt','r') as f:
                for num, line in enumerate(f):
                    self.info_header[num] = line
                    if '#Name' in line:
                        self.info_skiprows = num
                        break
            self._load_info()
        else:
            print("I need a feature and a redshift!")
            return
    
    def _load_info(self):
        info = pd.read_csv(self.folder + self.feature + '/info.txt', keep_default_na=False,\
                         delim_whitespace=True, skiprows=self.info_skiprows, engine='python')
        info.set_index('#Name', inplace=True)
        self.n_available_observation = len(info.index.values)
        if self.n_available_observation == 0:
            print("No available data of %s"%(self.feature))
            return
        else:
            self.available_observation = info.index.values
            self.info                    = info
            if not self.quiet:
                print("available data of %s includes:"%self.feature)
                print(*self.available_observation)
            self._target_observations()
            
    def _target_observations(self):
        flags = np.zeros(self.n_available_observation, dtype=bool)
        names = np.zeros(self.n_available_observation, dtype=object)
        for ii, observation in enumerate(self.available_observation):
            list_redshift = glob.glob1(self.folder + self.feature + '/' +\
                                       observation,'z*.dat')
            for fname in list_redshift:
                if np.abs(filename_to_redshift(fname) - self.z_target) < self.z_tol:
                    flags[ii] = True
                    names[ii] = fname
        
        self.n_target_observation = sum(flags)
        
        if self.n_target_observation == 0:
            print("No available data of %s at %.2f<z_target<%.2f"%(self.feature, \
                                self.z_target - self.z_tol, self.z_target + self.z_tol))
        else:
            self.target_observation = pd.DataFrame({'Name': self.available_observation[flags],
                                                    'DataType': np.zeros(self.n_target_observation, dtype=object),
                                                    'FileName': names[flags]})
            self._load_observational_data()
        
    def _load_observational_data(self):
        # load and convert observation to Salpeter IMF,
        # phi in linear and errors in upper, lower limits
             
        results = [[],] * self.n_target_observation
        for ii in range(self.n_target_observation):
            if not self.quiet:
                print("\nLoading observational data from %s..."%self.target_observation['Name'][ii])
                print("Filename %s"%(self.folder + self.feature + '/' +self.target_observation['FileName'][ii]))
            
            info = self.info.loc[self.target_observation['Name'][ii]]
            data = np.loadtxt(self.folder + self.feature + '/' +\
                              self.target_observation['Name'][ii] + '/' +\
                              self.target_observation['FileName'][ii])
            if data.ndim == 1:
                data = data.reshape([1,-1])
            
            if 'Type' in info.index:
                self.target_observation['DataType'][ii] = info.Type
            else:
                self.target_observation['DataType'][ii] = 'Unknown'
            results[ii] = self._convert_observational_data(data, info)
            if not self.quiet:
                print("..done")
        self.target_observation = self.target_observation.assign(Data=pd.Series(results).values)
        self.target_observation.set_index('Name', inplace=True)

    def _convert_observational_data(self, data, info):
        # check if this is data or function
        if 'Type' in info.index:
            if info.Type != 'data':
                if info.Type == 'Schechter':
                    if self.feature in ['GSMF', 'GSMF_Blue', 'GSMF_Red', 'GSMF_Bulge', 'GSMF_Disk', 'GSMF_Quiescent']:
                        data = self._Schechter(*data[0])
                    if self.feature in ['QLF_UV','QLF_optical','GLF_UV']:
                        data = self._SchechterMagnitude(*data[0])
                if info.Type == 'DoubleSchechter':
                    if self.feature in ['QLF_UV','QLF_optical','GLF_UV']:
                        data = self._DoubleSchechterMagnitude(*data[0])
                if info.Type == 'Hopkins2007':
                    if self.feature in ['QLF_bolometric',]:
                        data = self._Hopkins2007(*data[0])
        elif not self.quiet:
            print("WARNING! No Infomation about Type, default is data")
            
        # convert to Salpeter IMF
        if 'IMF' in info.index:
            if info.IMF in list_IMF:
                data = self._convert_IMF(data, info.IMF, self.IMF_out)
            else:
                print('Error! Unrecognized IMF (%s)! quit'%info.IMF)
                return -1
        elif not self.quiet:
            print("WARNING! No Infomation about IMF, dafault is Salpeter")
        
        # convert magnitude
        if 'MAG' in info.index:
            if info.MAG in list_MAG:
                data = self._convert_MAG(data, info.MAG)
            else:
                print('Error! Unrecognized MAG (%s)! quit'%info.MAG)
                return -1
        elif not self.quiet:
            print("WARNING! No Infomation about MAG, dafault is Salpeter")

        # convert to format of LogM, phi, phi_u_limit, phi_l_limit
        if 'Errors' in info.index:
            if info.Errors in list_Errors:
                data = self._convert_Errors(data, info.Errors)
            else:
                print('Error! Unrecognized Format of Errors (%s)! quit'%info.Errors)
                return -2
        elif not self.quiet:
            print("WARNING! No Infomation about Errors, dafault is ULLimits")
        
        # convert Log_phi to phi  
        if 'Log' in info.index:
            data = self._convert_phi_regime(data, int(info.Log))
        elif not self.quiet:
            print("WARNING! No Infomation about Log, dafault is 0 (linear)")
        
        # add normalization of phi
        if 'NormalizationOfPhi' in info.index:
            data = self._convert_phi_normalization(data, float(info.NormalizationOfPhi))
        elif not self.quiet:
            print("WARNING! No Infomation about NormalizationOfPhi, dafault is 1")
        
        # convert M to Msol with input h
        data = self._convert_h_x(data, info.h)
        
        # convert phi to Mpc**-3dex**-1 with input h
        if 'HinPhi' in info.index:
            data = self._convert_h_y(data, info.h, int(info.HinPhi))
        else:
            data[:,1:4] *= self.h**3
            if not self.quiet:
                print("WARNING! No Infomation about HinPhi, dafault is 1")
                print("Converting y axis from h=1.000 to h=%.3f with a power of 3"%(self.h))
        
        return data
    
    def _convert_IMF(self, data, IMF, IMF_out):
        if IMF == IMF_out:
            return data
        if IMF != 'Salpeter':
            if not self.quiet:
                print("Converting IMF from %s to Salpeter"%IMF)
            if IMF == 'DietSalpeter':
                data[:,0] += DietSal2Sal
            if IMF == 'Kroupa':
                data[:,0] += Kro2Sal
            if IMF == 'Chabrier':
                data[:,0] += Cha2Sal
        if IMF_out != 'Salpeter':
            if not self.quiet:
                print("Converting IMF from Salpeter to %s"%IMF_out)
            if IMF_out == 'DietSalpeter':
                data[:,0] -= DietSal2Sal
            if IMF_out == 'Kroupa':
                data[:,0] -= Kro2Sal
            if IMF_out == 'Chabrier':
                data[:,0] -= Cha2Sal
        return data
        
    def _convert_MAG(self, data, MAG):
        if MAG not in ['M1450', 'MABB']:
            if not self.quiet:
                print("Converting Magnituedes according to %s"%MAG)
            if MAG == 'MB2MABB':
                data[:,0] += _MB_to_MABB
            if MAG == 'LBsol2M1450':
                from magnitudes import _MBvega2M1450
                data[:,1] *= 0.4*data[:,0]*np.log(10)
                data[:,0]  = _MBvega2M1450(5.48-2.5*np.log10(data[:,0]))
            if MAG == 'Mg2MB':
                data[:,0] += _Mg_to_MB
            if MAG == 'M14502MB':
                from magnitudes import _M14502MB
                data[:,0] = _M14502MB(data[:,0])
        return data

    def _convert_Errors(self, data, Errors):
        if Errors == 'None' and data.shape[1]==2:
            data = np.concatenate((data, data[:,-1][:,None], data[:,-1][:,None]), axis=1)
            #data = np.concatenate((data, -data[:,-1][:,None], -data[:,-1][:,None]*2), axis=1) why did i do this?
        if Errors != 'ULLimits' and Errors != 'None':
            if not self.quiet:
                print("Converting the Errors from %s to Upper Lower Limits"%Errors)
            if Errors == 'Delta':
                data = np.concatenate((data, -data[:,-1][:,None]), axis=1)
            if Errors == 'LULimits' or Errors == 'LUDeltas':
                data[:,2:4] = data[:,3:1:-1]
            if Errors == 'ULDeltas' or Errors == 'LUDeltas' or Errors == 'Delta':
                if np.sum(data[:,3])  > 0:
                    data[:,3] = -data[:,3]
                data[:,2:4] += data[:,1][:,None]
        return data            
        
    def _convert_phi_regime(self, data, Log):
        if Log == 1:
            if not self.quiet:
                print("Converting Phi from logarithm to linear")
            data[:,1:4] = 10**data[:,1:4]
        return data
    
    def _convert_phi_normalization(self, data, NormalizationOfPhi):
        if NormalizationOfPhi!=1:
            if not self.quiet:
                print("Converting the Normalization of Phi from %.2e to 1"%NormalizationOfPhi)
            data[:,1:4] *= NormalizationOfPhi
        return data
        
    def _convert_h_x(self, data, h_obs):
        if self.feature in ['GSMF', 'GSMF_Blue', 'GSMF_Red', 'GSMF_Bulge', 'GSMF_Disk', 'QLF_bolometric', 'GSMF_Quiescent']:
            if not self.quiet:
                print("Converting x axis from h=%.3f to h=%.3f with a power of -2 because it is %s"%(h_obs,self.h, self.feature))
            data[:,0] -= 2 * np.log10(self.h/h_obs)
        elif self.feature in ['QLF_optical','QLF_UV', 'GLF_UV']:
            if not self.quiet:
                print("Converting x axis from h=%.3f to h=%.3f with a power of -5 because it is %s"%(h_obs,self.h, self.feature))
            data[:,0] -= 5 * np.log10(self.h/h_obs)
        elif self.feature in ['SFRF','BHMF']:
            if not self.quiet:
                print("Converting x axis from h=%.3f to h=%.3f with a power of 0 because it is %s"%(h_obs,self.h, self.feature))
        else:
            if not self.quiet:
                print("Warning! Don't know how to convert h for x axis of %s"%(self.feature))
        return data

    def _convert_h_y(self, data, h_obs, HinPhi):
        if self.feature in ['GSMF', 'GSMF_Blue', 'GSMF_Red', 'GSMF_Bulge', 'GSMF_Disk', 'GSMF_Quiescent',\
                            'BHMF', 'QLF_bolometric', 'QLF_optical', 'QLF_UV','SFRF', 'GLF_UV']:
            if HinPhi == 1:
                if not self.quiet:
                    print("Converting y axis from h=1.000 to h=%.3f with a power of 3"%(self.h))
                data[:,1:4] *= self.h**3
            else:
                if not self.quiet:
                    print("Converting y axis from h=%.3f to h=%.3f with a power of 3"%(h_obs,self.h))
                data[:,1:4] *= (self.h/h_obs)**3
        else:
            if not self.quiet:
                print("Warning! Don't know how to convert h for y axis of %s"%(self.feature))
        return data
    
    def _Schechter(self, xmin, xmax, logM_star, alpha1, phi_star1, alpha2=0.0, phi_star2=0.0):
        if not self.quiet:
            print("Using Schechter Function for Mass")
        data      = np.zeros([100,4])
        data[:,0] = np.linspace(xmin,xmax,100)
        data[:,1] = 2.302585093 * np.exp(-10 ** (data[:,0] - logM_star)) *\
                    (phi_star1 * (10 ** ((data[:,0] - logM_star)*(1. + alpha1))) +\
                     phi_star2 * (10 ** ((data[:,0] - logM_star)*(1. + alpha2))))
        data[:,2] = data[:,1]
        data[:,3] = data[:,1]
        return data

    def _SchechterMagnitude(self, xmin, xmax, mag_star, phi_star, alpha):
        if not self.quiet:
            print("Using Schechter Function for Magnitude")
        data      = np.zeros([100,4])
        data[:,0] = np.linspace(xmin,xmax,100)
        MStarMinM = 0.4 * (mag_star - data[:,0])
        data[:,1] = 2.302585093 * 0.4 * np.exp(-10.**(0.4 * (mag_star - data[:,0]))) *\
                    phi_star / (10 ** (-0.4 * (mag_star - data[:,0]) * (alpha + 1.)))
                   
        data[:,2] = data[:,1]
        data[:,3] = data[:,1]
        return data

    def _DoubleSchechterMagnitude(self, xmin, xmax, mag_star, phi_star, alpha1, alpha2):
        if not self.quiet:
            print("Using Double Schechter Function for Magnitude")
        data      = np.zeros([100,4])
        data[:,0] = np.linspace(xmin,xmax,100)
        data[:,1] = phi_star / (10 ** (-0.4 * (mag_star - data[:,0]) * (alpha1 + 1.)) +\
                                10 ** (-0.4 * (mag_star - data[:,0]) * (alpha2 + 1.)))

        data[:,2] = data[:,1]
        data[:,3] = data[:,1]
        return data

    def _Hopkins2007(self, log_x_min, log_x_max, log_phi_star, delta_log_phi_star, log_L_star,\
                     delta_log_L_star, gamma_1, delta_gamma_1, gamma_2, delta_gamma_2):
        if not self.quiet:
            print("Using Hopkins 2007 Function for Bolometric Luminosity")
        data      = np.zeros([100,4])
        data[:,0] = np.linspace(log_x_min,log_x_max,100)
        x         = 10 ** data[:,0]
        L_star    = 10.**log_L_star
        data[:,1] = log_phi_star - np.log10((x/L_star)**gamma_1 + (x/L_star)**gamma_2)
        delta_log_phi = delta_log_phi_star + 1. / 2.302585093 / ((x / L_star)**gamma_1 + (x / L_star)**gamma_2) *\
                        (((np.abs(np.log(x / L_star))) * delta_gamma_1 + gamma_1 * delta_log_L_star * 2.302585093) * (x / L_star)**gamma_1 +\
                         ((np.abs(np.log(x / L_star))) * delta_gamma_2 + gamma_2 * delta_log_L_star * 2.302585093) * (x / L_star)**gamma_2)
        data[:,2] = data[:,1] + delta_log_phi
        data[:,3] = data[:,1] - delta_log_phi
        return data
