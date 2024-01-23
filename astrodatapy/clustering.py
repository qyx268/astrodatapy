from __future__ import print_function
import os
import glob
import pandas as pd
import numpy  as np
from .misc_tools import *

fileDir = os.path.dirname(os.path.abspath(__file__))
class clustering:
    def __init__(self, feature=None, z_target=None, quiet=False, \
                 h = 1, z_tol=0.25, folder=fileDir+'/data/'):
        
        '''
        feature:  QC_2PTCF
        z_target: target redshift
        z_tol:    tolerance for the target redshift, default = 0.25
        h:        Hubble constant in units of 100 km/s/Mpc, default = 1.0
        quiet:    if True, no warning or processing is give, default = False
        folder:   where observational data is, default = ./data/
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
                         sep='\t\t', skiprows=self.info_skiprows, engine='python')
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
        if 'Type' in info.index:
            if info.Type != 'data':
                if info.Type == 'PowerLaw':
                    data = self._PowerLaw(*data[0])
                if info.Type == 'PowerLaw_2COMPONENTS':
                    data = self._PowerLaw_2COMPONENTS(*data[0])
        elif not self.quiet:
            print("WARNING! No Infomation about Type, dafault is data")
        
        return data

    def _PowerLaw(self, xmin, xmax, r_0, delta_r_0, gamma, delta_gamma):
        bins       = np.linspace(xmin, xmax, (xmax - xmin) * 100 + 1)
        width      = (bins[1] - bins[0]) / 2.
        middles    = bins[:-1] + width
        bins       = 10**bins
        middles    = 10**middles
        r_0       /= self.h
        delta_r_0 /= self.h

        y       = (middles/r_0)**-(gamma)
        delta_y = y * (np.abs(np.log(middles / r_0)) * delta_gamma +\
                       np.abs(gamma / r_0) * delta_r_0)
        data      = np.zeros([len(middles),4])
        data[:,0] = middles
        data[:,1] = y
        data[:,2] = y + delta_y
        data[:,3] = y - delta_y
        return data

    def _PowerLaw_2COMPONENTS(self, xmin, xmax, r_0_GG, delta_r_0_GG, gamma_GG, delta_gamma_GG,\
                                                r_0_QG, delta_r_0_QG, gamma_QG, delta_gamma_QG):
        bins          = np.linspace(xmin, xmax, int((xmax - xmin) * 100) + 1)
        width         = (bins[1] - bins[0]) / 2.
        middles       = bins[:-1] + width
        bins          = 10**bins
        middles       = 10**middles
        r_0_GG       /= self.h
        delta_r_0_GG /= self.h
        r_0_QG       /= self.h
        delta_r_0_QG /= self.h 

        y       = ((middles / r_0_QG)**-gamma_QG)**2. / ((middles / r_0_GG)**-gamma_GG)
        delta_y = y * (np.abs(np.log(middles / r_0_GG)) * delta_gamma_GG +\
                       np.abs(gamma_GG / r_0_GG) * delta_r_0_GG +\
                       np.abs(np.log(middles / r_0_QG)) * delta_gamma_QG +\
                       np.abs(gamma_QG / r_0_QG) * delta_r_0_QG)
        data      = np.zeros([len(middles),4])
        data[:,0] = middles
        data[:,1] = y
        data[:,2] = y + delta_y
        data[:,3] = y - delta_y
        return data
