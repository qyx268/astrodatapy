from __future__ import print_function
import os
import glob
import pandas as pd
import numpy  as np
from .misc_tools import *

# For DataFormat of correlation
List_DataFormat = ['xULLimitsyULLimits', 'xULLimitsyLULimits', 'xULLimitsyULDeltas', 'xULLimitsyLUDeltas', 'xULLimitsyDelta', 'xULLimitsy',\
                   'xLULimitsyULLimits', 'xLULimitsyLULimits', 'xLULimitsyULDeltas', 'xLULimitsyLUDeltas', 'xLULimitsyDelta', 'xLULimitsy',\
                   'xULDeltasyULLimits', 'xULDeltasyLULimits', 'xULDeltasyULDeltas', 'xULDeltasyLUDeltas', 'xULDeltasyDelta', 'xULDeltasy',\
                   'xLUDeltasyULLimits', 'xLUDeltasyLULimits', 'xLUDeltasyULDeltas', 'xLUDeltasyLUDeltas', 'xLUDeltasyDelta', 'xLUDeltasy',\
                   'xDeltayULLimits'   , 'xDeltayLULimits'   , 'xDeltayULDeltas'   , 'xDeltayLUDeltas'   , 'xDeltayDelta'   , 'xDeltay',\
                   'xyULLimits'        , 'xyLULimits'        , 'xyULDeltas'        , 'xyLUDeltas'        , 'xyDelta'        , 'xy']


fileDir = os.path.dirname(os.path.abspath(__file__))
class correlation:
    def __init__(self, feature=None, z_target=None, quiet=False, \
                 h = 1, z_tol=0.25, folder=fileDir+'/data/'):
        
        '''
        feature:  Magorrian, Tully_Fisher, DiskSize_StellarMass,
                  GasFraction_StellarMass, sSFR_StellarMass_Blue,
                  HaloMass_StellarMass, HaloMass_StellarMass_Blue,
                  HaloMass_StellarMass_Red
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
        info = pd.read_csv(self.folder + self.feature + '/info.txt', keep_default_na=False, \
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
        if 'DataFormat' in info.index:
            if info.DataFormat in List_DataFormat:
                data = self._convert_DataFormat(data, info.DataFormat)
            else:
                print('Error! Unrecognized Format of Data (%s)! quit'%info.DataFormat)
                return -2
        elif not self.quiet:
            print("WARNING! No Infomation about DataFormat, dafault is xULLimitsyULLimits")
        
        # convert Log_x(y) to x(y)  
        if 'Log' in info.index:
            data = self._convert_regime(data, info.Log)
        elif not self.quiet:
            print("WARNING! No Infomation about Log, dafault is None, both x and y are in linear regime")
        
        # add normalization of x(y)
        if 'NormalizationOfx' in info.index:
            data = self._convert_x_normalization(data, float(info.NormalizationOfx))
        elif not self.quiet:
            print("WARNING! No Infomation about NormalizationOfx, dafault is 1")
        if 'NormalizationOfy' in info.index:
            data = self._convert_y_normalization(data, float(info.NormalizationOfy))
        elif not self.quiet:
            print("WARNING! No Infomation about NormalizationOfy, dafault is 1")

        return data
        
    def _convert_DataFormat(self, data, DataFormat):
        if (DataFormat != 'xULLimitsyULLimits') and (not self.quiet):
            print("Converting the DataFormat from %s to xULLimitsyULLimits"%DataFormat)
        data_new = np.zeros([len(data),6])
        if DataFormat == 'xULLimitsyULLimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,1]
            data_new[:,2] = data[:,2]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,4]
            data_new[:,5] = data[:,5]
        if DataFormat == 'xULLimitsyLULimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,1]
            data_new[:,2] = data[:,2]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,5]
            data_new[:,5] = data[:,4]
        if DataFormat == 'xULLimitsyULDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,1]
            data_new[:,2] = data[:,2]
            data_new[:,3] = data[:,3]
            if data[0,5]  > 0:
                data[:,5] = -data[:,5]
            data_new[:,4] = data[:,3] + data[:,4]
            data_new[:,5] = data[:,3] + data[:,5]
        if DataFormat == 'xULLimitsyLUDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,1]
            data_new[:,2] = data[:,2]
            data_new[:,3] = data[:,3]
            if data[0,4]  > 0:
                data[:,4] = -data[:,4]
            data_new[:,4] = data[:,3] + data[:,5]
            data_new[:,5] = data[:,3] + data[:,4]
        if DataFormat == 'xULLimitsyDelta':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,1]
            data_new[:,2] = data[:,2]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,3] + data[:,4]
            data_new[:,5] = data[:,3] - data[:,4]
        if DataFormat == 'xULLimitsy':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,1]
            data_new[:,2] = data[:,2]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,3]
            data_new[:,5] = data[:,3]

        if DataFormat == 'xLULimitsyULLimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,2]
            data_new[:,2] = data[:,1]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,4]
            data_new[:,5] = data[:,5]
        if DataFormat == 'xLULimitsyLULimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,2]
            data_new[:,2] = data[:,1]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,5]
            data_new[:,5] = data[:,4]
        if DataFormat == 'xLULimitsyULDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,2]
            data_new[:,2] = data[:,1]
            data_new[:,3] = data[:,3]
            if data[0,5]  > 0:
                data[:,5] = -data[:,5]
            data_new[:,4] = data[:,3] + data[:,4]
            data_new[:,5] = data[:,3] + data[:,5]
        if DataFormat == 'xLULimitsyLUDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,2]
            data_new[:,2] = data[:,1]
            data_new[:,3] = data[:,3]
            if data[0,4]  > 0:
                data[:,4] = -data[:,4]
            data_new[:,4] = data[:,3] + data[:,5]
            data_new[:,5] = data[:,3] + data[:,4]
        if DataFormat == 'xLULimitsyDelta':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,2]
            data_new[:,2] = data[:,1]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,3] + data[:,4]
            data_new[:,5] = data[:,3] - data[:,4]
        if DataFormat == 'xLULimitsy':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,2]
            data_new[:,2] = data[:,1]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,3]
            data_new[:,5] = data[:,3]

        if DataFormat == 'xULDeltasyULLimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            if data[0,2] > 0:
                data[:,2] = - data[:,2]
            data_new[:,2] = data[:,0] + data[:,2]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,4]
            data_new[:,5] = data[:,5]
        if DataFormat == 'xULDeltasyLULimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            if data[0,2] > 0:
                data[:,2] = - data[:,2]
            data_new[:,2] = data[:,0] + data[:,2]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,5]
            data_new[:,5] = data[:,4]
        if DataFormat == 'xULDeltasyULDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            if data[0,2] > 0:
                data[:,2] = - data[:,2]
            data_new[:,2] = data[:,0] + data[:,2]
            data_new[:,3] = data[:,3]
            if data[0,5]  > 0:
                data[:,5] = -data[:,5]
            data_new[:,4] = data[:,3] + data[:,4]
            data_new[:,5] = data[:,3] + data[:,5]
        if DataFormat == 'xULDeltasyLUDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            if data[0,2] > 0:
                data[:,2] = - data[:,2]
            data_new[:,2] = data[:,0] + data[:,2]
            data_new[:,3] = data[:,3]
            if data[0,4]  > 0:
                data[:,4] = -data[:,4]
            data_new[:,4] = data[:,3] + data[:,5]
            data_new[:,5] = data[:,3] + data[:,4]
        if DataFormat == 'xULDeltasyDelta':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            if data[0,2] > 0:
                data[:,2] = - data[:,2]
            data_new[:,2] = data[:,0] + data[:,2]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,3] + data[:,4]
            data_new[:,5] = data[:,3] - data[:,4]
        if DataFormat == 'xULDeltasy':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            if data[0,2] > 0:
                data[:,2] = - data[:,2]
            data_new[:,2] = data[:,0] + data[:,2]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,3]
            data_new[:,5] = data[:,3]

        if DataFormat == 'xLUDeltasyULLimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,2]
            if data[0,1] > 0:
                data[:,1] = - data[:,1]
            data_new[:,2] = data[:,0] + data[:,1]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,4]
            data_new[:,5] = data[:,5]
        if DataFormat == 'xLUDeltasyLULimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,2]
            if data[0,1] > 0:
                data[:,1] = - data[:,1]
            data_new[:,2] = data[:,0] + data[:,1]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,5]
            data_new[:,5] = data[:,4]
        if DataFormat == 'xLUDeltasyULDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,2]
            if data[0,1] > 0:
                data[:,1] = - data[:,1]
            data_new[:,2] = data[:,0] + data[:,1]
            data_new[:,3] = data[:,3]
            if data[0,5]  > 0:
                data[:,5] = -data[:,5]
            data_new[:,4] = data[:,3] + data[:,4]
            data_new[:,5] = data[:,3] + data[:,5]
        if DataFormat == 'xLUDeltasyLUDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,2]
            if data[0,1] > 0:
                data[:,1] = - data[:,1]
            data_new[:,2] = data[:,0] + data[:,1]
            data_new[:,3] = data[:,3]
            if data[0,4]  > 0:
                data[:,4] = -data[:,4]
            data_new[:,4] = data[:,3] + data[:,5]
            data_new[:,5] = data[:,3] + data[:,4]
        if DataFormat == 'xLUDeltasyDelta':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,2]
            if data[0,1] > 0:
                data[:,1] = - data[:,1]
            data_new[:,2] = data[:,0] + data[:,1]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,3] + data[:,4]
            data_new[:,5] = data[:,3] - data[:,4]
        if DataFormat == 'xLUDeltasy':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,2]
            if data[0,1] > 0:
                data[:,1] = - data[:,1]
            data_new[:,2] = data[:,0] + data[:,1]
            data_new[:,3] = data[:,3]
            data_new[:,4] = data[:,3]
            data_new[:,5] = data[:,3]

        if DataFormat == 'xDeltayULLimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            data_new[:,2] = data[:,0] - data[:,1]
            data_new[:,3] = data[:,2]
            data_new[:,4] = data[:,3]
            data_new[:,5] = data[:,4]
        if DataFormat == 'xDeltayLULimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            data_new[:,2] = data[:,0] - data[:,1]
            data_new[:,3] = data[:,2]
            data_new[:,4] = data[:,4]
            data_new[:,5] = data[:,3]
        if DataFormat == 'xDeltayULDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            data_new[:,2] = data[:,0] - data[:,1]
            data_new[:,3] = data[:,2]
            if data[0,4]  > 0:
                data[:,4] = -data[:,4]
            data_new[:,4] = data[:,2] + data[:,3]
            data_new[:,5] = data[:,2] + data[:,4]
        if DataFormat == 'xDeltayLUDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            data_new[:,2] = data[:,0] - data[:,1]
            data_new[:,3] = data[:,2]
            if data[0,3]  > 0:
                data[:,3] = -data[:,3]
            data_new[:,4] = data[:,2] + data[:,4]
            data_new[:,5] = data[:,2] + data[:,3]
        if DataFormat == 'xDeltayDelta':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            data_new[:,2] = data[:,0] - data[:,1]
            data_new[:,3] = data[:,2]
            data_new[:,4] = data[:,2] + data[:,3]
            data_new[:,5] = data[:,2] - data[:,3]
        if DataFormat == 'xDeltay':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0] + data[:,1]
            data_new[:,2] = data[:,0] - data[:,1]
            data_new[:,3] = data[:,2]
            data_new[:,4] = data[:,2]
            data_new[:,5] = data[:,2]

        if DataFormat == 'xyULLimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0]
            data_new[:,2] = data[:,0]
            data_new[:,3] = data[:,1]
            data_new[:,4] = data[:,2]
            data_new[:,5] = data[:,3]
        if DataFormat == 'xyLULimits':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0]
            data_new[:,2] = data[:,0]
            data_new[:,3] = data[:,1]
            data_new[:,4] = data[:,3]
            data_new[:,5] = data[:,2]
        if DataFormat == 'xyULDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0]
            data_new[:,2] = data[:,0]
            data_new[:,3] = data[:,1]
            if data[0,3]  > 0:
                data[:,3] = -data[:,3]
            data_new[:,4] = data[:,1] + data[:,2]
            data_new[:,5] = data[:,1] + data[:,3]
        if DataFormat == 'xyLUDeltas':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0]
            data_new[:,2] = data[:,0]
            data_new[:,3] = data[:,1]
            if data[0,2]  > 0:
                data[:,2] = -data[:,2]
            data_new[:,4] = data[:,1] + data[:,3]
            data_new[:,5] = data[:,1] + data[:,2]
        if DataFormat == 'xyDelta':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0]
            data_new[:,2] = data[:,0]
            data_new[:,3] = data[:,1]
            data_new[:,4] = data[:,1] + data[:,2]
            data_new[:,5] = data[:,1] - data[:,2]
        if DataFormat == 'xy':
            data_new[:,0] = data[:,0]
            data_new[:,1] = data[:,0]
            data_new[:,2] = data[:,0]
            data_new[:,3] = data[:,1]
            data_new[:,4] = data[:,1]
            data_new[:,5] = data[:,1]

        return data_new

    def _convert_regime(self, data, Log):
        if Log == 'x':
            if not self.quiet:
                print("Converting x from logarithm to linear")
            data[:,0:3] = 10**data[:,0:3]
        elif Log == 'y':
            if not self.quiet:
                print("Converting y from logarithm to linear")
            data[:,3:6] = 10**data[:,3:6]
        elif Log == 'xy' or Log == 'yx':
            if not self.quiet:
                print("Converting x and y from logarithm to linear")
            data = 10**data
        return data

    def _convert_x_normalization(self, data, NormalizationOfx):
        if NormalizationOfx!=1:
            if not self.quiet:
                print("Converting the Normalization of y from %.2e to 1"%NormalizationOfx)
            data[:,0:3] *= NormalizationOfx
        return data

    def _convert_y_normalization(self, data, NormalizationOfy):
        if NormalizationOfy!=1:
            if not self.quiet:
                print("Converting the Normalization of y from %.2e to 1"%NormalizationOfy)
            data[:,3:6] *= NormalizationOfy
        return data
