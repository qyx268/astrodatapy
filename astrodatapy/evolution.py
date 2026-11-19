"""Quantities as a function of redshift, such as the ionizing emissivity, the H I photoionization
rate, the mean free path of ionizing photons and the IGM temperature."""
from __future__ import print_function
import os
import pandas as pd
import numpy  as np
from .correlation import correlation as _correlation, List_DataFormat as _List_DataFormat

fileDir = os.path.dirname(os.path.abspath(__file__))

# Units in which evolution returns each feature (the data files keep the published units)
evolution_units = {'IGM_Emissivity': 's^-1 cMpc^-3 (ionizing photons per second per comoving Mpc^3)',
                   'IGM_Gamma': '1e-12 s^-1 (H I photoionization rate in units of 10^-12 s^-1)',
                   'IGM_MFP': 'proper Mpc (mean free path of ionizing photons at 912 A)',
                   'IGM_T0': 'K (IGM temperature at the mean density)'}


class evolution(_correlation):
    """Load every dataset of a quantity that is tabulated as a function of redshift.

    Each dataset of ``feature`` is one file of redshifts and values (with the redshift bin and
    uncertainties when the paper gives them). The values are converted to the units of
    ``evolution_units`` and to the Hubble parameter ``h``: the ionizing emissivity in
    s^-1 cMpc^-3, the H I photoionization rate in **1e-12 s^-1**, the mean free path in
    **proper** Mpc and the temperature in K.

    Parameters
    ----------
    feature : str
        IGM_Emissivity, IGM_Gamma, IGM_MFP or IGM_T0.
    quiet : bool, optional
        If True, do not print what is loaded and converted. Default False.
    h : float, optional
        Hubble parameter, H0 / (100 km/s/Mpc), to convert the data to. Default 1.
    z_min, z_max : float, optional
        Keep only the rows with z_min <= z <= z_max. Default: all redshifts.
    folder : str, optional
        Directory holding the data. Default: the data shipped with the package.

    Attributes
    ----------
    units : str
        Units of the returned values.
    available_observation : array of str
        Names of all datasets of ``feature``.
    n_target_observation : int
        Number of datasets with rows in the redshift range.
    target_observation : pandas.DataFrame
        One row per dataset, indexed by dataset name, with columns ``DataType`` (``data``, or
        ``data_sim`` for a simulation/model curve), ``Method``, ``FileName``, ``Data``,
        ``ULimit`` and ``LLimit``. ``Data`` is an (N, 6) array of z, upper and lower edge of the
        redshift bin, the value, and its upper and lower bound. ``ULimit`` and ``LLimit`` are
        boolean arrays of length N marking upper and lower limits; in those rows the value and
        both bounds hold the limit.
    """

    def __init__(self, feature=None, quiet=False, h=1, z_min=None, z_max=None, folder=fileDir+'/data/'):
        self.feature                  = feature
        self.folder                   = folder+'/'
        self.h                        = h
        self.quiet                    = quiet
        self.z_min                    = z_min
        self.z_max                    = z_max
        self.IMF_out                  = None
        self.units                    = evolution_units.get(feature, '')
        self.info                     = None
        self.n_available_observation  = 0
        self.available_observation    = None
        self.n_target_observation     = 0
        self.target_observation       = None
        if self.feature is None:
            print("I need a feature!")
            return
        if not self.quiet:
            print("You are requesting %s (in %s) with h=%.3f"%(self.feature, self.units, self.h))
            print("quiet=True to silent")
        self._load_info()

    def _load_info(self):
        path = self.folder + self.feature + '/info.txt'
        with open(path, 'r') as f:
            skiprows = next(num for num, line in enumerate(f) if line.startswith('#Name'))
        info = pd.read_csv(path, keep_default_na=False, sep=r'\s+', skiprows=skiprows, engine='python')
        self.info = info.set_index('#Name')
        self.available_observation = self.info.index.values
        self.n_available_observation = len(self.available_observation)
        if not self.quiet:
            print("available data of %s includes:"%self.feature)
            print(*self.available_observation)
        self._load_observational_data()

    def _load_observational_data(self):
        rows = []
        for name in self.available_observation:
            info = self.info.loc[name]
            fname = name + '.dat'
            if not self.quiet:
                print("\nLoading observational data from %s..."%name)
                print("Filename %s"%(self.folder + self.feature + '/' + fname))
            data = np.loadtxt(self.folder + self.feature + '/' + fname, ndmin=2)
            data, ulimit, llimit = self._read_limits(data, info.DataFormat)
            data = self._convert_observational_data(data, info)
            if isinstance(data, int):
                continue
            keep = np.ones(len(data), dtype=bool)
            if self.z_min is not None:
                keep &= data[:,0] >= self.z_min
            if self.z_max is not None:
                keep &= data[:,0] <= self.z_max
            if not keep.any():
                continue
            rows.append((name, info.Type if 'Type' in info.index else 'data',
                         info.Method if 'Method' in info.index else '',
                         fname, data[keep], ulimit[keep], llimit[keep]))
            if not self.quiet:
                print("..done")
        self.n_target_observation = len(rows)
        if not rows:
            print("No available data of %s in the requested redshift range"%self.feature)
            return
        self.target_observation = pd.DataFrame({
            'Name': [r[0] for r in rows], 'DataType': [r[1] for r in rows], 'Method': [r[2] for r in rows],
            'FileName': [r[3] for r in rows], 'Data': [r[4] for r in rows],
            'ULimit': [r[5] for r in rows], 'LLimit': [r[6] for r in rows]}).set_index('Name')

    def _read_limits(self, data, DataFormat):
        # Upper limits are rows whose y error columns are all nan, lower limits rows whose y error
        # columns are all inf (as in the data files). Their errors are set to zero (deltas) or to
        # the value (limits) before the conversion, so that the value and both bounds are the limit.
        yformat = DataFormat[DataFormat.index('y') + 1:] if 'y' in DataFormat else ''
        n_yerr = {'': 0, 'Delta': 1}.get(yformat, 2)
        ulimit = np.zeros(len(data), dtype=bool)
        llimit = np.zeros(len(data), dtype=bool)
        if n_yerr:
            yerr = data[:,data.shape[1]-n_yerr:]
            ulimit = np.isnan(yerr).all(axis=1)
            llimit = np.isinf(yerr).all(axis=1)
            flagged = ulimit | llimit
            if flagged.any():
                fill = 0.0 if 'Delta' in yformat else data[flagged, data.shape[1]-n_yerr-1][:,None]
                data[np.ix_(flagged, np.arange(data.shape[1]-n_yerr, data.shape[1]))] = fill
        return data, ulimit, llimit

    def _convert_observational_data(self, data, info):
        if info.DataFormat not in _List_DataFormat:
            print('Error! Unrecognized Format of Data (%s)! quit'%info.DataFormat)
            return -2
        data = self._convert_DataFormat(data, info.DataFormat)
        if 'Log' in info.index:
            data = self._convert_regime(data, info.Log)
        if 'NormalizationOfy' in info.index:
            data = self._convert_y_normalization(data, float(info.NormalizationOfy))
        # lengths tabulated in comoving units are converted to proper ones
        if 'Comoving' in info.index and str(info.Comoving) == '1':
            if not self.quiet:
                print("Converting y from comoving to proper units, / (1+z)")
            data[:,3:6] /= (1. + data[:,0:1])
        # y scales as h^hPowery (e.g. -1 for lengths tabulated in h^-1 Mpc, with h = 1 in info.txt)
        if 'h' in info.index and info.h != 'None':
            data = self._convert_h(data, float(info.h), 0., float(info.hPowery))
        return data
