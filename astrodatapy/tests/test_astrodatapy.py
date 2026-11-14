from astrodatapy.misc_tools import available_redshifts

def test_number_density():
    from astrodatapy.number_density import number_density as function
    features = ['GSMF', 'GSMF_Blue', 'GSMF_Red', 'GSMF_Bulge', 'GSMF_Disk',\
                'GSMF_Quiescent', 'BHM', 'BHMF', 'QLF_bolometric',\
                'QLF_optical', 'QLF_UV', 'QLF_Xray', 'SFRF', 'GLF_UV', 'GLF_OIII']

    for feature in features:
        for z_target in available_redshifts(feature):
            assert function(feature, z_target)

def test_bhm_masses_in_msun():
    # Supermassive black holes; Caramete2010 is tabulated in 10^6 Msun and must be rescaled
    import numpy as np
    from astrodatapy.number_density import number_density
    obs = number_density('BHM', 0.0, quiet=True, h=0.7)
    mass = obs.target_observation.loc['Caramete2010', 'Data'][:, 1]
    assert 1e7 < np.median(mass) < 1e9

CORRELATION_FEATURES = ['Magorrian', 'Tully_Fisher', 'DiskSize_StellarMass',
                        'GasFraction_StellarMass', 'sSFR_StellarMass_Blue',
                        'HaloMass_StellarMass', 'HaloMass_StellarMass_Blue',
                        'HaloMass_StellarMass_Red']

def test_correlation():
    from astrodatapy.correlation import correlation as function
    for feature in CORRELATION_FEATURES:
        for z_target in available_redshifts(feature):
            assert function(feature, z_target)

def test_correlation_info_has_h_and_imf():
    # every correlation dataset says which h and IMF its values assume and how x and y scale with h
    from astrodatapy.correlation import correlation, list_IMF
    for feature in CORRELATION_FEATURES:
        info = correlation(feature, available_redshifts(feature)[0], quiet=True).info
        for name, row in info.iterrows():
            assert 0.4 < float(row['h']) <= 1, (feature, name)
            assert -3 <= float(row['hPowerx']) <= 0, (feature, name)
            assert -3 <= float(row['hPowery']) <= 0, (feature, name)
            assert row['IMF'] in list(list_IMF) + ['None'], (feature, name)

def test_correlation_h_and_imf_conversion():
    # between h=0.5 with Salpeter and h=1 with Chabrier, each axis changes by 2^hPower, and the
    # stellar-mass axis (and inversely the gas fraction) also by the Salpeter -> Chabrier factor
    import numpy as np
    from astrodatapy.correlation import (correlation, Cha2Sal, _stellar_mass_axis,
                                         _inverse_stellar_mass_axis, _axis_columns)
    for feature in CORRELATION_FEATURES:
        for z_target in set(available_redshifts(feature)):
            a = correlation(feature, z_target, quiet=True, z_tol=0.01, h=0.5, IMF_out='Salpeter')
            b = correlation(feature, z_target, quiet=True, z_tol=0.01, h=1.0, IMF_out='Chabrier')
            for name in a.target_observation.index:
                info = a.info.loc[name]
                expected = a.target_observation.loc[name, 'Data'].copy()
                expected[:, 0:3] *= 2.0**float(info['hPowerx'])
                expected[:, 3:6] *= 2.0**float(info['hPowery'])
                if info['IMF'] != 'None':
                    expected[:, _axis_columns[_stellar_mass_axis[feature]]] *= 10**-Cha2Sal
                    if feature in _inverse_stellar_mass_axis:
                        expected[:, _axis_columns[_inverse_stellar_mass_axis[feature]]] /= 10**-Cha2Sal
                got = b.target_observation.loc[name, 'Data']
                assert np.allclose(got, expected, rtol=1e-12, atol=0, equal_nan=True), (feature, name)

def test_clustering():
    from astrodatapy.clustering import clustering as function
    features = ['QC_2PTCF',]
    for feature in features:
        for z_target in available_redshifts(feature):
            assert function(feature, z_target)
