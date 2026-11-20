"""Tests of astrodatapy.evolution and the IGM features (quantities tabulated against redshift)."""
import numpy as np

IGM_FEATURES = ['IGM_Emissivity', 'IGM_Gamma', 'IGM_MFP', 'IGM_T0', 'IGM_xHI']
# value = the MAP model, bounds = 16th/84th percentiles: the MAP may lie outside them
MAP_WITH_PERCENTILES = {'Qin2025'}


def _obs(feature, **kw):
    from astrodatapy.evolution import evolution
    return evolution(feature, quiet=True, **kw).target_observation


def test_every_igm_dataset_loads():
    # every dataset loads as (N, 6) rows of z, z_hi, z_lo, y, y_hi, y_lo with consistent bounds,
    # and limits have the value in both bounds
    from astrodatapy.evolution import evolution
    for feature in IGM_FEATURES:
        o = evolution(feature, quiet=True)
        assert o.n_target_observation == o.n_available_observation > 0, feature
        for name, row in o.target_observation.iterrows():
            d, ul, ll = row['Data'], row['ULimit'], row['LLimit']
            assert d.ndim == 2 and d.shape[1] == 6 and np.isfinite(d).all(), (feature, name)
            assert ul.dtype == bool and ll.dtype == bool and len(ul) == len(ll) == len(d), (feature, name)
            assert not (ul & ll).any(), (feature, name)
            lim = ul | ll
            assert np.array_equal(d[lim, 4], d[lim, 3]) and np.array_equal(d[lim, 5], d[lim, 3]), (feature, name)
            if name in MAP_WITH_PERCENTILES:
                assert (d[:, 4] >= d[:, 5]).all(), (feature, name)
            else:
                assert (d[:, 4] >= d[:, 3]).all() and (d[:, 5] <= d[:, 3]).all(), (feature, name)
            assert (d[:, 1] >= d[:, 0]).all() and (d[:, 2] <= d[:, 0]).all(), (feature, name)
            assert (d[:, 3] > 0).all(), (feature, name)
            assert row['DataType'] in ('data', 'data_sim'), (feature, name)


def test_igm_units():
    # returned units: Ndot in s^-1 cMpc^-3, Gamma in 1e-12 s^-1, MFP in proper Mpc, T0 in K, x_HI a fraction
    ranges = {'IGM_Emissivity': (1e48, 1e53), 'IGM_Gamma': (1e-3, 10.), 'IGM_MFP': (1e-2, 1e3), 'IGM_T0': (1e3, 1e5),
              'IGM_xHI': (0., 1.01)}
    for feature, (lo, hi) in ranges.items():
        for name, row in _obs(feature, h=0.7).iterrows():
            assert ((row['Data'][:, 3] > lo) & (row['Data'][:, 3] < hi)).all(), (feature, name)


def test_igm_conversions():
    # h-scaled comoving lengths, comoving lengths, h scaling, log values and limits
    mfp = _obs('IGM_MFP', h=0.678)
    gaikwad = mfp.loc['Gaikwad2023']
    assert np.isclose(gaikwad['Data'][-1, 3], 8.318 / 0.678 / 7.0)          # h^-1 cMpc -> pMpc at z = 6
    assert np.array_equal(gaikwad['LLimit'], np.arange(12) < 3)            # z < 5.2 lower limits
    davies = mfp.loc['Davies2024']
    assert np.isclose(davies['Data'][4, 3], 29.0 / 6.4)                    # cMpc -> pMpc at z = 5.4
    assert list(davies['LLimit'].nonzero()[0]) == [0, 1, 2] and list(davies['ULimit'].nonzero()[0]) == [10, 11]
    cain = mfp.loc['Cain2024_FlexRT']                                       # simulation curve, h^-1 cMpc
    assert cain['DataType'] == 'data_sim' and np.isclose(cain['Data'][0, 3], 81.012 / 0.678 / 5.7998)
    w07 = _obs('IGM_MFP', h=0.7).loc['Worseck2014']['Data']
    w10 = _obs('IGM_MFP', h=1.0).loc['Worseck2014']['Data']
    assert np.isclose(w07[-1, 3], 10.3) and np.allclose(w10[:, 3:6], w07[:, 3:6] * 0.7)   # h70^-1 pMpc
    gamma = _obs('IGM_Gamma')
    assert np.allclose(gamma.loc['Calverley2011']['Data'][:, 3], 10**np.array([-0.15, -0.84]))  # log Gamma -> 1e-12 s^-1
    assert list(gamma.loc['Bolton2007']['ULimit']) == [False, True]
    assert np.allclose(gamma.loc['Becker2013']['Data'][0, 3:], 10**np.array([0.015, 0.147, -0.131]))  # log, total errors
    assert np.allclose(gamma.loc['Qin2025']['Data'][0, 3:], [0.66, 0.70, 0.64])   # z = 4.90: MAP, 84th, 16th
    assert np.allclose(mfp.loc['Qin2025']['Data'][6, 3:], [3.32, 3.89, 3.46])     # z = 5.64: MAP below the 16th
    gao07 = _obs('IGM_MFP', h=0.7).loc['Gao2025']['Data']
    gao10 = _obs('IGM_MFP', h=1.0).loc['Gao2025']['Data']
    assert np.isclose(gao07[0, 3], 65.4) and np.allclose(gao10[:, 3:6], gao07[:, 3:6] * 0.7)  # h70^-1 pMpc
    emis = _obs('IGM_Emissivity')
    assert np.allclose(emis.loc['Becker2021']['Data'][0, 3:], 1e51 * 10**np.array([0., 0.16, -0.17]))  # log Ndot/1e51
    assert np.isclose(emis.loc['Qin2026_Fiducial']['Data'][-1, 3], 4.5926e50)  # Ndot/1e50 at z = 5.0065
    assert np.allclose(emis.loc['Cain2025']['Data'][2, 3:], [8.70e50, 9.85e50, 6.15e50])  # z = 5, asymmetric errors
    assert emis.loc['Cain2025']['DataType'] == 'data' and emis.loc['Cain2025_LateStartLateEnd']['DataType'] == 'data_sim'
    assert np.isclose(emis.loc['Asthana2024_Fiducial']['Data'][0, 3], 4.3431e50)  # z = 5, read off Fig. 3(c)
    gangolli = mfp.loc['Gangolli2025_NionLUV']                               # simulation curve, h^-1 cMpc
    assert np.isclose(gangolli['Data'][0, 3], 81.06 / 0.678 / 5.8007)
    xhi = _obs('IGM_xHI')
    q25 = xhi.loc['Qin2025']
    assert list(q25['ULimit'].nonzero()[0]) == [0, 1, 2, 3, 4] and np.allclose(q25['Data'][0, 3:], 0.005)  # rounds to 0.00
    assert np.allclose(q25['Data'][19, 3:], [0.47, 0.50, 0.48])               # z = 7.60: MAP below the 16th
    assert np.allclose(xhi.loc['Gaikwad2023']['Data'][0, 3:], [2.534e-5, 2.534e-5 + 5.85e-6, 2.534e-5 - 4.44e-6])
    assert xhi.loc['Qin2026_Fiducial']['DataType'] == 'data_sim'
    assert (xhi.loc['Cain2025_EarlyStartLateEnd']['Data'][:, 3] >= 1e-3).all()  # linear-axis curve, cut at 1e-3
    assert list(xhi.loc['Bolton2007']['LLimit']) == [False] * 4 + [True]   # z = 6: lower limit (Gamma upper limit)
    t0 = _obs('IGM_T0')
    assert np.isclose(t0.loc['Bolton2012']['Data'][0, 3], 10**4.21)
    assert np.isclose(t0.loc['Walther2019']['Data'][14, 3], 5330.)         # z = 5.0, 0.533e4 K
    assert np.allclose(t0.loc['Hiss2018']['Data'][0, 3:], [13721., 13721. + 1694., 13721. - 2152.])  # asymmetric


def test_igm_redshift_range():
    o = _obs('IGM_MFP', z_min=5.5, z_max=6.0)
    for name, row in o.iterrows():
        assert ((row['Data'][:, 0] >= 5.5) & (row['Data'][:, 0] <= 6.0)).all(), name
    assert 'Worseck2014' not in o.index        # its points end at z = 5.16


def test_available_redshifts_ignores_evolution_files():
    from astrodatapy.misc_tools import available_redshifts
    for feature in IGM_FEATURES:
        assert len(available_redshifts(feature)) == 0
