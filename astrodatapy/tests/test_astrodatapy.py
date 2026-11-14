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

def test_correlation():
    from astrodatapy.correlation import correlation as function
    features = ['Magorrian', 'Tully_Fisher', 'DiskSize_StellarMass',
                'GasFraction_StellarMass', 'sSFR_StellarMass_Blue',
                'HaloMass_StellarMass', 'HaloMass_StellarMass_Blue',
                'HaloMass_StellarMass_Red']
    for feature in features:
        for z_target in available_redshifts(feature):
            assert function(feature, z_target)

def test_clustering():
    from astrodatapy.clustering import clustering as function
    features = ['QC_2PTCF',]
    for feature in features:
        for z_target in available_redshifts(feature):
            assert function(feature, z_target)
