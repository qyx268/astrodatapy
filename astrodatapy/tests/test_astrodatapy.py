from astrodatapy.misc_tools import available_redshifts

def test_number_density():
    from astrodatapy.number_density import number_density as function
    features = ['GSMF', 'GSMF_Blue', 'GSMF_Red', 'GSMF_Bulge', 'GSMF_Disk',\
                'GSMF_Quiescent', 'BHM', 'BHMF', 'QLF_bolometric',\
                'QLF_optical', 'QLF_UV', 'SFRF', 'GLF_UV']

    for feature in features:
        for z_target in available_redshifts(feature):
            assert function(feature, z_target)

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
