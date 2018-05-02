def _Mbh2Ledd(Mass, EddingtonRatio=1.0):
    '''Black hole mass in Msol to Eddington luminosity in Lsol'''
    return EddingtonRatio*Mass*32731.9788953923

def _Lbol2Mbol(Lbol):
    '''bolometric luminosity to bolometric magnitude'''
    return 4.74 - 2.5*np.log10(Lbol)

def _kBvega(Lbol):
    '''bolometric correction for B band in Vega'''
    return 6.25 * (Lbol/1e10)**-0.37 + 9.00 * (Lbol/1e10)**-0.012

def _dkBdLbol(Lbol):
    return -2.3125e-10 * (Lbol/1e10)**-1.37 - 0.108e-10 * (Lbol/1e10)**-1.012

def _dlogkBdlogLbol(Lbol):
    return Lbol/_kBvega(Lbol)*_dkBdLbol(Lbol)

def _Lbol2MBvega(Lbol):
    '''bolometric luminosity to B magnitude in Vega'''
    return _Lbol2Mbol(Lbol) + 2.5*np.log10(_kBvega(Lbol))

def _k15um(Lbol):
    '''bolometric correction for 15 um'''
    return 7.40 * (Lbol/1e10)**-0.37 + 10.66 * (Lbol/1e10)**-0.014

def _Lbol2M15um(Lbol):
    '''bolometric luminosity to 15um in Vega'''
    return _Lbol2Mbol(Lbol) + 2.5*np.log10(_k15um(Lbol))

def _ksoftX(Lbol):
    '''bolometric correction for 0.5-2keV'''
    return 17.87 * (Lbol/1e10)**0.28 + 10.03 * (Lbol/1e10)**-0.020

def _Lbol2MsoftX(Lbol):
    '''bolometric luminosity to 0.5-2keV'''
    return _Lbol2Mbol(Lbol) + 2.5*np.log10(_ksoftX(Lbol))

def _khardX(Lbol):
    '''bolometric correction for 2-10keV'''
    return 10.83 * (Lbol/1e10)**0.28 + 6.08 * (Lbol/1e10)**-0.020

def _Lbol2MhardX(Lbol):
    '''bolometric luminosity to 2-10keV'''
    return _Lbol2Mbol(Lbol) + 2.5*np.log10(_khardX(Lbol))

def _MBvega2MB(MBvega):
    '''convert B magnitude from Vega to AB system'''
    return MBvega - 0.09

def _MB2M1450(MB):
    '''convert B magnitede to UV 1450 magnitude (all AB system) assuming spetral index is 0.44'''
    return MB + 0.524

def _M14502MB(M1450):
    '''convert UV 1450 magnitede to B magnitude (all AB system) assuming spetral index is 0.44'''
    return M1450 - 0.524

def _MBvega2M1450(MBvega):
    '''convert B magnitede in Vega to UV 1450 magnitude assuming spetral index is 0.44'''
    return _MB2M1450(_MBvega2MB(MBvega))

def _Lbol2MB(Lbol):
    '''bolometric luminosity to B magnitude in AB system'''
    return _MBvega2MB(_Lbol2MBvega(Lbol))

def _Lbol2M1450(Lbol):
    '''bolometric luminosity to UV 1450 magnitude in AB system assuming spetral index is 0.44'''
    return _MB2M1450(_Lbol2MB(Lbol))

def _sumMag(M1,M2):
    '''add two absolute AB magnitude together'''
    return -2.5 * np.log10(10.**(-M1 / 2.5) + 10.**(-M2 / 2.5))
