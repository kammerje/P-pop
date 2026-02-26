"""
# =============================================================================
# P-POP
# A Monte-Carlo tool to simulate exoplanet populations
# =============================================================================
"""


# =============================================================================
# STAR
# =============================================================================

class Star():
    
    def __init__(self,
                 Name,
                 Dist, # pc
                 Stype,
                 Rad, # Rsun
                 Teff, # K
                 Mass, # Msun
                 RA, # deg
                 Dec, # deg
                 Vmag=None, # mag
                 Jmag=None, # mag
                 Hmag=None, # mag
                 WDSsep=None, # au
                 WDSdmag=None, # mag
                 lGal=None, # deg
                 bGal=None): # deg
        """
        Parameters
        ----------
        Name: str
            Name of the star.
        Dist: float
            Distance (pc) of the star.
        Stype: char
            Spectral type of the star.
        Rad: float
            Radius (Rsun) of the star.
        Teff: float
            Effective temperature (K) of the star.
        Mass: float
            Mass (Msun) of the star.
        RA: float
            Right ascension (deg) of the star.
        Dec: float
            Declination (deg) of the star.
        Vmag: None, float
            V band magnitude (mag) of the star.
        Jmag: None, float
            J band magnitude (mag) of the star.
        Hmag: None, float
            H band magnitude (mag) of the star.
        WDSsep: None, float
            Separation (au) of stellar companion.
        WDSdmag: None, float
            Delta magnitude (mag) of stellar companion.
        lGal: None, float
            Galactic longitude (deg) of the star.
        bGal: None, float
            Galactic latitude (deg) of the star.
        """
        
        self.Name = Name
        self.Dist = Dist # pc
        self.Stype = Stype
        self.Rad = Rad # Rsun
        self.Teff = Teff # K
        self.Mass = Mass # Msun
        self.RA = RA # deg
        self.Dec = Dec # deg
        self.Vmag = Vmag # mag
        self.Jmag = Jmag # mag
        self.Hmag = Hmag # mag
        self.WDSsep = WDSsep # au
        self.WDSdmag = WDSdmag # mag
        self.lGal = lGal # deg
        self.bGal = bGal # deg
        
        pass
