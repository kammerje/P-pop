"""
# =============================================================================
# P-pop
# A Monte-Carlo tool to simulate exoplanet populations
# =============================================================================
"""


# =============================================================================
# IMPORTS
# =============================================================================

import astropy.table as at
import csv
import matplotlib.pyplot as plt
import numpy as np
import os

from astropy.coordinates import SkyCoord
import astropy.units as u


# =============================================================================
# HPIC
# =============================================================================

class StarCatalog():
    
    def __init__(self,
                 Stypes=['B', 'A', 'F', 'G', 'K', 'M'],
                 Dist_range=[0, 30], # pc
                 Dec_range=[-90, 90], # deg
                 Teff_range=[0., np.inf], # K
                 Path=os.path.join(os.path.split(os.path.abspath(__file__))[0], 'HPIC.csv')):
        """
        Parameters
        ----------
        Stypes: list
            Spectral types to be included.
        Dist_range: list
            Distance range (pc) to be included.
        Dec_range: list
            Declination range (deg) to be included.
        Teff_range: list
            Effective temperature range (K) to be included.
        Path: str
            Path of HPIC.
        """
        
        self.SC = self.read(Stypes,
                            Dist_range,
                            Dec_range,
                            Teff_range,
                            Path)
        
        pass
    
    def read(self,
             Stypes=['B', 'A', 'F', 'G', 'K', 'M'],
             Dist_range=[0, 30], # pc
             Dec_range=[-90, 90], # deg
             Teff_range=[0., np.inf], # K
             Path='StarCatalogs/HPIC.csv'):
        """
        Parameters
        ----------
        Stypes: list
            Spectral types to be included.
        Dist_range: list
            Distance range (pc) to be included.
        Dec_range: list
            Declination range (deg) to be included.
        Teff_range: list
            Effective temperature range (K) to be included.
        Path: str
            Path of HPIC.
        
        Returns
        -------
        SC: astropy table
            Table of stars with the following columns (if provided):
            - Name: name,
            - Dist: distance (pc),
            - Stype: spectral type,
            - Rad: radius (Rsun),
            - Teff: effective temperature (K),
            - Mass: mass (Msun),
            - RA: right ascension (deg),
            - Dec: declination (deg),
            - Vmag: V band magnitude (mag),
            - Jmag: J band magnitude (mag),
            - Hmag: H band magnitude (mag),
            - WDSsep: separation of stellar companion (au),
            - WDSdmag: delta magnitude of stellar companion (mag).
            - lGal: Galactic longitude (deg).
            - bGal: Galactic latitude (deg).
        """
        
        # Print.
        print('--> Reading star catalog HPIC.csv')
        
        # Read input star catalog.
        Nin = 0
        Name = []
        Dist = [] # pc
        Stype = []
        Rad = [] # Rsun
        Teff = [] # K
        Mass = [] # Msun
        Age = [] # Gyr
        RA = [] # deg
        Dec = [] # deg
        Vmag = [] # mag
        Jmag = [] # mag
        Hmag = [] # mag
        Kmag = [] # mag
        WDSsep = [] # au
        WDSdmag = [] # mag
        lGal = [] # deg
        bGal = [] # deg
        with open(Path) as csvfile:
            SC_in = csv.reader(csvfile, delimiter=',')
            for i, row in enumerate(SC_in):
                if (i == 0):
                    
                    #
                    row = np.array(row)
                    ColName = np.where(row == 'star_name')[0][0]
                    ColDist = np.where(row == 'sy_dist')[0][0]
                    ColStype = np.where(row == 'st_spectype')[0][0]
                    ColRad = np.where(row == 'st_rad')[0][0]
                    ColTeff = np.where(row == 'st_teff')[0][0]
                    ColMass = np.where(row == 'st_mass')[0][0]
                    ColAge = np.where(row == 'st_age')[0][0]
                    ColRA = np.where(row == 'ra')[0][0]
                    ColDec = np.where(row == 'dec')[0][0]
                    ColVmag = np.where(row == 'sy_vmag')[0][0]
                    ColJmag = np.where(row == 'sy_jmag')[0][0]
                    ColHmag = np.where(row == 'sy_hmag')[0][0]
                    ColKmag = np.where(row == 'sy_kmag')[0][0]
                    ColBinWDS = np.where(row == 'known_binary_fl')[0][0]
                    ColBinGaia = np.where(row == 'gaia_binary_fl')[0][0]
                    ColWDSsep = np.where(row == 'wds_sep')[0][0]
                    ColWDSdmag = np.where(row == 'wds_deltamag')[0][0]
                
                elif (i > 0):
                    
                    #
                    Nin += 1
                    tempName = str(row[ColName])
                    if row[ColDist] == 'null':
                        Nin -= 1
                        continue
                    else:
                        tempDist = float(row[ColDist].replace(',', '.')) # pc
                    if row[ColStype] == 'null' or ('V' in row[ColStype] and 'IV' not in row[ColStype]):
                        tempStype = str(row[ColStype])
                    else:
                        Nin -= 1
                        continue
                    if row[ColRad] == 'null':
                        Nin -= 1
                        continue
                    else:
                        tempRad = float(row[ColRad].replace(',', '.')) # Rsun
                    tempTeff = float(row[ColTeff].replace(',', '.')) # K
                    tempMass = float(row[ColMass].replace(',', '.')) # Msun
                    if row[ColAge] == 'null':
                        tempAge = np.nan
                    else:
                        tempAge = float(row[ColAge].replace(',', '.')) # Gyr
                    tempRA = float(row[ColRA].replace(',', '.')) # deg
                    tempDec = float(row[ColDec].replace(',', '.')) # deg
                    if row[ColVmag] == 'null':
                        tempVmag = np.nan
                    else:
                        tempVmag = float(row[ColVmag].replace(',', '.')) # mag
                    if row[ColJmag] == 'null':
                        tempJmag = np.nan
                    else:
                        tempJmag = float(row[ColJmag].replace(',', '.')) # mag
                    if row[ColHmag] == 'null':
                        tempHmag = np.nan
                    else:
                        tempHmag = float(row[ColHmag].replace(',', '.')) # mag
                    if row[ColKmag] == 'null':
                        tempKmag = np.nan
                    else:
                        tempKmag = float(row[ColKmag].replace(',', '.')) # mag
                    if row[ColWDSsep] == 'null':
                        tempWDSsep = np.inf
                    else:
                        # tempWDSsep = float(row[ColWDSsep].replace(',', '.')) # arcsec
                        tempWDSsep = float(row[ColWDSsep].replace(',', '.')) * tempDist # au
                    if row[ColWDSdmag] == 'null':
                        tempWDSdmag = np.nan
                    else:
                        tempWDSdmag = float(row[ColWDSdmag].replace(',', '.')) # mag
                    c_icrs = SkyCoord(ra=tempRA * u.deg, dec=tempDec * u.deg, frame='icrs')
                    c_galactic = c_icrs.transform_to('galactic')
                    templGal = float(c_galactic.l.value) # deg
                    tempbGal = float(c_galactic.b.value) # deg
                    
                    #
                    Name += [tempName]
                    Dist += [tempDist] # pc
                    Stype += [tempStype[0]]
                    Rad += [tempRad] # Rsun
                    Teff += [tempTeff] # K
                    Mass += [tempMass] # Msun
                    Age += [tempAge] # Gyr
                    RA += [tempRA] # deg
                    Dec += [tempDec] # deg
                    Vmag += [tempVmag] # mag
                    Jmag += [tempJmag] # mag
                    Hmag += [tempHmag] # mag
                    Kmag += [tempKmag] # mag
                    WDSsep += [tempWDSsep] # au
                    WDSdmag += [tempWDSdmag] # mag
                    lGal += [templGal] # deg
                    bGal += [tempbGal] # deg
        
        # Convert lists to arrays.
        Name = np.array(Name)
        Dist = np.array(Dist) # pc
        Stype = np.array(Stype)
        Rad = np.array(Rad) # Rsun
        Teff = np.array(Teff) # K
        Mass = np.array(Mass) # Msun
        Age = np.array(Age) # Gyr
        RA = np.array(RA) # deg
        Dec = np.array(Dec) # deg
        Vmag = np.array(Vmag) # mag
        Jmag = np.array(Jmag) # mag
        Hmag = np.array(Hmag) # mag
        Kmag = np.array(Kmag) # mag
        WDSsep = np.array(WDSsep) # au
        WDSdmag = np.array(WDSdmag) # mag
        lGal = np.array(lGal) # deg
        bGal = np.array(bGal) # deg
        
        # Create output star catalog.
        SC_out = at.Table(names=('Name', 'Dist', 'Stype', 'Rad', 'Teff', 'Mass', 'Age', 'RA', 'Dec', 'Vmag', 'Jmag', 'Hmag', 'Kmag', 'WDSsep', 'WDSdmag', 'lGal', 'bGal'),
                          dtype = ('S32', 'd', 'c', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd'))
        
        for i in range(len(Name)):
            
            # Check whether stellar type fits.
            temp1 = (Stypes is None) or (Stype[i] in Stypes)
            
            # Check whether distance fits.
            temp2 = (Dist_range is None) or (Dist_range[0] <= Dist[i] <= Dist_range[1])
            
            # Check whether declination fits.
            temp3 = (Dec_range is None) or (Dec_range[0] <= Dec[i] <= Dec_range[1])
            
            # Check whether effective temperature fits.
            temp4 = (Teff_range is None) or (Teff_range[0] <= Teff[i] <= Teff_range[1])
            
            # Fill output star catalog.
            if (temp1 and temp2 and temp3 and temp4):
                SC_out.add_row([Name[i], Dist[i], Stype[i], Rad[i], Teff[i], Mass[i], Age[i], RA[i], Dec[i], Vmag[i], Jmag[i], Hmag[i], Kmag[i], WDSsep[i], WDSdmag[i], lGal[i], bGal[i]])
        
        # Print.
        text1 = len(SC_out)/float(Nin)*100.
        print('--> Including %.0f = %.2f%% stars' % (len(SC_out), text1))
        text2 = np.unique(np.array(Stype))
        if (Stypes is None):
            print('--> Including spectral types '+str(text2))
        else:
            text3 = [f for f in text2 if f not in Stypes]
            print('--> Including spectral types '+str(Stypes)+', excluding spectral types '+str(text3))
        text4 = np.min(SC_out['Dist'])
        text5 = np.max(SC_out['Dist'])
        if (Dist_range is None):
            print('--> Distance in [%.2f, %.2f] pc' % (text4, text5))
        else:
            print('--> Distance in [%.2f, %.2f] pc' % (text4, text5)+', distance limits [%.2f, %.2f] pc' % (Dist_range[0], Dist_range[1]))
        text6 = np.min(SC_out['Dec'])
        text7 = np.max(SC_out['Dec'])
        if (Dec_range is None):
            print('--> Declination in [%.2f, %.2f] deg' % (text6, text7))
        else:
            print('--> Declination in [%.2f, %.2f] deg' % (text6, text7)+', declination limits [%.2f, %.2f] deg' % (Dec_range[0], Dec_range[1]))
        text8 = np.min(SC_out['Teff'])
        text9 = np.max(SC_out['Teff'])
        if (Teff_range is None):
            print('--> Effective temperature in [%.1f, %.1f] K' % (text8, text9))
        else:
            print('--> Effective temperature in [%.1f, %.1f] K' % (text8, text9)+', effective temperature limits [%.1f, %.1f] K' % (Teff_range[0], Teff_range[1]))
        
        # # Constants.
        # G = 6.674e-11 # m^3/kg/s^2
        # Rsun = 695700000. # m
        # Msun = 1.989e30 # kg
        # au = 149597870700. # m
        # sigma = 5.670e-8 # W/m^2/K^4
        
        # # Variables.
        # Porb = 20000.
        # ep = 0.
        # thetap = 0.
        
        # # Computations.
        # ap = ((G*SC_out['Mass']*Msun*(Porb*86400.)**2)/(4.*np.pi**2))**(1./3.)/au
        # rp = ap*(1.-ep**2)/(1.+ep*np.cos(thetap))
        # Fp = sigma*SC_out['Teff']**4*(SC_out['Rad']*Rsun)**2/(rp*au)**2/1361.

        # import pdb; pdb.set_trace()
        
        return SC_out
    
    def SummaryPlot(self,
                    FigDir=None,
                    block=True):
        """
        Parameters
        ----------
        FigDir: str
            Directory to which summary plots are saved.
        block: bool
            If True, blocks plots when showing.
        """
        
        f, ax = plt.subplots(2, 3)
        ax[0, 0].hist(self.SC['Dist'], bins=25)
        ax[0, 0].grid(axis='y')
        ax[0, 0].set_xlabel('Distance [pc]')
        ax[0, 0].set_ylabel('Number')
        ax[0, 1].hist(self.SC['Mass'], bins=25)
        ax[0, 1].grid(axis='y')
        ax[0, 1].set_xlabel('Mass [$M_\odot$]')
        ax[0, 1].set_ylabel('Number')
        ax[0, 2].hist(self.SC['Rad'], bins=25)
        ax[0, 2].grid(axis='y')
        ax[0, 2].set_xlabel('Radius [$R_\odot$]')
        ax[0, 2].set_ylabel('Number')
        ax[1, 0].hist(self.SC['Teff'], bins=25)
        ax[1, 0].grid(axis='y')
        ax[1, 0].set_xlabel('Effective temperature [K]')
        ax[1, 0].set_ylabel('Number')
        ax[1, 1].scatter(self.SC['Teff'], self.SC['Mass'], c=self.SC['Teff'], cmap='jet_r', s=2)
        ax[1, 1].invert_xaxis()
        ax[1, 1].grid()
        ax[1, 1].set_xlabel('Effective temperature [K]')
        ax[1, 1].set_ylabel('Mass [$M_\odot$]')
        plt.subplot(236, projection='aitoff')
        plt.grid()
        plt.scatter(((self.SC['RA']*np.pi/180.+np.pi) % 2.*np.pi)-np.pi, self.SC['Dec']*np.pi/180., s=2)
        plt.xlabel('Right ascension [deg]')
        plt.ylabel('Declination [deg]')
        plt.tight_layout()
        if (FigDir is not None):
            plt.savefig(FigDir+'StarCatalog.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
