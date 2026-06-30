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
import matplotlib.pyplot as plt
import numpy as np
import os


# =============================================================================
# CROSSFIELDBRIGHTSAMPLE
# =============================================================================

class StarCatalog():
    """
    https://ui.adsabs.harvard.edu/abs/2013A%26A...551A..99C/abstract
    """
    
    def __init__(self,
                 Stypes=['A', 'F', 'G', 'K', 'M'],
                 Dist_range=[0, 20], # pc
                 Dec_range=[-90, 90], # deg
                 Teff_range=[0., np.inf], # K
                 Path=os.path.join(os.path.split(os.path.abspath(__file__))[0], 'CrossfieldBrightSample.tbl')):
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
            Path of CrossfieldBrightSample.
        """
        
        self.SC = self.read(Stypes,
                            Dist_range,
                            Dec_range,
                            Teff_range,
                            Path)
        
        pass
    
    def read(self,
             Stypes=['A', 'F', 'G', 'K', 'M'],
             Dist_range=[0, 20], # pc
             Dec_range=[-90, 90], # deg
             Teff_range=[0., np.inf], # K
             Path='StarCatalogs/CrossfieldBrightSample.tbl'):
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
            Path of CrossfieldBrightSample.
        
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
            - WDSsep: separation of stellar companion (arcsec),
            - WDSdmag: delta magnitude of stellar companion (mag).
        """
        
        # Print.
        print('--> Reading star catalog CrossfieldBrightSample.tbl')
        
        # Read input star catalog.
        SC_in = at.Table.read(Path, format='ipac')
        
        # Create output star catalog.
        SC_out = at.Table(names=('Name', 'Dist', 'Stype', 'Rad', 'Teff', 'Mass', 'RA', 'Dec'),\
                          dtype=('S64', 'd', 'c', 'd', 'd', 'd', 'd', 'd'))
        
        Stype = []
        Dist = [] # pc
        Dec = [] # deg
        for i in range(len(SC_in)):
            
            # Check whether stellar type fits.
            Stype += [str(SC_in[i]['Sptype'])[0]]
            temp1 = (Stypes is None) or (Stype[-1] in Stypes)
            
            # Check whether distance fits.
            Dist += [1./(1e-3*float(SC_in[i]['pi_trig']))] # pc
            temp2 = (Dist_range is None) or (Dist_range[0] <= Dist[-1] <= Dist_range[1])
            
            # Check whether declination fits.
            temp = int(SC_in[i]['DEC_2000'][0:3]) # deg
            if (temp >= 0):
                Dec += [float(temp)+int(SC_in[i]['DEC_2000'][4:6])/60.+float(SC_in[i]['DEC_2000'][7:14])/3600.] # deg
            else:
                Dec += [float(temp)-int(SC_in[i]['DEC_2000'][4:6])/60.-float(SC_in[i]['DEC_2000'][7:14])/3600.] # deg
            temp3 = (Dec_range is None) or (Dec_range[0] <= Dec[-1] <= Dec_range[1])
            
            # Check whether effective temperature fits.
            temp4 = (Teff_range is None) or (Teff_range[0] <= SC_in[i]['teff'] <= Teff_range[1])
            
            # Fill output star catalog.
            if (temp1 and temp2 and temp3 and temp4):
                RA = (int(SC_in[i]['RA_2000'][0:2])+int(SC_in[i]['RA_2000'][3:5])/60.+float(SC_in[i]['RA_2000'][6:13])/3600.)*360./24. # deg
                SC_out.add_row([SC_in[i]['DiscoveryName'], Dist[-1], Stype[-1], SC_in[i]['r'], SC_in[i]['teff'], SC_in[i]['m'], RA, Dec[-1]])
        
        # Print.
        text1 = len(SC_out)/float(i+1)*100.
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
