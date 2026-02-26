"""
# =============================================================================
# P-POP
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


# =============================================================================
# LTC (VERSION 3)
# =============================================================================

class StarCatalog():
    
    def __init__(self,
                 Stypes=['B', 'A', 'F', 'G', 'K', 'M', 'D'],
                 Dist_range=[0, 30], # pc
                 Dec_range=[-90, 90], # deg
                 Teff_range=[0., np.inf], # K
                 Path=os.path.join(os.path.split(os.path.abspath(__file__))[0], 'LTC3.csv')):
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
            Path of LTC (version 3).
        """
        
        self.SC = self.read(Stypes,
                            Dist_range,
                            Dec_range,
                            Teff_range,
                            Path)
        
        pass
    
    def read(self,
             Stypes=['B', 'A', 'F', 'G', 'K', 'M', 'D'],
             Dist_range=[0, 30], # pc
             Dec_range=[-90, 90], # deg
             Teff_range=[0., np.inf], # K
             Path='StarCatalogs/LTC3.csv'):
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
            Path of LTC (version 3).
        
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
        print('--> Reading star catalog LTC3.csv')
        
        # Read input star catalog.
        Nin = 0
        Name = []
        Dist = [] # pc
        Stype = []
        Rad = [] # Rsun
        Teff = [] # K
        Mass = [] # Msun
        RA = [] # deg
        Dec = [] # deg
        Vmag = [] # mag
        Jmag = [] # mag
        Hmag = [] # mag
        WDSsep = [] # au
        WDSdmag = [] # mag
        lGal = [] # deg
        bGal = [] # deg
        with open(Path) as csvfile:
            SC_in = csv.reader(csvfile, delimiter=',')
            for i, row in enumerate(SC_in):
                if (i == 39):
                    
                    #
                    row = np.array(row)
                    ColName = np.where(row == 'sim_name')[0][0]
                    ColDist = np.where(row == 'distance')[0][0]
                    ColStype = np.where(row == 'sim_sptype')[0][0]
                    ColRad = np.where(row == 'mod_R')[0][0]
                    ColTeff = np.where(row == 'mod_Teff')[0][0]
                    ColMass = np.where(row == 'mod_M')[0][0]
                    ColRA = np.where(row == 'sim_ra')[0][0]
                    ColDec = np.where(row == 'sim_dec')[0][0]
#                    ColVmag = np.where(row == 'st_vmag')[0][0]
#                    ColJmag = np.where(row == 'st_j2m')[0][0]
#                    ColHmag = np.where(row == 'st_h2m')[0][0]
                    ColWDSsep = np.where(row == 'sep_phys')[0][0]
#                    ColWDSdmag = np.where(row == 'wds_deltamag')[0][0]
                    CollGal = np.where(row == 'gal_coord_l')[0][0]
                    ColbGal = np.where(row == 'gal_coord_b')[0][0]
                
                elif (i > 39):
                    
                    #
                    Nin += 1
                    tempName = str(row[ColName])
                    tempDist = float(row[ColDist]) # pc
                    tempStype = str(row[ColStype])
                    if (len(tempStype) == 0):
                        continue
                    tempRad = float(row[ColRad]) # Rsun
                    tempTeff = float(row[ColTeff]) # K
                    tempMass = float(row[ColMass]) # Msun
#                    tempRA = str(row[ColRA]) # deg
#                    tempRA = (int(tempRA[0:2])+int(tempRA[3:5])/60.+float(tempRA[6:11])/3600.)*360./24. # deg
                    tempRA = float(row[ColRA]) # deg
#                    tempDec = str(row[ColDec]) # deg
#                    temp = int(tempDec[0:3])
#                    if (temp >= 0.):
#                        tempDec = float(temp)+int(tempDec[4:6])/60.+float(tempDec[7:11])/3600. # deg
#                    else:
#                        tempDec = float(temp)-int(tempDec[4:6])/60.-float(tempDec[7:11])/3600. # deg
                    tempDec = float(row[ColDec]) # deg
#                    tempVmag = float(row[ColVmag]) # mag
                    tempVmag = np.nan
#                    tempJmag = float(row[ColJmag]) # mag
                    tempJmag = np.nan
#                    tempHmag = float(row[ColHmag]) # mag
                    tempHmag = np.nan
#                    try:
#                        tempWDSsep = float(row[ColWDSsep]) # arcsec
#                    except:
#                        tempWDSsep = np.inf
                    if (row[ColWDSsep] != 'nan'):
                        tempWDSsep = float(row[ColWDSsep]) # au
                    else:
                        tempWDSsep = np.inf
#                    try:
#                        tempWDSdmag = float(row[ColWDSdmag]) # mag
#                    except:
#                        tempWDSdmag = np.inf
                    tempWDSdmag = np.nan
                    templGal = float(row[CollGal]) # deg
                    tempbGal = float(row[ColbGal]) # deg
                    
                    #
                    Name += [tempName]
                    Dist += [tempDist] # pc
                    Stype += [tempStype[0]]
                    Rad += [tempRad] # Rsun
                    Teff += [tempTeff] # K
                    Mass += [tempMass] # Msun
                    RA += [tempRA] # deg
                    Dec += [tempDec] # deg
                    Vmag += [tempVmag] # mag
                    Jmag += [tempJmag] # mag
                    Hmag += [tempHmag] # mag
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
        RA = np.array(RA) # deg
        Dec = np.array(Dec) # deg
        Vmag = np.array(Vmag) # mag
        Jmag = np.array(Jmag) # mag
        Hmag = np.array(Hmag) # mag
        WDSsep = np.array(WDSsep) # au
        WDSdmag = np.array(WDSdmag) # mag
        lGal = np.array(lGal) # deg
        bGal = np.array(bGal) # deg
        
        # Create output star catalog.
        SC_out = at.Table(names=('Name', 'Dist', 'Stype', 'Rad', 'Teff', 'Mass', 'RA', 'Dec', 'Vmag', 'Jmag', 'Hmag', 'WDSsep', 'WDSdmag', 'lGal', 'bGal'),
                          dtype = ('S32', 'd', 'c', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd'))
        
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
                SC_out.add_row([Name[i], Dist[i], Stype[i], Rad[i], Teff[i], Mass[i], RA[i], Dec[i], Vmag[i], Jmag[i], Hmag[i], WDSsep[i], WDSdmag[i], lGal[i], bGal[i]])
        
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
