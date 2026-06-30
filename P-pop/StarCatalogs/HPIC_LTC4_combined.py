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


# =============================================================================
# LTC (VERSION 4)
# =============================================================================

class StarCatalog():
    
    def __init__(self,
                 Stypes=['B', 'A', 'F', 'G', 'K', 'M', 'D'],
                 Dist_range=[0, 30], # pc
                 Dec_range=[-90, 90], # deg
                 Teff_range=[0., np.inf], # K
                 Path=os.path.join(os.path.split(os.path.abspath(__file__))[0], 'HPIC_LTC4_combined.txt')):
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
            Path of LTC (version 4).
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
             Path='StarCatalogs/HPIC_LTC4_combined.txt'):
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
            Path of LTC (version 4).
        
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
        print('--> Reading star catalog HPIC_LTC4_combined.txt')
        
        # Read input star catalog.
        Nin = 0
        Name_all = []
        NiceName_all = []
        RejectReason = []
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
        Fp_optHZin = [] # d
        Fp_optHZout = [] # d
        Fp_conHZin = [] # d
        Fp_conHZout = [] # d
        Porb_optHZin = [] # d
        Porb_optHZout = [] # d
        Porb_conHZin = [] # d
        Porb_conHZout = [] # d
        with open(Path) as txtfile:
            lines = txtfile.readlines()
            for i, row in enumerate(lines):
                row = row.split('|')
                if (i == 0):
                    
                    #
                    row = np.array(row)
                    ColName = np.where(row == 'star_name')[0][0]
                    ColNiceName = np.where(row == 'simbad_name')[0][0]
                    ColDist = np.where(row == 'sy_dist')[0][0]
                    ColStype = np.where(row == 'st_spectype')[0][0]
                    ColRad = np.where(row == 'st_rad')[0][0]
                    ColTeff = np.where(row == 'st_teff')[0][0]
                    ColMass = np.where(row == 'st_mass')[0][0]
                    ColRA = np.where(row == 'ra')[0][0]
                    ColDec = np.where(row == 'dec')[0][0]
                    ColWDSsep = np.where(row == 'fl_known')[0][0]
                
                elif (i > 0):
                    
                    #
                    Nin += 1
                    Name_all += [str(row[ColName])]
                    NiceName_all += [str(row[ColNiceName])]
                    RejectReason += [None]
                    tempName = str(row[ColName])
                    tempDist = float(row[ColDist]) # pc
                    tempStype = str(row[ColStype])
                    if (len(tempStype) == 0) or (tempStype == 'nan'):
                        RejectReason[-1] = 'No SpT'
                        continue
                    if (tempStype[0] == 'D'):
                        RejectReason[-1] = 'White dwarf'
                        continue
                    if (tempStype[0] == 'S'):
                        RejectReason[-1] = 'Giant star'
                        continue
                    if row[ColRad] == '':
                        RejectReason[-1] = 'No radius'
                        continue
                    tempRad = float(row[ColRad]) # Rsun
                    tempTeff = float(row[ColTeff]) # K
                    tempMass = float(row[ColMass]) # Msun
                    tempRA = float(row[ColRA]) # deg
                    tempDec = float(row[ColDec]) # deg
                    tempVmag = np.nan
                    tempJmag = np.nan
                    tempHmag = np.nan
                    tempWDSsep = np.inf # arcsec
                    tempWDSdmag = np.nan
                    
                    #
                    Name += [tempName]
                    Dist += [tempDist] # pc
                    if (tempStype[0] == 'd') or (tempStype[0] == 'k'):
                        Stype += [tempStype[1]]
                    elif (tempStype[0:2] == 'sd'):
                        Stype += [tempStype[2]]
                    elif (tempStype[0:3] == 'sd:'):
                        Stype += [tempStype[3]]
                    elif (tempStype[0:3] == 'esd'):
                        Stype += [tempStype[3]]
                    elif (tempStype[0:4] == 'd/sd'):
                        Stype += [tempStype[4]]
                    else:
                        Stype += [tempStype[0]]
                    Rad += [tempRad] # Rsun
                    Teff += [tempTeff] # K
                    Mass += [tempMass] # Msun
                    RA += [tempRA] # deg
                    Dec += [tempDec] # deg
                    WDSsep += [tempWDSsep] # au
                    WDSdmag += [tempWDSdmag] # mag
                    
                    G = 6.674e-11
                    sigma = 5.67e-8
                    Msun = 1.989e30
                    Rsun = 696340000.
                    
                    Lum = 4.*np.pi*(Rad[-1]*Rsun)**2*sigma*Teff[-1]**4
                    Tmod = Teff[-1]-5780.
                    Fmin = 0.320 + 5.547e-5 * Tmod + 1.526e-9 * Tmod**2 - 2.874e-12 * Tmod**3 - 5.011e-16 * Tmod**4
                    Fmax = 1.776 + 2.136e-4 * Tmod + 2.533e-8 * Tmod**2 - 1.332e-11 * Tmod**3 - 3.097e-15 * Tmod**4
                    Fp_optHZout += [Fmin]
                    Fp_optHZin += [Fmax]
                    amin = np.sqrt(Lum/(4.*np.pi)/Fmax/1361.)
                    amax = np.sqrt(Lum/(4.*np.pi)/Fmin/1361.)
                    Pmin = np.sqrt(4.*np.pi**2/(G*Mass[-1]*Msun)*amin**3) / 3600. / 24.
                    Pmax = np.sqrt(4.*np.pi**2/(G*Mass[-1]*Msun)*amax**3) / 3600. / 24.
                    Porb_optHZin += [Pmin]
                    Porb_optHZout += [Pmax]
                    Fmin = 0.356 + 6.171e-5 * Tmod + 1.698e-9 * Tmod**2 - 3.198e-12 * Tmod**3 - 5.575e-16 * Tmod**4
                    Fmax = 1.107 + 1.332e-4 * Tmod + 1.580e-8 * Tmod**2 - 8.308e-12 * Tmod**3 - 1.931e-15 * Tmod**4
                    Fp_conHZout += [Fmin]
                    Fp_conHZin += [Fmax]
                    amin = np.sqrt(Lum/(4.*np.pi)/Fmax/1361.)
                    amax = np.sqrt(Lum/(4.*np.pi)/Fmin/1361.)
                    Pmin = np.sqrt(4.*np.pi**2/(G*Mass[-1]*Msun)*amin**3) / 3600. / 24.
                    Pmax = np.sqrt(4.*np.pi**2/(G*Mass[-1]*Msun)*amax**3) / 3600. / 24.
                    Porb_conHZin += [Pmin]
                    Porb_conHZout += [Pmax]
        
        # Convert lists to arrays.
        Name_all = np.array(Name_all)
        NiceName_all = np.array(NiceName_all)
        RejectReason = np.array(RejectReason)
        Name = np.array(Name)
        Dist = np.array(Dist) # pc
        Stype = np.array(Stype)
        Rad = np.array(Rad) # Rsun
        Teff = np.array(Teff) # K
        Mass = np.array(Mass) # Msun
        RA = np.array(RA) # deg
        Dec = np.array(Dec) # deg
        WDSsep = np.array(WDSsep) # au
        WDSdmag = np.array(WDSdmag) # mag
        Fp_optHZin = np.array(Fp_optHZin)
        Fp_optHZout = np.array(Fp_optHZout)
        Fp_conHZin = np.array(Fp_conHZin)
        Fp_conHZout = np.array(Fp_conHZout)
        Porb_optHZin = np.array(Porb_optHZin)
        Porb_optHZout = np.array(Porb_optHZout)
        Porb_conHZin = np.array(Porb_conHZin)
        Porb_conHZout = np.array(Porb_conHZout)
        
        # with open(os.path.join(os.path.split(os.path.abspath(__file__))[0], 'HPIC_LTC4_combined_proc.csv'), 'w', newline='') as csvfile:
        #     csvwriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        #     csvwriter.writerow(['Name', 'NiceName', 'Rejected?', 'Binary?', 'SpT', 'Dist', 'Teff', 'Mass', 'Rad', 'Pmin optHZ', 'Pmax optHZ', 'Pmin conHZ', 'Pmax conHZ', 'Fmin optHZ', 'Fmax optHZ', 'Fmin conHZ', 'Fmax conHZ', 'SAG13 optHZ?', 'SAG13 conHZ?', 'SAG13Extrap optHZ?', 'SAG13Extrap conHZ?', 'Bryson2021 optHZ?', 'Bryson2021 conHZ?', 'Dressing2015 optHZ?', 'Dressing2015 conHZ?', 'Dressing2015Extrap optHZ?', 'Dressing2015Extrap conHZ?', 'Kaminski2025 optHZ?', 'Kaminski2025 conHZ?', 'Kaminski2025Extrap optHZ?', 'Kaminski2025Extrap conHZ?'])
        #     j = 0
        #     for i in range(len(Name_all)):
        #         if RejectReason[i] is None:
        #             if (Porb_optHZin[j] > 0.5) and (Porb_optHZout[j] < 500.):
        #                 SAG13_optHZ = 'TRUE'
        #             else:
        #                 SAG13_optHZ = 'FALSE'
        #             if (Porb_conHZin[j] > 0.5) and (Porb_conHZout[j] < 500.):
        #                 SAG13_conHZ = 'TRUE'
        #             else:
        #                 SAG13_conHZ = 'FALSE'
        #             if (Porb_optHZin[j] > 0.5) and (Porb_optHZout[j] < 20000.):
        #                 SAG13Extrap_optHZ = 'TRUE'
        #             else:
        #                 SAG13Extrap_optHZ = 'FALSE'
        #             if (Porb_conHZin[j] > 0.5) and (Porb_conHZout[j] < 20000.):
        #                 SAG13Extrap_conHZ = 'TRUE'
        #             else:
        #                 SAG13Extrap_conHZ = 'FALSE'
        #             if (Fp_optHZout[j] > 0.2) and (Fp_optHZin[j] < 2.2):
        #                 Bryson2021_optHZ = 'TRUE'
        #             else:
        #                 Bryson2021_optHZ = 'FALSE'
        #             if (Fp_conHZout[j] > 0.2) and (Fp_conHZin[j] < 2.2):
        #                 Bryson2021_conHZ = 'TRUE'
        #             else:
        #                 Bryson2021_conHZ = 'FALSE'
        #             if Stype[j] == 'M':
        #                 if (Porb_optHZin[j] > 0.5) and (Porb_optHZout[j] < 200.):
        #                     Dressing2015_optHZ = 'TRUE'
        #                 else:
        #                     Dressing2015_optHZ = 'FALSE'
        #                 if (Porb_conHZin[j] > 0.5) and (Porb_conHZout[j] < 200.):
        #                     Dressing2015_conHZ = 'TRUE'
        #                 else:
        #                     Dressing2015_conHZ = 'FALSE'
        #                 if (Porb_optHZin[j] > 0.5) and (Porb_optHZout[j] < 400.):
        #                     Dressing2015Extrap_optHZ = 'TRUE'
        #                 else:
        #                     Dressing2015Extrap_optHZ = 'FALSE'
        #                 if (Porb_conHZin[j] > 0.5) and (Porb_conHZout[j] < 400.):
        #                     Dressing2015Extrap_conHZ = 'TRUE'
        #                 else:
        #                     Dressing2015Extrap_conHZ = 'FALSE'
        #             else:
        #                 Dressing2015_optHZ = 'n/a'
        #                 Dressing2015_conHZ = 'n/a'
        #                 Dressing2015Extrap_optHZ = 'n/a'
        #                 Dressing2015Extrap_conHZ = 'n/a'
        #             if Stype[j] == 'M':
        #                 if (Porb_optHZin[j] > 1.) and (Porb_optHZout[j] < 100.):
        #                     Kaminski2025_optHZ = 'TRUE'
        #                 else:
        #                     Kaminski2025_optHZ = 'FALSE'
        #                 if (Porb_conHZin[j] > 1.) and (Porb_conHZout[j] < 100.):
        #                     Kaminski2025_conHZ = 'TRUE'
        #                 else:
        #                     Kaminski2025_conHZ = 'FALSE'
        #                 if (Porb_optHZin[j] > 1.) and (Porb_optHZout[j] < 1000.):
        #                     Kaminski2025Extrap_optHZ = 'TRUE'
        #                 else:
        #                     Kaminski2025Extrap_optHZ = 'FALSE'
        #                 if (Porb_conHZin[j] > 1.) and (Porb_conHZout[j] < 1000.):
        #                     Kaminski2025Extrap_conHZ = 'TRUE'
        #                 else:
        #                     Kaminski2025Extrap_conHZ = 'FALSE'
        #             else:
        #                 Kaminski2025_optHZ = 'n/a'
        #                 Kaminski2025_conHZ = 'n/a'
        #                 Kaminski2025Extrap_optHZ = 'n/a'
        #                 Kaminski2025Extrap_conHZ = 'n/a'
        #             csvwriter.writerow([Name_all[i], NiceName_all[i], 'n/a', str(np.isfinite(WDSsep[j])), Stype[j], ('%.1f' % Dist[j]).replace('.', ','), '%.0f' % Teff[j], ('%.2f' % Mass[j]).replace('.', ','), ('%.2f' % Rad[j]).replace('.', ','), ('%.2f' % Porb_optHZin[j]).replace('.', ','), ('%.2f' % Porb_optHZout[j]).replace('.', ','), ('%.2f' % Porb_conHZin[j]).replace('.', ','), ('%.2f' % Porb_conHZout[j]).replace('.', ','), ('%.2f' % Fp_optHZout[j]).replace('.', ','), ('%.2f' % Fp_optHZin[j]).replace('.', ','), ('%.2f' % Fp_conHZout[j]).replace('.', ','), ('%.2f' % Fp_conHZin[j]).replace('.', ','), SAG13_optHZ, SAG13_conHZ, SAG13Extrap_optHZ, SAG13Extrap_conHZ, Bryson2021_optHZ, Bryson2021_conHZ, Dressing2015_optHZ, Dressing2015_conHZ, Dressing2015Extrap_optHZ, Dressing2015Extrap_conHZ, Kaminski2025_optHZ, Kaminski2025_conHZ, Kaminski2025Extrap_optHZ, Kaminski2025Extrap_conHZ])
        #             j += 1
        #         else:
        #             csvwriter.writerow([Name_all[i], NiceName_all[i], RejectReason[i]])
        
        # import pdb; pdb.set_trace()
        
        # Create output star catalog.
        SC_out = at.Table(names=('Name', 'Dist', 'Stype', 'Rad', 'Teff', 'Mass', 'RA', 'Dec', 'WDSsep', 'WDSdmag'),
                          dtype = ('S32', 'd', 'c', 'd', 'd', 'd', 'd', 'd', 'd', 'd'))
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
                SC_out.add_row([Name[i], Dist[i], Stype[i], Rad[i], Teff[i], Mass[i], RA[i], Dec[i], WDSsep[i], WDSdmag[i]])
        
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
