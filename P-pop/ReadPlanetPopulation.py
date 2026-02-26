"""
# =============================================================================
# P-POP
# A Monte-Carlo tool to simulate exoplanet populations
# =============================================================================
"""


# =============================================================================
# IMPORTS
# =============================================================================

import matplotlib.pyplot as plt
import numpy as np
import sys


# =============================================================================
# PLANET POPULATION
# =============================================================================

class PlanetPopulation():
    
    def __init__(self,
                 PathPlanetTable):
        """
        Parameters
        ----------
        PathPlanetTable: str
            Path of the planet table to be read.
        """
        
        # Print.
        print('--> Reading planet table '+PathPlanetTable)
        
        # Read the planet table.
        Table = open(PathPlanetTable, 'r')
        TableLines = Table.readlines()
        
        Nuniverse = []
        Rp = [] # Rearth
        Porb = [] # d
        Mp = [] # Mearth
        ep = []
        ip = [] # rad
        Omegap = [] # rad
        omegap = [] # rad
        thetap = [] # rad
        Abond = []
        AgeomVIS = []
        AgeomMIR = []
        z = []
        ap = [] # au
        rp = [] # au
        AngSep = [] # arcsec
        maxAngSep = [] # arcsec
        Fp = [] # Searth
        fp = []
        Tp = [] # K
        Nstar = []
        Rs = [] # Rsun
        Ms = [] # Msun
        Ts = [] # K
        Ds = [] # pc
        Stype = []
        RA = [] # deg
        Dec = [] # deg
        lGal = [] # deg
        bGal = [] # deg
        WDSsep = [] # au
        name = []
        
        Nlines = len(TableLines)
        for i, Line in enumerate(TableLines):
            
            if (i % 10000 == 0):
                sys.stdout.write('\rProcessed line %.0f of %.0f' % (i, Nlines))
                sys.stdout.flush()
            
            tempLine = Line.split('\t')
            if (i == 0):
                if ('nMC' in tempLine):
                    self.isold = True
                else:
                    self.isold = False
            if (i == 1):
                if ('Nuniverse' in tempLine):
                    self.isold = False
            if ((self.isold == True) and (i == 0)):
                # The second line (i = 1) contains the column names of the new
                # P-pop while the first line (i = 0) contains the column names
                # of the old P-pop.
                ColNuniverse = np.where(np.array(tempLine) == 'nMC')[0][0]
                ColRp = np.where(np.array(tempLine) == 'Rp')[0][0]
                ColPorb = np.where(np.array(tempLine) == 'Porb')[0][0]
                ColMp = np.where(np.array(tempLine) == 'Mp')[0][0]
                Colep = np.where(np.array(tempLine) == 'ecc')[0][0]
                Colip = np.where(np.array(tempLine) == 'inc')[0][0]
                ColOmegap = np.where(np.array(tempLine) == 'Omega')[0][0]
                Colomegap = np.where(np.array(tempLine) == 'omega')[0][0]
                Colthetap = np.where(np.array(tempLine) == 'theta')[0][0]
                ColAbond = np.where(np.array(tempLine) == 'Abond')[0][0]
                ColAgeomVIS = np.where(np.array(tempLine) == 'AgeomVIS')[0][0]
                ColAgeomMIR = np.where(np.array(tempLine) == 'AgeomMIR')[0][0]
                Colz = np.where(np.array(tempLine) == 'zodis')[0][0]
                Colap = np.where(np.array(tempLine) == 'a')[0][0]
                Colrp = np.where(np.array(tempLine) == 'rp')[0][0]
                ColAngSep = np.where(np.array(tempLine) == 'ang_sep')[0][0]
                ColmaxAngSep = np.where(np.array(tempLine) == 'ang_sep_max')[0][0]
                ColFp = np.where(np.array(tempLine) == 'Finc')[0][0]
                Colfp = np.where(np.array(tempLine) == 'f')[0][0]
                ColTp = np.where(np.array(tempLine) == 'Tp')[0][0]
                ColNstar = np.where(np.array(tempLine) == 'nstar')[0][0]
                ColRs = np.where(np.array(tempLine) == 'Rs')[0][0]
                ColMs = np.where(np.array(tempLine) == 'Ms')[0][0]
                ColTs = np.where(np.array(tempLine) == 'Ts')[0][0]
                ColDs = np.where(np.array(tempLine) == 'dist')[0][0]
                ColStype = np.where(np.array(tempLine) == 'stype')[0][0]
                ColRA = np.where(np.array(tempLine) == 'ra')[0][0]
                ColDec = np.where(np.array(tempLine) == 'dec')[0][0]
            elif ((self.isold == False) and (i == 1)):
                # The second line (i = 1) contains the column names of the new
                # P-pop while the first line (i = 0) contains the column names
                # of the old P-pop.
                ColNuniverse = np.where(np.array(tempLine) == 'Nuniverse')[0][0]
                ColRp = np.where(np.array(tempLine) == 'Rp')[0][0]
                ColPorb = np.where(np.array(tempLine) == 'Porb')[0][0]
                ColMp = np.where(np.array(tempLine) == 'Mp')[0][0]
                Colep = np.where(np.array(tempLine) == 'ep')[0][0]
                Colip = np.where(np.array(tempLine) == 'ip')[0][0]
                ColOmegap = np.where(np.array(tempLine) == 'Omegap')[0][0]
                Colomegap = np.where(np.array(tempLine) == 'omegap')[0][0]
                Colthetap = np.where(np.array(tempLine) == 'thetap')[0][0]
                ColAbond = np.where(np.array(tempLine) == 'Abond')[0][0]
                ColAgeomVIS = np.where(np.array(tempLine) == 'AgeomVIS')[0][0]
                ColAgeomMIR = np.where(np.array(tempLine) == 'AgeomMIR')[0][0]
                Colz = np.where(np.array(tempLine) == 'z')[0][0]
                Colap = np.where(np.array(tempLine) == 'ap')[0][0]
                Colrp = np.where(np.array(tempLine) == 'rp')[0][0]
                ColAngSep = np.where(np.array(tempLine) == 'AngSep')[0][0]
                ColmaxAngSep = np.where(np.array(tempLine) == 'maxAngSep')[0][0]
                ColFp = np.where(np.array(tempLine) == 'Fp')[0][0]
                Colfp = np.where(np.array(tempLine) == 'fp')[0][0]
                ColTp = np.where(np.array(tempLine) == 'Tp')[0][0]
                ColNstar = np.where(np.array(tempLine) == 'Nstar')[0][0]
                ColRs = np.where(np.array(tempLine) == 'Rs')[0][0]
                ColMs = np.where(np.array(tempLine) == 'Ms')[0][0]
                ColTs = np.where(np.array(tempLine) == 'Ts')[0][0]
                ColDs = np.where(np.array(tempLine) == 'Ds')[0][0]
                ColStype = np.where(np.array(tempLine) == 'Stype')[0][0]
                ColRA = np.where(np.array(tempLine) == 'RA')[0][0]
                ColDec = np.where(np.array(tempLine) == 'Dec')[0][0]
                CollGal = np.where(np.array(tempLine) == 'lGal')[0][0]
                ColbGal = np.where(np.array(tempLine) == 'bGal')[0][0]
                ColWDSsep = np.where(np.array(tempLine) == 'WDSsep')[0][0]
                Colname = np.where(np.array(tempLine) == 'name')[0][0]
            elif (((self.isold == True) and (i > 0)) or ((self.isold == False) and (i > 1))):
                Nuniverse += [int(tempLine[ColNuniverse])]
                Rp += [float(tempLine[ColRp])] # Rearth
                Porb += [float(tempLine[ColPorb])] # d
                Mp += [float(tempLine[ColMp])] # Mearth
                ep += [float(tempLine[Colep])]
                ip += [float(tempLine[Colip])] # rad
                Omegap += [float(tempLine[ColOmegap])] # rad
                omegap += [float(tempLine[Colomegap])] # rad
                thetap += [float(tempLine[Colthetap])] # rad
                Abond += [float(tempLine[ColAbond])]
                AgeomVIS += [float(tempLine[ColAgeomVIS])]
                AgeomMIR += [float(tempLine[ColAgeomMIR])]
                z += [float(tempLine[Colz])]
                ap += [float(tempLine[Colap])] # au
                rp += [float(tempLine[Colrp])] # au
                AngSep += [float(tempLine[ColAngSep])] # arcsec
                maxAngSep += [float(tempLine[ColmaxAngSep])] # arcsec
                Fp += [float(tempLine[ColFp])] # Searth
                fp += [float(tempLine[Colfp])]
                Tp += [float(tempLine[ColTp])] # K
                Nstar += [int(tempLine[ColNstar])]
                Rs += [float(tempLine[ColRs])] # Rsun
                Ms += [float(tempLine[ColMs])] # Msun
                Ts += [float(tempLine[ColTs])] # K
                Ds += [float(tempLine[ColDs])] # pc
                Stype += [str(tempLine[ColStype])]
                RA += [float(tempLine[ColRA])] # deg
                Dec += [float(tempLine[ColDec])] # deg
                try:
                    lGal += [float(tempLine[CollGal])] # deg
                    bGal += [float(tempLine[ColbGal])] # deg
                except:
                    lGal += [None] # deg
                    bGal += [None] # deg
                try:
                    WDSsep += [float(tempLine[ColWDSsep])] # au
                except:
                    WDSsep += [np.inf] # au
                try:
                    name += [str(tempLine[Colname][:-1])]
                except:
                    name += [None]
        sys.stdout.write('\rProcessed line %.0f of %.0f' % (Nlines, Nlines))
        sys.stdout.flush()
        print('')
        
        self.Nuniverse = np.array(Nuniverse)
        self.Rp = np.array(Rp) # Rearth
        self.Porb = np.array(Porb) # d
        self.Mp = np.array(Mp) # Mearth
        self.ep = np.array(ep)
        self.ip = np.array(ip) # rad
        self.Omegap = np.array(Omegap) # rad
        self.omegap = np.array(omegap) # rad
        self.thetap = np.array(thetap) # rad
        self.Abond = np.array(Abond)
        self.AgeomVIS = np.array(AgeomVIS)
        self.AgeomMIR = np.array(AgeomMIR)
        self.z = np.array(z)
        self.ap = np.array(ap) # au
        self.rp = np.array(rp) # au
        self.AngSep = np.array(AngSep) # arcsec
        self.maxAngSep = np.array(maxAngSep) # arcsec
        self.Fp = np.array(Fp) # Searth
        self.fp = np.array(fp)
        self.Tp = np.array(Tp) # K
        self.Nstar = np.array(Nstar)
        self.Rs = np.array(Rs) # Rsun
        self.Ms = np.array(Ms) # Msun
        self.Ts = np.array(Ts) # K
        self.Ds = np.array(Ds) # pc
        self.Stype = np.array(Stype)
        self.RA = np.array(RA) # deg
        self.Dec = np.array(Dec) # deg
        self.lGal = np.array(lGal) # deg
        self.bGal = np.array(bGal) # deg
        self.WDSsep = np.array(WDSsep) # au
        self.name = np.array(name)
        
        self.Phot = {}
        
        pass
    
    def ComputeHZ(self,
                  Model='MS'):
        """
        Parameters
        ----------
        Model: MS, POST-MS
            Model for computing the habitable zone (au).
        """
        
        # Print.
        print('--> Computing the habitable zone')
        
        # Compute the habitable zone (au).
        if (Model == 'MS'):
            S0in, S0out = 1.7665, 0.3240
            Ain, Aout = 1.3351E-4, 5.3221E-5
            Bin, Bout = 3.1515E-9, 1.4288E-9
            Cin, Cout = -3.3488E-12, -1.1049E-12
            T = self.Ts-5780. # K
            self.HZin = S0in+Ain*T+Bin*T**2+Cin*T**3 # au
            self.HZout = S0out+Aout*T+Bout*T**2+Cout*T**3 # au
            # S0in, S0out = 1.0140, 0.3438
            # Ain, Aout = 8.1774E-5, 5.8942E-5
            # Bin, Bout = 1.7063E-9, 1.6558E-9
            # Cin, Cout = -4.3241E-12, -3.0045E-12
            # Din, Dout = -6.6462e-16, -5.2982e-16
            # T = self.Ts-5780. # K
            # self.HZin = S0in+Ain*T+Bin*T**2+Cin*T**3+Din*T**4 # au
            # self.HZout = S0out+Aout*T+Bout*T**2+Cout*T**3+Dout*T**4 # au
        elif (Model == 'POST-MS'):
            S0in, S0out = 1.1066, 0.3240
            Ain, Aout = 1.2181E-4, 5.3221E-5
            Bin, Bout = 1.5340E-8, 1.4288E-9
            Cin, Cout = -1.5018E-12, -1.1049E-12
            T = self.Ts-5780. # K
            self.HZin = S0in+Ain*T+Bin*T**2+Cin*T**3 # au
            self.HZout = S0out+Aout*T+Bout*T**2+Cout*T**3 # au
        else:
            print('--> WARNING: '+str(Model)+' is an unknown model')
            Model = 'MS'
            S0in, S0out = 1.7665, 0.3240
            Ain, Aout = 1.3351E-4, 5.3221E-5
            Bin, Bout = 3.1515E-9, 1.4288E-9
            Cin, Cout = -3.3488E-12, -1.1049E-12
            T = self.Ts-5780. # K
            self.HZin = S0in+Ain*T+Bin*T**2+Cin*T**3 # au
            self.HZout = S0out+Aout*T+Bout*T**2+Cout*T**3 # au
        print('--> Using model '+str(Model))
        
        pass
    
    def appendPhotometry(self,
                         PathPhotometryTable,
                         Tag):
        """
        Parameters
        ----------
        PathPhotometryTable: str
            Path of the photometry table to be read.
        Tag: str
            Tag for the self.Phot dictionary under which the photometry is
            appended.
        """
        
        # Print.
        print('--> Reading photometry table '+PathPhotometryTable)
        
        # Read the photometry table.
        Table = open(PathPhotometryTable, 'r')
        TableLines = Table.readlines()
        
        self.Phot[Tag] = {}
        self.Phot[Tag]['HEAD'] = []
        self.Phot[Tag]['DATA'] = []
        
        Nlines = len(TableLines)
        for i, Line in enumerate(TableLines):
            
            if (i % 10000 == 0):
                sys.stdout.write('\rProcessed line %.0f of %.0f' % (i, Nlines))
                sys.stdout.flush()
            
            tempLine = Line.split('\t')
            if ((self.isold == True) and (i == 0)):
                # The second line (i = 1) contains the column names of the new
                # P-pop while the first line (i = 0) contains the column names
                # of the old P-pop.
                for j in range(len(tempLine)-1):
                    self.Phot[Tag]['HEAD'] += [tempLine[j]]
                    self.Phot[Tag]['DATA'] += [[]]
            elif ((self.isold == False) and (i == 1)):
                # The second line (i = 1) contains the column names of the new
                # P-pop while the first line (i = 0) contains the column names
                # of the old P-pop.
                for j in range(len(tempLine)-1):
                    self.Phot[Tag]['HEAD'] += [tempLine[j]]
                    self.Phot[Tag]['DATA'] += [[]]
            elif (((self.isold == True) and (i > 0)) or ((self.isold == False) and (i > 1))):
                for j in range(len(tempLine)-1):
                    self.Phot[Tag]['DATA'][j] += [float(tempLine[j])]
        sys.stdout.write('\rProcessed line %.0f of %.0f' % (Nlines, Nlines))
        sys.stdout.flush()
        print('')
        
        for i in range(len(self.Phot[Tag]['DATA'])):
            self.Phot[Tag]['DATA'][i] = np.array(self.Phot[Tag]['DATA'][i])
        
        pass
