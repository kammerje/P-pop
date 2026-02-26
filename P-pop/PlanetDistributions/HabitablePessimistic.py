from __future__ import division
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

from mpl_toolkits.mplot3d import Axes3D

import Star


# =============================================================================
# HABITABLEPESSIMISTIC
# =============================================================================

class PlanetDistribution():
    """
    eta = 0.15 (AFGK)
    eta = 0.25 (M)
    """
    
    def __init__(self,
                 Scenario):
        """
        Parameters
        ----------
        Scenario: 'baseline', 'pessimistic', 'optimistic', 'mc'
            Scenario for planet occurrence rates.
        """
        
        # Print.
        print('--> Initializing HabitablePessimistic planet distribution')
        
        # Constants.
        self.G = 6.674e-11 # m^3/kg/s^2
        self.Rsun = 695700000. # m
        self.Msun = 1.989e30 # kg
        self.sigma = 5.670e-8 # W/m^2/K^4
        
        # Model parameters.
        self.returns = ['Rp', 'Porb']
        self.Rates = {'A': 0.15, 'F': 0.15, 'G': 0.15, 'K': 0.15, 'M': 0.25}
        self.BinsRp = np.log(np.array([0.5, 1.5])) # Rearth
        self.BinsFp = np.log(np.array([0.35, 1.75])) # Searth
        
        pass
    
    def draw(self,
             Rp_range=None, # Rearth
             Porb_range=None, # d
             Nplanets=None,
             Scale=1.,
             Star=None):
        """
        Parameters
        ----------
        Rp_range: list
            Requested planet radius range (Rearth).
        Porb_range: list
            Requested planet orbital period range (d).
        Nplanets: None, int
            Number of planets to be drawn.
        Scale: float
            Scaling factor for the planet occurrence rates.
        Star: instance
            Instance of class Star.
        
        Returns
        -------
        Rp: array
            Radius (Rearth) of drawn planets.
        Porb: array
            Orbital period (d) of drawn planets.
        """
        
        if (Star is None):
            raise UserWarning()
        
        Rp = [] # Rearth
        Porb = [] # d
        
        # Get the planet occurrence rates.
        tempF0 = self.Rates[Star.Stype]
        
        # Apply scaling for the planet occurrence rates.
        tempF0 *= Scale
        
        # If the number of planets is not given, draw it from a Poisson
        # distribution. Note that the final number of drawn planets might be
        # smaller than the drawn number because of clipping to the requested Rp
        # and Porb range.
        if (Nplanets is None):
            Nplanets = np.random.poisson(tempF0)
            for i in range(Nplanets):
                
                # Randomly select bin.
                tempRp = np.exp(self.BinsRp[0]+(self.BinsRp[1]-self.BinsRp[0])*np.random.rand()) # Rearth
                tempFp = np.exp(self.BinsFp[0]+(self.BinsFp[1]-self.BinsFp[0])*np.random.rand()) # Searth
                
                tempPorb = self.FtoP(tempFp,
                                     Star)
                
                Rp += [tempRp] # Rearth
                Porb += [tempPorb] # d
        
        # If the number of planets is given, draw exactly this number of
        # planets in the requested Rp and Porb range.
        else:
            while (len(Rp) < Nplanets):
                
                # Randomly select bin.
                tempRp = np.exp(self.BinsRp[0]+(self.BinsRp[1]-self.BinsRp[0])*np.random.rand()) # Rearth
                tempFp = np.exp(self.BinsFp[0]+(self.BinsFp[1]-self.BinsFp[0])*np.random.rand()) # Searth
                
                tempPorb = self.FtoP(tempFp,
                                     Star)
                
                Rp += [tempRp] # Rearth
                Porb += [tempPorb] # d
        
        return np.array(Rp), np.array(Porb)
    
    def FtoP(self,
             Fp, # Searth
             Star):
        """
        Parameters
        ----------
        Fp: float
            Incident host star flux (Searth) of drawn planet.
        Star: instance
            Instance of class Star.
        
        Returns
        -------
        Porb: float
            Orbital period (d) of drawn planet.
        """
        
        rp = np.sqrt(self.sigma*Star.Teff**4*(Star.Rad*self.Rsun)**2/(1361.*Fp)) # m
        
        Porb = np.sqrt((4.*np.pi**2*rp**3)/(self.G*Star.Mass*self.Msun))/86400. # d
        
        return Porb
    
    def SummaryPlot(self,
                    Ntest=100000,
                    Rp_range=None, # Rearth
                    Porb_range=None, # d
                    FigDir=None,
                    block=True):
        """
        Parameters
        ----------
        Ntest: int
            Number of test draws for summary plot.
        Rp_range: list
            Requested planet radius range (Rearth).
        Porb_range: list
            Requested planet orbital period range (d).
        FigDir: str
            Directory to which summary plots are saved.
        block: bool
            If True, blocks plots when showing.
        """
        
        Rp = []
        Porb = []
        Sun = Star.Star('Sun',
                        10., # pc
                        'G',
                        1., # Rsun
                        5780., # K
                        1., # Msun
                        0., # deg
                        0.) # deg
        for i in range(Ntest):
            tempRp, tempPorb = self.draw(Rp_range,
                                         Porb_range,
                                         Star=Sun)
            Rp += [tempRp]
            Porb += [tempPorb]
        Rp = np.concatenate(Rp)
        Porb = np.concatenate(Porb)
        
        print('--> HabitablePessimistic:\n%.2f planets/star in [%.1f, %.1f] Rearth and [%.1f, %.1f] d' % (len(Rp)/float(Ntest), np.min(Rp), np.max(Rp), np.min(Porb), np.max(Porb)))
        
        Weight = 1./len(Rp)
        f, ax = plt.subplots(1, 2)
        ax[0].hist(Rp, bins=25, weights=np.ones_like(Rp)*Weight)
        ax[0].grid(axis='y')
        ax[0].set_xlabel('Planet radius [$R_\oplus$]')
        ax[0].set_ylabel('Fraction')
        ax[1].hist(Porb, bins=25, weights=np.ones_like(Porb)*Weight)
        ax[1].grid(axis='y')
        ax[1].set_xlabel('Planet orbital period [d]')
        ax[1].set_ylabel('Fraction')
        plt.suptitle('HabitablePessimistic')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'PlanetDistribution_HabitablePessimistic.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
