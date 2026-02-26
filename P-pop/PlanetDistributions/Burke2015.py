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


# =============================================================================
# BURKE2015
# =============================================================================

class PlanetDistribution():
    """
    https://ui.adsabs.harvard.edu/abs/2015ApJ...809....8B/abstract
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
        print('--> Initializing Burke2015 planet distribution')
        
        # Model parameters.
        self.returns = ['Rp', 'Porb']
        if (Scenario == 'baseline'):
            self.F0 = 0.731
            self.alpha1 = 19.684
            self.alpha2 = -1.779
            self.Rbrk = 0.941 # Rearth
            self.beta = -0.655
        elif (Scenario == 'pessimistic'):
            self.F0 = 0.502
            self.alpha1 = 27.317
            self.alpha2 = -1.334
            self.Rbrk = 0.940 # Rearth
            self.beta = -0.773
        elif (Scenario == 'optimistic'):
            self.F0 = 1.121
            self.alpha1 = 29.886
            self.alpha2 = -2.183
            self.Rbrk = 0.940 # Rearth
            self.beta = -0.576
        elif (Scenario == 'mc'):
            raise UserWarning('Scenario mc is not implemented yet')
        else:
            print('--> WARNING: '+str(Scenario)+' is an unknown scenario')
            Scenario = 'baseline'
            self.F0 = 0.731
            self.alpha1 = 19.684
            self.alpha2 = -1.779
            self.Rbrk = 0.941 # Rearth
            self.beta = -0.655
        print('--> Using scenario '+str(Scenario))
        self.Rp_lims = [0.75, 2.5] # Rearth
        self.Porb_lims = [50., 300.] # d
        self.R0 = (self.Rp_lims[0]+self.Rp_lims[1])/2. # Rearth
        self.P0 = (self.Porb_lims[0]+self.Porb_lims[1])/2. # d
        self.CR = 1./(((self.Rbrk**(self.alpha1+1.)-self.Rp_lims[0]**(self.alpha1+1.))/self.R0**self.alpha1/(self.alpha1+1.))+((self.Rp_lims[1]**(self.alpha2+1.)-self.Rbrk**(self.alpha2+1.))/self.R0**self.alpha1/(self.alpha2+1.)*self.Rbrk**(self.alpha1-self.alpha2)))
        self.CP = 1./((self.Porb_lims[1]**(self.beta+1.)-self.Porb_lims[0]**(self.beta+1.))/self.P0**self.beta/(self.beta+1.))
        
        pass
    
    def iCDF_R(self,
               x):
        """
        Parameters
        ----------
        x: float, 0 <= x <= 1
            Uniformly distributed random number.
        
        Returns
        -------
        Rp: float
            Planet radius (Rearth) distributed according to Burke2015 planet
            distribution.
        """
        
        # Check whether Rp < Rbrk or Rp > Rbrk.
        temp = self.CR*(self.Rbrk**(self.alpha1+1.)-self.Rp_lims[0]**(self.alpha1+1.))/self.R0**self.alpha1/(self.alpha1+1.)
        if (x < temp):
            return (((self.alpha1+1.)*self.R0**self.alpha1/self.CR)*x+self.Rp_lims[0]**(self.alpha1+1.))**(1./(self.alpha1+1.))
        else:
            return (((self.alpha2+1.)*self.R0**self.alpha1*self.Rbrk**(self.alpha2-self.alpha1)/self.CR)*(x-temp)+self.Rbrk**(self.alpha2+1.))**(1./(self.alpha2+1.))
    
    def iCDF_P(self,
               x):
        """
        Parameters
        ----------
        x: float, 0 <= x <= 1
            Uniformly distributed random number.
        
        Returns
        -------
        Porb: float
            Planet orbital period (d) distributed according to Burke2015 planet
            distribution.
        """
        
        return (((self.beta+1.)*self.P0**self.beta/self.CP)*x+self.Porb_lims[0]**(self.beta+1.))**(1./(self.beta+1.))
    
    def draw(self,
             Rp_range=[0.75, 2.5], # Rearth
             Porb_range=[50., 300.], # d
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
        
        Rp = [] # Rearth
        Porb = [] # d
        
        # Apply scaling for the planet occurrence rates.
        tempF0 = self.F0*Scale
        
        # If the number of planets is not given, draw it from a Poisson
        # distribution. Note that the final number of drawn planets might be
        # smaller than the drawn number because of clipping to the requested Rp
        # and Porb range.
        if (Nplanets is None):
            Nplanets = np.random.poisson(tempF0)
            for i in range(Nplanets):
                tempRp = self.iCDF_R(np.random.rand()) # Rearth
                tempPorb = self.iCDF_P(np.random.rand()) # d
                if (Rp_range[0] <= tempRp <= Rp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Rp += [tempRp] # Rearth
                    Porb += [tempPorb] # d
        
        # If the number of planets is given, draw exactly this number of
        # planets in the requested Rp and Porb range.
        else:
            while (len(Rp) < Nplanets):
                tempRp = self.iCDF_R(np.random.rand()) # Rearth
                tempPorb = self.iCDF_P(np.random.rand()) # d
                if (Rp_range[0] <= tempRp <= Rp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Rp += [tempRp] # Rearth
                    Porb += [tempPorb] # d
        
        return np.array(Rp), np.array(Porb)
    
    def SummaryPlot(self,
                    Ntest = 100000,
                    Rp_range=[0.75, 2.5], # Rearth
                    Porb_range=[50., 300.], # d
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
        for i in range(Ntest):
            tempRp, tempPorb = self.draw(Rp_range,
                                         Porb_range)
            Rp += [tempRp]
            Porb += [tempPorb]
        Rp = np.concatenate(Rp)
        Porb = np.concatenate(Porb)
        
        print('--> Burke2015:\n%.2f planets/star in [%.1f, %.1f] Rearth and [%.1f, %.1f] d' % (len(Rp)/float(Ntest), Rp_range[0], Rp_range[1], Porb_range[0], Porb_range[1]))
        
        Weight = 1./len(Rp)
        f, ax = plt.subplots(1, 2)
        ax[0].hist(Rp, bins=np.logspace(np.log10(np.min(Rp)), np.log10(np.max(Rp)), 25), weights=np.ones_like(Rp)*Weight)
        ax[0].set_xscale('log')
        ax[0].grid(axis='y')
        ax[0].set_xlabel('Planet radius [$R_\oplus$]')
        ax[0].set_ylabel('Fraction')
        ax[1].hist(Porb, bins=np.logspace(np.log10(np.min(Porb)), np.log10(np.max(Porb)), 25), weights=np.ones_like(Porb)*Weight)
        ax[1].set_xscale('log')
        ax[1].grid(axis='y')
        ax[1].set_xlabel('Planet orbital period [d]')
        ax[1].set_ylabel('Fraction')
        plt.suptitle('Burke2015')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'PlanetDistribution_Burke2015.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
