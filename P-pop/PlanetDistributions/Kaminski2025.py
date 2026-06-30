"""
# =============================================================================
# P-pop
# A Monte-Carlo tool to simulate exoplanet populations
# =============================================================================
"""


# =============================================================================
# IMPORTS
# =============================================================================

import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# SAG13
# =============================================================================

class PlanetDistribution():
    """
    https://ui.adsabs.harvard.edu/abs/2025A%26A...696A.101K/abstract
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
        print('--> Initializing Kaminski2025 planet distribution')
        
        # Model parameters.
        self.returns = ['Mp', 'Porb']
        if (Scenario == 'baseline'):
            self.Rates = np.array([[0.72, 0.63], [0.11, 0.06]])
        elif (Scenario == 'pessimistic'):
            self.Rates = np.array([[0.49, 0.36], [0.04, 0.02]])
        elif (Scenario == 'optimistic'):
            self.Rates = np.array([[1.01, 1.01], [0.22, 0.14]])
        elif (Scenario == 'mc'):
            raise UserWarning('Scenario mc is not implemented yet')
        else:
            print('--> WARNING: '+str(Scenario)+' is an unknown scenario')
            Scenario = 'baseline'
            self.Rates = np.array([[0.72, 0.63], [0.11, 0.06]])
        print('--> Using scenario '+str(Scenario))
        msini2m = 1./np.sin(np.deg2rad(60.))
        self.BinsMp = np.log(np.array([0.5, 3., 10.])*msini2m) # Mearth
        self.BinsPorb = np.log(np.array([1., 10., 100.])) # d
        
        self.F0 = np.sum(self.Rates)
        
        # Ndraw = 100000
        # Mps = []
        # Porbs = []
        # for i in range(Ndraw):
        #     Mp, Porb = self.draw()
        #     Mps += [Mp]
        #     Porbs += [Porb]
        # Mps = np.concatenate(Mps)
        # Porbs = np.concatenate(Porbs)
        # eta = len(Mps) / Ndraw
        
        # import pdb; pdb.set_trace()
        
        pass
    
    def draw(self,
             Mp_range=[0.57735027, 11.54700538], # Mearth
             Porb_range=[1., 100.], # d
             Nplanets=None,
             Scale=1.,
             Star=None):
        """
        Parameters
        ----------
        Mp_range: list
            Requested planet mass range (Mearth).
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
        Mp: array
            Mass (Mearth) of drawn planets.
        Porb: array
            Orbital period (d) of drawn planets.
        """
        
        Mp = [] # Mearth
        Porb = [] # d
        
        # Apply scaling for the planet occurrence rates.
        tempF0 = self.F0*Scale
        tempRates = self.Rates.copy()*Scale
        
        # If the number of planets is not given, draw it from a Poisson
        # distribution. Note that the final number of drawn planets might be
        # smaller than the drawn number because of clipping to the requested Rp
        # and Porb range.
        if (Nplanets is None):
            Nplanets = np.random.poisson(tempF0)
            for i in range(Nplanets):
                
                # Randomly select bin.
                temp = np.random.choice(len(tempRates.flatten()), p=tempRates.flatten()/tempF0)
                ww = np.unravel_index(temp, tempRates.shape)
                tempMp = np.exp(self.BinsMp[ww[0]]+(self.BinsMp[ww[0]+1]-self.BinsMp[ww[0]])*np.random.rand()) # Mearth
                tempPorb = np.exp(self.BinsPorb[ww[1]]+(self.BinsPorb[ww[1]+1]-self.BinsPorb[ww[1]])*np.random.rand()) # d
                if (Mp_range[0] <= tempMp <= Mp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Mp += [tempMp] # Mearth
                    Porb += [tempPorb] # d
        
        # If the number of planets is given, draw exactly this number of
        # planets in the requested Rp and Porb range.
        else:
            while (len(Mp) < Nplanets):
                
                # Randomly select bin.
                temp = np.random.choice(len(tempRates.flatten()), p=tempRates.flatten()/tempF0)
                ww = np.unravel_index(temp, tempRates.shape)
                tempMp = np.exp(self.BinsMp[ww[0]]+(self.BinsMp[ww[0]+1]-self.BinsMp[ww[0]])*np.random.rand()) # Mearth
                tempPorb = np.exp(self.BinsPorb[ww[1]]+(self.BinsPorb[ww[1]+1]-self.BinsPorb[ww[1]])*np.random.rand()) # d
                if (Mp_range[0] <= tempMp <= Mp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Mp += [tempMp] # Mearth
                    Porb += [tempPorb] # d
        
        return np.array(Mp), np.array(Porb)
    
    def SummaryPlot(self,
                    Ntest=100000,
                    Mp_range=[0.57735027, 11.54700538], # Mearth
                    Porb_range=[1., 100.], # d
                    FigDir=None,
                    block=True):
        """
        Parameters
        ----------
        Ntest: int
            Number of test draws for summary plot.
        Mp_range: list
            Requested planet mass range (Mearth).
        Porb_range: list
            Requested planet orbital period range (d).
        FigDir: str
            Directory to which summary plots are saved.
        block: bool
            If True, blocks plots when showing.
        """
        
        Mp = []
        Porb = []
        for i in range(Ntest):
            tempMp, tempPorb = self.draw(Mp_range,
                                         Porb_range)
            Mp += [tempMp]
            Porb += [tempPorb]
        Mp = np.concatenate(Mp)
        Porb = np.concatenate(Porb)
        
        print('--> Kaminski2025:\n%.2f planets/star in [%.1f, %.1f] Mearth and [%.1f, %.1f] d' % (len(Mp)/float(Ntest), Mp_range[0], Mp_range[1], Porb_range[0], Porb_range[1]))
        
        Weight = 1./len(Mp)
        f, ax = plt.subplots(1, 2)
        ax[0].hist(Mp, bins=np.logspace(np.log10(np.min(Mp)), np.log10(np.max(Mp)), 25), weights=np.ones_like(Mp)*Weight)
        ax[0].set_xscale('log')
        ax[0].grid(axis='y')
        ax[0].set_xlabel('Planet mass [$M_\oplus$]')
        ax[0].set_ylabel('Fraction')
        ax[1].hist(Porb, bins=np.logspace(np.log10(np.min(Porb)), np.log10(np.max(Porb)), 25), weights=np.ones_like(Porb)*Weight)
        ax[1].set_xscale('log')
        ax[1].grid(axis='y')
        ax[1].set_xlabel('Planet orbital period [d]')
        ax[1].set_ylabel('Fraction')
        plt.suptitle('Kaminski2025')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'PlanetDistribution_Kaminski2025.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
