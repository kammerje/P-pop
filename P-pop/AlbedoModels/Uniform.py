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
# UNIFORM
# =============================================================================

class AlbedoModel():
    
    def __init__(self):
        
        pass
    
    def getAbond(self,
                 Rp): # Rearth
        """
        Parameters
        ----------
        Rp: array
            Radius (Rearth) of drawn planets.
        
        Returns
        -------
        Abond: array
            Bond albedo of drawn planets.
        """
        
        return 0.8*np.random.rand(len(Rp))
    
    def getAgeomVIS(self,
                    Rp): # Rearth
        """
        Parameters
        ----------
        Rp: array
            Radius (Rearth) of drawn planets.
        
        Returns
        -------
        AgeomVIS: array
            Geometric albedo in the visible of drawn planets.
        """
        
        return 0.6*np.random.rand(len(Rp))
    
    def getAgeomMIR(self,
                    Rp): # Rearth
        """
        Parameters
        ----------
        Rp: array
            Radius (Rearth) of drawn planets.
        
        Returns
        -------
        AgeomMIR: array
            Geometric albedo in the mid-infrared of drawn planets.
        """
        
        return 0.1*np.random.rand(len(Rp))
    
    def SummaryPlot(self,
                    Ntest=100000,
                    FigDir=None,
                    block=True):
        """
        Parameters
        ----------
        Ntest: int
            Number of test draws for summary plot.
        FigDir: str
            Directory to which summary plots are saved.
        block: bool
            If True, blocks plots when showing.
        """
        
        Abond = self.getAbond(np.zeros(Ntest))
        AgeomVIS = self.getAgeomVIS(np.zeros(Ntest))
        AgeomMIR = self.getAgeomMIR(np.zeros(Ntest))
        
        print('--> Uniform:\n%.2f+-%.2f planet Bond albedo\n%.2f+-%.2f planet geometric albedo (VIS)\n%.2f+-%.2f planet geometric albedo (MIR)' % (np.mean(Abond), np.std(Abond), np.mean(AgeomVIS), np.std(AgeomVIS), np.mean(AgeomMIR), np.std(AgeomMIR)))
        
        Weight = 1./len(Abond)
        f, ax = plt.subplots(1, 3)
        ax[0].hist(Abond, bins=25, weights=np.ones_like(Abond)*Weight)
        ax[0].grid(axis='y')
        ax[0].set_xlabel('Bond albedo')
        ax[0].set_ylabel('Fraction')
        ax[1].hist(AgeomVIS, bins=25, weights=np.ones_like(AgeomVIS)*Weight)
        ax[1].grid(axis='y')
        ax[1].set_xlabel('Geom. albedo (VIS)')
        ax[1].set_ylabel('Fraction')
        ax[2].hist(AgeomMIR, bins=25, weights=np.ones_like(AgeomMIR)*Weight)
        ax[2].grid(axis='y')
        ax[2].set_xlabel('Geom. albedo (MIR)')
        ax[2].set_ylabel('Fraction')
        plt.suptitle('Uniform')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'AlbedoModel.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
