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
# CIRCULAR
# =============================================================================

class EccentricityModel():
    
    def __init__(self):
        
        pass
    
    def getEccentricity(self,
                        Porb): # d
        """
        Parameters
        ----------
        Porb: array
            Orbital period (d) of drawn planets.
        
        Returns
        -------
        ep: array
            Eccentricity of drawn planets.
        """
        
        return np.zeros_like(Porb)
    
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
        
        ep = self.getEccentricity(np.zeros(Ntest))
        
        print('--> Circular:\n%.2f+-%.2f planet eccentricity' % (np.mean(ep), np.std(ep)))
        
        Weight = 1./len(ep)
        f, ax = plt.subplots(1, 1)
        ax.hist(ep, bins=25, weights=np.ones_like(ep)*Weight)
        ax.grid(axis='y')
        ax.set_xlabel('Planet eccentricity')
        ax.set_ylabel('Fraction')
        plt.suptitle('Circular')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'EccentricityModel.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
