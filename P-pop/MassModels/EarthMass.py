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
import sys

sys.path.append('MassModels/Forecaster')
import mr_forecast


# =============================================================================
# CHEN2017
# =============================================================================

class MassModel():
    """
    1 Earth mass/Earth radius for all planets.
    """
    
    def __init__(self):
        
        pass
    
    def RadiusToMass(self,
                     Rp): # Rearth
        """
        Parameters
        ----------
        Rp: array
            Radius (Rearth) of drawn planets.
        
        Returns
        -------
        Mp: array
            Mass (Mearth) of drawn planets.
        """
        
        # If there are no planets in the system return an empty array.
        if (len(Rp) == 0):
            return np.array([])
        else:
            
            # Forecast planet mass.
            Mp = np.array([1.]*len(Rp)) # Mearth
            return Mp
    
    def MassToRadius(self,
                     Mp): # Mearth
        """
        Parameters
        ----------
        Mp: array
            Mass (Mearth) of drawn planets.
        
        Returns
        -------
        Rp: array
            Radius (Rearth) of drawn planets.
        """
        
        # If there are no planets in the system return an empty array.
        if (len(Mp) == 0):
            return np.array([])
        else:
            
            # Forecast planet radius.
            Rp = np.array([1.]*len(Mp)) # Rearth
            return Rp
    
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
        
        Ntest = Ntest//10
        Rp = np.ones((Ntest))
        Mp = self.RadiusToMass(Rp)
        
        print('--> EarthMass:\n%.2f Rearth mean input planet radius\n%.2f Mearth mean output planet mass' % (np.mean(Rp), np.mean(Mp)))
        
        Weight = 1./len(Rp)
        f, ax = plt.subplots(1, 2)
        ax[0].hist(Rp, bins=25, weights=np.ones_like(Rp)*Weight)
        ax[0].grid(axis='y')
        ax[0].set_xlabel('Planet radius [$R_\oplus$]')
        ax[0].set_ylabel('Fraction')
        ax[1].hist(Mp, bins=25, weights=np.ones_like(Mp)*Weight)
        ax[1].grid(axis='y')
        ax[1].set_xlabel('Planet mass [$M_\oplus$]')
        ax[1].set_ylabel('Fraction')
        plt.suptitle('Chen2017')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'MassModel.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
