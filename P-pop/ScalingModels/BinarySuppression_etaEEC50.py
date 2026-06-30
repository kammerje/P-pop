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
# BINARYSUPPRESSION
# =============================================================================

class ScalingModel():
    
    def __init__(self):
        
        pass
    
    def getScale(self,
                 Star):
        """
        Parameters
        ----------
        Star: instance
            Instance of class Star.
        
        Returns
        -------
        Scale: float
            Scaling factor for the planet occurrence rates.
        """
        
        fact = 0.500/0.408
        try:
            if (Star.WDSsep >= 50.):
                return 1. * fact
            else:
                return 0.3 * fact
        except:
            return 1. * fact
    
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
        
        pass
