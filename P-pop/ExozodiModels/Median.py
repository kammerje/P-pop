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
import scipy.stats as stats


# =============================================================================
# MEDIAN
# =============================================================================

class ExozodiModel():
    """
    https://ui.adsabs.harvard.edu/abs/2020AJ....159..177E/abstract
    """
    
    def __init__(self,
                 Scenario):
        """
        Parameters
        ----------
        Scenario: 'baseline', 'pessimistic', 'optimistic', 'mc'
            Scenario for exozodi level.
        """
        
        # Print.
        print('--> Initializing Median exozodi distribution')
        
        # Model parameters.
        self.zodi = 3.
        
        pass
    
    def getExozodiLevel(self):
        """
        Returns
        -------
        z: float
            Exozodiacal dust level of drawn system.
        """
        
        z = self.zodi
        
        return z
    
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
        z = []
        for i in range(Ntest):
            z += [self.getExozodiLevel()]
        z = np.array(z)
        
        print('--> Median:\n%.2f/%.2f median/mean system exozodiacal dust level' % (np.median(z), np.mean(z)))
        
        Colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        f, ax = plt.subplots(1, 1)
        Weight = 1./len(z)
        y, x, _ = ax.hist(z, bins=np.logspace(-5, 5, 50), weights=np.ones_like(z)*Weight, color=Colors[0], alpha=0.5, label='Median')
        ax.set_xscale('log')
        ax.grid(axis='y')
        ax.set_xlabel('System exozodiacal dust level')
        ax.set_ylabel('Fraction')
        ax.legend()
        plt.suptitle('Median')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'ExozodiModel.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
