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
# ERTEL2018
# =============================================================================

class ExozodiModel():
    """
    https://ui.adsabs.harvard.edu/abs/2018AJ....155..194E/abstract
    """
    
    def __init__(self,
                 Scenario):
        """
        Parameters
        ----------
        Scenario: 'baseline', 'pessimistic', 'optimistic', 'mc'
            Scenario for exozodi level.
        """
        
        pass
    
    def getExozodiLevel(self):
        """
        Returns
        -------
        z: float
            Exozodiacal dust level of drawn system.
        """
        
        # Clean Sun-like stars.
        mu = np.log(26.); zeta = 1.2
        
        # Clean stars.
#        mu = np.log(13.); zeta = 1.5
        
        return np.random.lognormal(mean=mu, sigma=zeta)
    
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
        
        z = []
        for i in range(Ntest):
            z += [self.getExozodiLevel()]
        z = np.array(z)
        
        print('--> Ertel2018:\n%.2f/%.2f median/mean system exozodiacal dust level' % (np.median(z), np.mean(z)))
        
        Weight = 1./len(z)
        f, ax = plt.subplots(1, 1)
        ax.hist(z, bins=np.logspace(np.log10(np.min(z)), np.log10(np.max(z)), 25), weights=np.ones_like(z)*Weight)
        ax.set_xscale('log')
        ax.grid(axis='y')
        ax.set_xlabel('System exozodiacal dust level')
        ax.set_ylabel('Fraction')
        plt.suptitle('Ertel2018')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'ExozodiModel.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
