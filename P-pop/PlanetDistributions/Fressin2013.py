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


# =============================================================================
# FRESSIN2013
# =============================================================================

class PlanetDistribution():
    """
    https://ui.adsabs.harvard.edu/abs/2013ApJ...766...81F/abstract
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
        print('--> Initializing Fressin2013 planet distribution')
        
        # Model parameters.
        self.returns = ['Rp', 'Porb']
        if (Scenario == 'baseline'):
            self.Rates = np.array([[0.18, 0.61, 1.72, 2.7, 2.7, 2.93, 4.08, 3.46, 0., 0., 0.], [0.17, 0.74, 1.49, 2.9, 4.3, 4.49, 5.29, 3.66, 6.54, 0., 0.], [0.035, 0.18, 0.73, 1.93, 3.67, 5.29, 6.45, 5.25, 4.31, 3.09, 0.], [0.004, 0.006, 0.11, 0.091, 0.29, 0.32, 0.49, 0.66, 0.43, 0.53, 0.24], [0.015, 0.067, 0.17, 0.18, 0.27, 0.23, 0.35, 0.71, 1.25, 0.94, 1.05]])*1e-2
        elif (Scenario == 'pessimistic'):
            self.Rates = np.array([[0.18, 0.61, 1.72, 2.7, 2.7, 2.93, 4.08, 3.46, 0., 0., 0.], [0.17, 0.74, 1.49, 2.9, 4.3, 4.49, 5.29, 3.66, 6.54, 0., 0.], [0.035, 0.18, 0.73, 1.93, 3.67, 5.29, 6.45, 5.25, 4.31, 3.09, 0.], [0.004, 0.006, 0.11, 0.091, 0.29, 0.32, 0.49, 0.66, 0.43, 0.53, 0.24], [0.015, 0.067, 0.17, 0.18, 0.27, 0.23, 0.35, 0.71, 1.25, 0.94, 1.05]])*1e-2
            self.Rates -= np.array([[0.04, 0.15, 0.43, 0.6, 0.83, 1.05, 1.88, 2.81, 0., 0., 0.], [0.03, 0.13, 0.23, 0.56, 0.73, 1., 1.48, 1.21, 2.2, 0., 0.], [0.011, 0.03, 0.09, 0.19, 0.39, 0.64, 1.01, 1.05, 1.03, 0.9, 0.], [0.003, 0.006, 0.03, 0.03, 0.07, 0.08, 0.12, 0.16, 0.17, 0.21, 0.15], [0.007, 0.018, 0.03, 0.04, 0.06, 0.06, 0.1, 0.17, 0.29, 0.28, 0.3]])*1e-2
        elif (Scenario == 'optimistic'):
            self.Rates = np.array([[0.18, 0.61, 1.72, 2.7, 2.7, 2.93, 4.08, 3.46, 0., 0., 0.], [0.17, 0.74, 1.49, 2.9, 4.3, 4.49, 5.29, 3.66, 6.54, 0., 0.], [0.035, 0.18, 0.73, 1.93, 3.67, 5.29, 6.45, 5.25, 4.31, 3.09, 0.], [0.004, 0.006, 0.11, 0.091, 0.29, 0.32, 0.49, 0.66, 0.43, 0.53, 0.24], [0.015, 0.067, 0.17, 0.18, 0.27, 0.23, 0.35, 0.71, 1.25, 0.94, 1.05]])*1e-2
            self.Rates += np.array([[0.04, 0.15, 0.43, 0.6, 0.83, 1.05, 1.88, 2.81, 0., 0., 0.], [0.03, 0.13, 0.23, 0.56, 0.73, 1., 1.48, 1.21, 2.2, 0., 0.], [0.011, 0.03, 0.09, 0.19, 0.39, 0.64, 1.01, 1.05, 1.03, 0.9, 0.], [0.003, 0.006, 0.03, 0.03, 0.07, 0.08, 0.12, 0.16, 0.17, 0.21, 0.15], [0.007, 0.018, 0.03, 0.04, 0.06, 0.06, 0.1, 0.17, 0.29, 0.28, 0.3]])*1e-2
        elif (Scenario == 'mc'):
            raise UserWarning('Scenario mc is not implemented yet')
        else:
            print('--> WARNING: '+str(Scenario)+' is an unknown scenario')
            Scenario = 'baseline'
            self.Rates = np.array([[0.18, 0.61, 1.72, 2.7, 2.7, 2.93, 4.08, 3.46, 0., 0., 0.], [0.17, 0.74, 1.49, 2.9, 4.3, 4.49, 5.29, 3.66, 6.54, 0., 0.], [0.035, 0.18, 0.73, 1.93, 3.67, 5.29, 6.45, 5.25, 4.31, 3.09, 0.], [0.004, 0.006, 0.11, 0.091, 0.29, 0.32, 0.49, 0.66, 0.43, 0.53, 0.24], [0.015, 0.067, 0.17, 0.18, 0.27, 0.23, 0.35, 0.71, 1.25, 0.94, 1.05]])*1e-2
        print('--> Using scenario '+str(Scenario))
        self.BinsRp = np.log(np.array([0.8, 1.25, 2., 4., 6., 22.])) # Rearth
        self.BinsPorb = np.log(np.array([0.8, 2., 3.4, 5.9, 10., 17., 29., 50., 85., 145., 245., 418.])) # d
        
        # Extrapolate empty bins.
        self.Rates_extrap = self.Extrapolate()
        self.F0 = np.sum(self.Rates_extrap)
        
        pass
    
    def Extrapolate(self):
        """
        Returns
        -------
        Rates_extrap: array
            Planet occurrence rates extrapolated with a bivariate quadratic
            function of the form
            a*Rp^2 + b*Porb^2 + c*Rp + d*Porb + e*Rp*Porb + f.
        """
        
        # Print.
        print('--> Extrapolating empty bins with bivariate quadratic function')
        
        x = (self.BinsPorb[:-1]+self.BinsPorb[1:])/2.
        y = (self.BinsRp[:-1]+self.BinsRp[1:])/2.
        self.xx, self.yy = np.meshgrid(x, y)
        self.zz = np.log(self.Rates)
        
        WW = np.isinf(self.zz) == False
        XX = self.xx[WW]
        YY = self.yy[WW]
        ZZ = self.zz[WW]
        
        # Compute best fit bivariate quadratic function.
        A_bqf = np.array([XX**2, YY**2, XX, YY, XX*YY, XX*0.+1.]).T
        p_bqf = np.linalg.lstsq(A_bqf, ZZ, rcond=None)
        def bqf(p, xx, yy):
            return p[0]*xx**2+p[1]*yy**2+p[2]*xx+p[3]*yy+p[4]*xx*yy+p[5]*(xx*0.+1.)
        self.r_bqf = bqf(p_bqf[0], self.xx.flatten(), self.yy.flatten()).reshape(self.zz.shape)
        
        # Extrapolate empty bins.
        ww = np.isinf(self.zz) == True
        Rates_extrap = self.Rates.copy()
        Rates_extrap[ww] = np.exp(bqf(p_bqf[0], self.xx[ww], self.yy[ww]))
        Rates_extrap[Rates_extrap < 0.] = 0.
        
        return Rates_extrap
    
    def draw(self,
             Rp_range=[0.8, 22.], # Rearth
             Porb_range=[0.8, 418.], # d
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
        tempRates_extrap = self.Rates_extrap.copy()*Scale
        
        # If the number of planets is not given, draw it from a Poisson
        # distribution. Note that the final number of drawn planets might be
        # smaller than the drawn number because of clipping to the requested Rp
        # and Porb range.
        if (Nplanets is None):
            Nplanets = np.random.poisson(tempF0)
            for i in range(Nplanets):
                
                # Randomly select bin.
                temp = np.random.choice(len(tempRates_extrap.flatten()), p=tempRates_extrap.flatten()/tempF0)
                ww = np.unravel_index(temp, tempRates_extrap.shape)
                tempRp = np.exp(self.BinsRp[ww[0]]+(self.BinsRp[ww[0]+1]-self.BinsRp[ww[0]])*np.random.rand()) # Rearth
                tempPorb = np.exp(self.BinsPorb[ww[1]]+(self.BinsPorb[ww[1]+1]-self.BinsPorb[ww[1]])*np.random.rand()) # d
                if (Rp_range[0] <= tempRp <= Rp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Rp += [tempRp] # Rearth
                    Porb += [tempPorb] # d
        
        # If the number of planets is given, draw exactly this number of
        # planets in the requested Rp and Porb range.
        else:
            while (len(Rp) < Nplanets):
                
                # Randomly select bin.
                temp = np.random.choice(len(tempRates_extrap.flatten()), p=tempRates_extrap.flatten()/tempF0)
                ww = np.unravel_index(temp, tempRates_extrap.shape)
                tempRp = np.exp(self.BinsRp[ww[0]]+(self.BinsRp[ww[0]+1]-self.BinsRp[ww[0]])*np.random.rand()) # Rearth
                tempPorb = np.exp(self.BinsPorb[ww[1]]+(self.BinsPorb[ww[1]+1]-self.BinsPorb[ww[1]])*np.random.rand()) # d
                if (Rp_range[0] <= tempRp <= Rp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Rp += [tempRp] # Rearth
                    Porb += [tempPorb] # d
        
        return np.array(Rp), np.array(Porb)
    
    def SummaryPlot(self,
                    Ntest=100000,
                    Rp_range=[0.8, 22.], # Rearth
                    Porb_range=[0.8, 418.], # d
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
        
        Rp_extrap = []
        Porb_extrap = []
        for i in range(Ntest):
            tempRp, tempPorb = self.draw(Rp_range,
                                         Porb_range)
            Rp_extrap += [tempRp]
            Porb_extrap += [tempPorb]
        Rp_extrap = np.concatenate(Rp_extrap)
        Porb_extrap = np.concatenate(Porb_extrap)
        dummy = self.Rates_extrap.copy()
        self.Rates_extrap = self.Rates.copy()
        self.F0 = np.sum(self.Rates_extrap)
        Rp = []
        Porb = []
        for i in range(Ntest):
            tempRp, tempPorb = self.draw(Rp_range,
                                         Porb_range)
            Rp += [tempRp]
            Porb += [tempPorb]
        Rp = np.concatenate(Rp)
        Porb = np.concatenate(Porb)
        self.Rates_extrap = dummy
        self.F0 = np.sum(self.Rates_extrap)
        
        print('--> Fressin2013:\n%.2f planets/star in [%.1f, %.1f] Rearth and [%.1f, %.1f] d' % (len(Rp_extrap)/float(Ntest), Rp_range[0], Rp_range[1], Porb_range[0], Porb_range[1]))
        
        Colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        Weight = 1./len(Rp_extrap)
        f, ax = plt.subplots(2, 2)
        temp = np.hstack((Rp_extrap, Rp))
        ax[0, 0].hist(Rp_extrap, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), weights=np.ones_like(Rp_extrap)*Weight, color=Colors[0], alpha=0.5, label='Extrap.')
        ax[0, 0].hist(Rp, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), weights=np.ones_like(Rp)*Weight, color=Colors[1], alpha=0.5, label='Original')
        ax[0, 0].set_xscale('log')
        ax[0, 0].grid(axis='y')
        ax[0, 0].set_xlabel('Planet radius [$R_\oplus$]')
        ax[0, 0].set_ylabel('Fraction')
        ax[0, 0].legend()
        temp = np.hstack((Porb_extrap, Porb))
        ax[0, 1].hist(Porb_extrap, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), weights=np.ones_like(Porb_extrap)*Weight, color=Colors[0], alpha=0.5, label='Extrap.')
        ax[0, 1].hist(Porb, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), weights=np.ones_like(Porb)*Weight, color=Colors[1], alpha=0.5, label='Original')
        ax[0, 1].set_xscale('log')
        ax[0, 1].grid(axis='y')
        ax[0, 1].set_xlabel('Planet orbital period [d]')
        ax[0, 1].set_ylabel('Fraction')
        ax[0, 1].legend()
        plt.subplot(212, projection='3d')
        plt.gca().plot_surface(self.xx, self.yy, self.zz)
        plt.gca().plot_wireframe(self.xx, self.yy, self.r_bqf, color='black')
        plt.gca().set_xlabel('Log(planet orbital period)')
        plt.gca().set_ylabel('Log(planet radius)')
        plt.gca().set_zlabel('Log(rates)')
        plt.gca().view_init(elev=30., azim=115.)
        plt.suptitle('Fressin2013')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'PlanetDistribution_Fressin2013.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
