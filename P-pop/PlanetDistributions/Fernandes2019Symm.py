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
# FERNANDES2019SYMM
# =============================================================================

class PlanetDistribution():
    """
    https://ui.adsabs.harvard.edu/abs/2019ApJ...874...81F/abstract
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
        print('--> Initializing Fernandes2019Symm planet distribution')
        
        # Model parameters.
        self.returns = ['Mp', 'Porb']
        if (Scenario == 'baseline'):
            self.m1 = -0.45
            self.p1 = 0.65
            self.p2 = -0.65
            self.Pbrk = 1581. # d
            self.c0 = 0.84 / 10. # normalization
        elif (Scenario == 'pessimistic'):
            self.m1 = -0.45
            self.p1 = 0.65
            self.p2 = -0.65
            self.Pbrk = 1581. # d
            self.c0 = (0.84-0.15) / 10. # normalization
        elif (Scenario == 'optimistic'):
            self.m1 = -0.45
            self.p1 = 0.65
            self.p2 = -0.65
            self.Pbrk = 1581. # d
            self.c0 = (0.84+0.18) / 10. # normalization
        elif (Scenario == 'mc'):
            raise UserWarning('Scenario mc is not implemented yet')
        else:
            print('--> WARNING: '+str(Scenario)+' is an unknown scenario')
            Scenario = 'baseline'
            self.m1 = -0.45
            self.p1 = 0.65
            self.p2 = -0.65
            self.Pbrk = 1581. # d
            self.c0 = 0.84
        print('--> Using scenario '+str(Scenario))
        self.Mp_lims = [3., 600.] # 10 Mearth
        self.Porb_lims = [10., 10000.] # d
        self.CM = 1./(self.Mp_lims[1]**self.m1/self.m1-self.Mp_lims[0]**self.m1/self.m1)
        self.CP1 = self.Pbrk**self.p1/(self.Pbrk**self.p1/self.p1-self.Porb_lims[0]**self.p1/self.p1)
        self.CP2 = self.Pbrk**self.p2/(self.Porb_lims[1]**self.p2/self.p2-self.Pbrk**self.p2/self.p2)
        self.F1 = self.c0/self.Pbrk**self.p1*(self.Mp_lims[1]**self.m1/self.m1-self.Mp_lims[0]**self.m1/self.m1)*(self.Pbrk**self.p1/self.p1-self.Porb_lims[0]**self.p1/self.p1)
        self.F2 = self.c0/self.Pbrk**self.p2*(self.Mp_lims[1]**self.m1/self.m1-self.Mp_lims[0]**self.m1/self.m1)*(self.Porb_lims[1]**self.p2/self.p2-self.Pbrk**self.p2/self.p2)
        
        # from scipy.integrate import dblquad
        
        # # Msun = 1.989e30
        # # G = 6.674e-11
        
        # # M = Msun
        
        # # Mmin = 0.1 * 317.8
        # # Mmax = 13. * 317.8
        
        # # amin = 10. * 1.496e11
        # # Pmin = np.sqrt(4.*np.pi**2/(G*M)*amin**3) / 3600. / 24.
        # # amax = 100. * 1.496e11
        # # Pmax = np.sqrt(4.*np.pi**2/(G*M)*amax**3) / 3600. / 24.
        
        # def func_0(M, P):
        #     return self.c0 * (P / self.Pbrk)**self.p1 * (M / 10.)**self.m1 / P / M
        # F0_0 = dblquad(func_0, 10., self.Pbrk, 30., 6000.)
        
        # def func_1(M, P):
        #     return self.c0 * (P / self.Pbrk)**self.p2 * (M / 10.)**self.m1 / P / M
        # F0_1 = dblquad(func_1, self.Pbrk, 10000., 30., 6000.)
        
        # Ndraw = 100000
        # Mps = []
        # Porbs = []
        # for i in range(Ndraw):
        #     Mp, Porb = self.draw(Mp_range=[30., 6000.], Porb_range=[10., 10000.])
        #     Mps += [Mp]
        #     Porbs += [Porb]
        # Mps = np.concatenate(Mps)
        # Porbs = np.concatenate(Porbs)
        # eta = len(Mps) / Ndraw
        
        # # def dNdlogP(P):
        # #     if P < self.Pbrk:
        # #         return self.c0 * (P / self.Pbrk)**self.p1 * (1. / 10.)**self.m1 * (6000.**self.m1 / self.m1 - 30.**self.m1 / self.m1)
        # #     else:
        # #         return self.c0 * (P / self.Pbrk)**self.p2 * (1. / 10.)**self.m1 * (6000.**self.m1 / self.m1 - 30.**self.m1 / self.m1)
        
        # # plt.close()
        # # xx = np.linspace(1., 4., 1000)
        # # yy = np.array([dNdlogP(10**x) for x in xx])
        # # plt.plot(xx, yy)
        # # plt.yscale('log')
        # # plt.show()
        
        # import pdb; pdb.set_trace()
        
        pass
    
    def iCDF_M(self,
               x):
        """
        Parameters
        ----------
        x: float, 0 <= x <= 1
            Uniformly distributed random number.
        
        Returns
        -------
        Mp: float
            Planet mass (Mearth) distributed according to Fernandes2019Symm
            planet distribution.
        """
        
        return (self.m1*x/self.CM+self.Mp_lims[0]**self.m1)**(1./self.m1)*10.
    
    def iCDF_P1(self,
                x):
        """
        Parameters
        ----------
        x: float, 0 <= x <= 1
            Uniformly distributed random number.
        
        Returns
        -------
        Porb: float
            Planet orbital period (d) given Porb < Pbrk distributed according
            to Fernandes2019Symm planet distribution.
        """
        
        return (self.Pbrk**self.p1*self.p1*x/self.CP1+self.Porb_lims[0]**self.p1)**(1./self.p1)
    
    def iCDF_P2(self,
                x):
        """
        Parameters
        ----------
        x: float, 0 <= x <= 1
            Uniformly distributed random number.
        
        Returns
        -------
        Porb: float
            Planet orbital period (d) given Porb > Pbrk distributed according
            to Fernandes2019Symm planet distribution.
        """
        
        return (self.Pbrk**self.p2*self.p2*x/self.CP2+self.Pbrk**self.p2)**(1./self.p2)
    
    def draw(self,
             Mp_range=[30., 6000.], # Mearth
             Porb_range=[10., 10000.], # d
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
        tempF0 = np.array([self.F1*Scale, self.F2*Scale])
        
        # If the number of planets is not given, draw it from a Poisson
        # distribution. Note that the final number of drawn planets might be
        # smaller than the drawn number because of clipping to the requested Rp
        # and Porb range.
        if (Nplanets is None):
            Nplanets = np.random.poisson(np.sum(tempF0))
            for i in range(Nplanets):
                
                # Randomly select whether Porb < Pbrk or Porb > Pbrk.
                temp = np.random.choice(len(tempF0), p=tempF0/np.sum(tempF0))
                if (temp == 0):
                    tempMp = self.iCDF_M(np.random.rand()) # Mearth
                    tempPorb = self.iCDF_P1(np.random.rand()) # d
                elif (temp == 1):
                    tempMp = self.iCDF_M(np.random.rand()) # Mearth
                    tempPorb = self.iCDF_P2(np.random.rand()) # d
                if (Mp_range[0] <= tempMp <= Mp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Mp += [tempMp] # Mearth
                    Porb += [tempPorb] # d
        
        # If the number of planets is given, draw exactly this number of
        # planets in the requested Rp and Porb range.
        else:
            while (len(Mp) < Nplanets):
                
                # Randomly select whether Rp < Rbrk or Rp > Rbrk.
                temp = np.random.choice(len(tempF0), p=tempF0/np.sum(tempF0))
                if (temp == 0):
                    tempMp = self.iCDF_M(np.random.rand()) # Mearth
                    tempPorb = self.iCDF_P1(np.random.rand()) # d
                elif (temp == 1):
                    tempMp = self.iCDF_M(np.random.rand()) # Mearth
                    tempPorb = self.iCDF_P2(np.random.rand()) # d
                if (Mp_range[0] <= tempMp <= Mp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Mp += [tempMp] # Mearth
                    Porb += [tempPorb] # d
        
        return np.array(Mp), np.array(Porb)
    
    def SummaryPlot(self,
                    Ntest=100000,
                    Mp_range=[30., 6000.], # Mearth
                    Porb_range=[10., 10000.], # d
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
        
        print('--> Fernandes2019Symm:\n%.2f planets/star in [%.1f, %.1f] Mearth and [%.1f, %.1f] d' % (len(Mp)/float(Ntest), Mp_range[0], Mp_range[1], Porb_range[0], Porb_range[1]))
        
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
        plt.suptitle('Fernandes2019Symm')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'PlanetDistribution_Fernandes2019Symm.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
