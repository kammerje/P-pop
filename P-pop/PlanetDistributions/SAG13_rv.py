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
from scipy.integrate import quad


# =============================================================================
# SAG13
# =============================================================================

class rv_R(stats.rv_continuous):
    
    def __init__(self, **kwargs):
        
        
        super(rv_R, self).__init__(**kwargs)
        
        pass
    
    def dNdR(self, Rp, Rbrk, alpha):
        
        return Rp**(alpha - 1.)
    
    def _pdf(self,  Rp, Rbrk, alpha, norm):
        
        return self.dNdR(Rp, Rbrk, alpha) / norm
    
    def _argcheck(self, *args):
        
        return 1

class rv_P(stats.rv_continuous):
    
    def __init__(self, **kwargs):
        
        
        super(rv_P, self).__init__(**kwargs)
        
        pass
    
    def dNdP(self, Porb, beta):
        
        return Porb**(beta - 1.)
    
    def _pdf(self,  Porb, beta, norm):
        
        return self.dNdP(Porb, beta) / norm
    
    def _argcheck(self, *args):
        
        return 1

class PlanetDistribution():
    """
    https://ui.adsabs.harvard.edu/abs/2018ApJ...856..122K/abstract
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
        print('--> Initializing SAG13 planet distribution')
        
        # Model parameters.
        self.returns = ['Rp', 'Porb']
        if (Scenario == 'baseline'):
            self.F0 = [2.42, 0.25]
            self.Gamma = [0.38, 0.73]
            self.alpha = [-0.19, -1.18]
            self.beta = [0.26, 0.59]
        elif (Scenario == 'pessimistic'):
            self.F0 = [1.14, 0.14]
            self.Gamma = [0.138, 0.72]
            self.alpha = [0.277, -1.56]
            self.beta = [0.204, 0.51]
        elif (Scenario == 'optimistic'):
            self.F0 = [5.60, 0.46]
            self.Gamma = [1.06, 0.78]
            self.alpha = [-0.68, -0.82]
            self.beta = [0.32, 0.67]
        elif (Scenario == 'mc'):
            raise UserWarning('Scenario mc is not implemented yet')
        else:
            print('--> WARNING: '+str(Scenario)+' is an unknown scenario')
            Scenario = 'baseline'
            self.F0 = [2.42, 0.25]
            self.Gamma = [0.38, 0.73]
            self.alpha = [-0.19, -1.18]
            self.beta = [0.26, 0.59]
        print('--> Using scenario '+str(Scenario))
        self.Rbrk = [0., 3.4, np.inf] # Rearth
        self.ytod = 365.24
        self.Rp_lims = [0.5, 16.] # Rearth
        self.Porb_lims = [0.5, 500.] # d
        self.CR0 = 1./(self.Rbrk[1]**self.alpha[0]/self.alpha[0]-self.Rp_lims[0]**self.alpha[0]/self.alpha[0])
        self.CP0 = 1./((self.Porb_lims[1]/self.ytod)**self.beta[0]/self.beta[0]-(self.Porb_lims[0]/self.ytod)**self.beta[0]/self.beta[0])
        self.CR1 = 1./(self.Rp_lims[1]**self.alpha[1]/self.alpha[1]-self.Rbrk[1]**self.alpha[1]/self.alpha[1])
        self.CP1 = 1./((self.Porb_lims[1]/self.ytod)**self.beta[1]/self.beta[1]-(self.Porb_lims[0]/self.ytod)**self.beta[1]/self.beta[1])
        
        self.CR_less = quad(self.dNdR, self.Rp_lims[0], self.Rbrk[1], args=(self.Rbrk[1], self.alpha[0]))
        self.CR_more = quad(self.dNdR, self.Rbrk[1], self.Rp_lims[1], args=(self.Rbrk[1], self.alpha[1]))
        self.CP_less = quad(self.dNdP, self.Porb_lims[0] / self.ytod, self.Porb_lims[1] / self.ytod, args=(self.beta[0]))
        self.CP_more = quad(self.dNdP, self.Porb_lims[0] / self.ytod, self.Porb_lims[1] / self.ytod, args=(self.beta[1]))
        
        pass
    
    def dNdR(self, Rp, Rbrk, alpha):
        
        return Rp**(alpha - 1.)
    
    def dNdP(self, Porb, beta):
        
        return Porb**(beta - 1.)
    
    def iCDF_R0(self,
                x):
        """
        Parameters
        ----------
        x: float, 0 <= x <= 1
            Uniformly distributed random number.
        
        Returns
        -------
        Rp: float
            Planet radius (Rearth) given Rp < Rbrk distributed according to
            SAG13 planet distribution.
        """
        
        return (self.alpha[0]*x/self.CR0+self.Rp_lims[0]**self.alpha[0])**(1./self.alpha[0])
    
    def iCDF_P0(self,
                x):
        """
        Parameters
        ----------
        x: float, 0 <= x <= 1
            Uniformly distributed random number.
        
        Returns
        -------
        Porb: float
            Planet orbital period (d) given Rp < Rbrk distributed according to
            SAG13 planet distribution.
        """
        
        return (self.beta[0]*x/self.CP0+(self.Porb_lims[0]/self.ytod)**self.beta[0])**(1./self.beta[0])
    
    def iCDF_R1(self,
                x):
        """
        Parameters
        ----------
        x: float, 0 <= x <= 1
            Uniformly distributed random number.
        
        Returns
        -------
        Rp: float
            Planet radius (Rearth) given Rp > Rbrk distributed according to
            SAG13 planet distribution.
        """
        
        return (self.alpha[1]*x/self.CR1+self.Rbrk[1]**self.alpha[1])**(1./self.alpha[1])
    
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
            Planet orbital period (d) given Rp > Rbrk distributed according to
            SAG13 planet distribution.
        """
        
        return (self.beta[1]*x/self.CP1+(self.Porb_lims[0]/self.ytod)**self.beta[1])**(1./self.beta[1])
    
    def draw_rv(self,
                Rp_range=[0.5, 16.], # Rearth
                Porb_range=[0.5, 500.], # d
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
        tempF0 = np.array(self.F0).copy()*Scale
        
        # If the number of planets is not given, draw it from a Poisson
        # distribution. Note that the final number of drawn planets might be
        # smaller than the drawn number because of clipping to the requested Rp
        # and Porb range.
        if (Nplanets is None):
            Nplanets = np.random.poisson(np.sum(tempF0))
            for i in range(Nplanets):
                
                # Randomly select whether Rp < Rbrk or Rp > Rbrk.
                temp = np.random.choice(len(tempF0), p=tempF0/np.sum(tempF0))
                if (temp == 0):
                    rv = rv_R(a=self.Rp_lims[0], b=self.Rbrk[1])
                    tempRp = rv.rvs(Rbrk=self.Rbrk[1], alpha=self.alpha[0], norm=self.CR_less)[0] # Rearth
                    rv = rv_P(a=self.Porb_lims[0] / self.ytod, b=self.Porb_lims[1] / self.ytod)
                    tempPorb = rv.rvs(beta=self.beta[0], norm=self.CP_less)[0] * self.ytod # d
                elif (temp == 1):
                    rv = rv_R(a=self.Rbrk[1], b=self.Rp_lims[1])
                    tempRp = rv.rvs(Rbrk=self.Rbrk[1], alpha=self.alpha[1], norm=self.CR_more)[0] # Rearth
                    rv = rv_P(a=self.Porb_lims[0] / self.ytod, b=self.Porb_lims[1] / self.ytod)
                    tempPorb = rv.rvs(beta=self.beta[1], norm=self.CP_more)[0] * self.ytod  # d
                if (Rp_range[0] <= tempRp <= Rp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Rp += [tempRp] # Rearth
                    Porb += [tempPorb] # d
        
        # If the number of planets is given, draw exactly this number of
        # planets in the requested Rp and Porb range.
        else:
            while (len(Rp) < Nplanets):
                
                # Randomly select whether Rp < Rbrk or Rp > Rbrk.
                temp = np.random.choice(len(tempF0), p=tempF0/np.sum(tempF0))
                if (temp == 0):
                    tempRp = self.iCDF_R0(np.random.rand()) # Rearth
                    tempPorb = self.iCDF_P0(np.random.rand())*self.ytod # d
                elif (temp == 1):
                    tempRp = self.iCDF_R1(np.random.rand()) # Rearth
                    tempPorb = self.iCDF_P1(np.random.rand())*self.ytod # d
                if (Rp_range[0] <= tempRp <= Rp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Rp += [tempRp] # Rearth
                    Porb += [tempPorb] # d
        
        return np.array(Rp), np.array(Porb)
    
    def draw(self,
             Rp_range=[0.5, 16.], # Rearth
             Porb_range=[0.5, 500.], # d
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
        tempF0 = np.array(self.F0).copy()*Scale
        
        # If the number of planets is not given, draw it from a Poisson
        # distribution. Note that the final number of drawn planets might be
        # smaller than the drawn number because of clipping to the requested Rp
        # and Porb range.
        if (Nplanets is None):
            Nplanets = np.random.poisson(np.sum(tempF0))
            for i in range(Nplanets):
                
                # Randomly select whether Rp < Rbrk or Rp > Rbrk.
                temp = np.random.choice(len(tempF0), p=tempF0/np.sum(tempF0))
                if (temp == 0):
                    tempRp = self.iCDF_R0(np.random.rand()) # Rearth
                    tempPorb = self.iCDF_P0(np.random.rand())*self.ytod # d
                elif (temp == 1):
                    tempRp = self.iCDF_R1(np.random.rand()) # Rearth
                    tempPorb = self.iCDF_P1(np.random.rand())*self.ytod # d
                if (Rp_range[0] <= tempRp <= Rp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Rp += [tempRp] # Rearth
                    Porb += [tempPorb] # d
        
        # If the number of planets is given, draw exactly this number of
        # planets in the requested Rp and Porb range.
        else:
            while (len(Rp) < Nplanets):
                
                # Randomly select whether Rp < Rbrk or Rp > Rbrk.
                temp = np.random.choice(len(tempF0), p=tempF0/np.sum(tempF0))
                if (temp == 0):
                    tempRp = self.iCDF_R0(np.random.rand()) # Rearth
                    tempPorb = self.iCDF_P0(np.random.rand())*self.ytod # d
                elif (temp == 1):
                    tempRp = self.iCDF_R1(np.random.rand()) # Rearth
                    tempPorb = self.iCDF_P1(np.random.rand())*self.ytod # d
                if (Rp_range[0] <= tempRp <= Rp_range[1] and Porb_range[0] <= tempPorb <= Porb_range[1]):
                    Rp += [tempRp] # Rearth
                    Porb += [tempPorb] # d
        
        return np.array(Rp), np.array(Porb)
    
    def SummaryPlot(self,
                    Ntest=100000,
                    Rp_range=[0.5, 16.], # Rearth
                    Porb_range=[0.5, 500.], # d
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
        
        Rp_rv = []
        Porb_rv = []
        for i in range(Ntest):
            tempRp, tempPorb = self.draw_rv(Rp_range,
                                            Porb_range)
            Rp_rv += [tempRp]
            Porb_rv += [tempPorb]
        Rp_rv = np.concatenate(Rp_rv)
        Porb_rv = np.concatenate(Porb_rv)
        
        print('--> SAG13:\n%.2f planets/star in [%.1f, %.1f] Rearth and [%.1f, %.1f] d' % (len(Rp)/float(Ntest), Rp_range[0], Rp_range[1], Porb_range[0], Porb_range[1]))
        
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
        plt.suptitle('SAG13')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'PlanetDistribution_SAG13.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
