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
import os

import scipy.stats as stats

import sys 
sys.path.append('..')
import Star


# =============================================================================
# BRYSON2021MODEL1HAB2HIGH
# =============================================================================

class PlanetDistribution():
    """
    https://ui.adsabs.harvard.edu/abs/2021AJ....161...36B/abstract
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
        print('--> Initializing Bryson2021Model1Hab2High planet distribution')
        
        # Constants.
        self.G = 6.674e-11 # m^3/kg/s^2
        self.Msun = 1.989e30 # kg
        self.Tsun = 5778. # K
        self.au = 149597870700. # m
        
        # Model parameters.
        self.returns = ['Rp', 'Porb']
        self.Rp_lims = [0.5, 2.5] # Rearth
        self.Fp_lims = [0.2, 2.2] # Searth
        # Model 1 hab stars low.
        # self.Teff_lims = [4800., 6300.] # K
        # self.F0s = np.array([1.08, 1.56, 0.57])
        # self.alphas = np.array([-1.05, 1.41, 1.20])
        # self.betas = np.array([-0.56, 0.48, 0.42])
        # self.gammas = np.array([-1.84, 3.33, 3.39])
        # Model 1 hab stars high.
        # self.Teff_lims = [4800., 6300.] # K
        # self.F0s = np.array([1.97, 3.73, 1.17])
        # self.alphas = np.array([-1.09, 1.36, 1.18])
        # self.betas = np.array([-1.18, 0.60, 0.56])
        # self.gammas = np.array([0.91, 3.87, 3.88])
        # Model 1 hab2 stars low.
        # self.Teff_lims = [3900., 6300.] # K
        # self.F0s = np.array([1.11, 0.88, 0.44])
        # self.alphas = np.array([-1.08, 0.94, 0.85])
        # self.betas = np.array([-0.84, 0.32, 0.30])
        # self.gammas = np.array([-2.67, 1.59, 1.57])
        # Model 1 hab2 stars high.
        self.Teff_lims = [3900., 6300.] # K
        self.F0s = np.array([1.59, 1.56, 0.70])
        self.alphas = np.array([-1.18, 0.96, 0.87])
        self.betas = np.array([-1.19, 0.37, 0.36])
        self.gammas = np.array([-1.38, 1.84, 1.78])
        if (Scenario == 'baseline'):
            self.F0 = self.F0s[0].copy()
            self.alpha = self.alphas[0].copy()
            self.beta = self.betas[0].copy()
            self.gamma = self.gammas[0].copy()
        elif (Scenario == 'optimistic'):
            self.F0 = self.F0s[0]+self.F0s[1]
            self.alpha = self.alphas[0].copy()
            self.beta = self.betas[0].copy()
            self.gamma = self.gammas[0].copy()
        elif (Scenario == 'pessimistic'):
            self.F0 = self.F0s[0]-self.F0s[2]
            self.alpha = self.alphas[0].copy()
            self.beta = self.betas[0].copy()
            self.gamma = self.gammas[0].copy()
        elif (Scenario == 'mc'):
            self.F0 = self.F0s[0].copy()
            self.alpha = self.alphas[0].copy()
            self.beta = self.betas[0].copy()
            self.gamma = self.gammas[0].copy()
            # self.mcpost = np.load(os.path.join(os.path.split(os.path.abspath(__file__))[0], 'BrysonPosteriorModel1Hab2Low.npy'))
            self.mcpost = np.load(os.path.join(os.path.split(os.path.abspath(__file__))[0], 'BrysonPosteriorModel1Hab2High.npy'))
        else:
            print('--> WARNING: '+str(Scenario)+' is an unknown scenario')
            Scenario = 'baseline'
            self.F0 = self.F0s[0].copy()
            self.alpha = self.alphas[0].copy()
            self.beta = self.betas[0].copy()
            self.gamma = self.gammas[0].copy()
            print('--> Using scenario '+str(Scenario))
        self.CR = 1./(self.Rp_lims[1]**(self.alpha+1.)/(self.alpha+1.)-self.Rp_lims[0]**(self.alpha+1.)/(self.alpha+1.))
        self.CF = 1./(self.Fp_lims[1]**(self.beta+1.)/(self.beta+1.)-self.Fp_lims[0]**(self.beta+1.)/(self.beta+1.))
        self.CT = 1./(((self.Teff_lims[1]**(self.gamma+5.49)/(self.gamma+5.49)-5117.**(self.gamma+5.49)/(self.gamma+5.49))*10.**(-16.77)/self.Tsun**self.gamma+(5117.**(self.gamma+4.16)/(self.gamma+4.16)-self.Teff_lims[0]**(self.gamma+4.16)/(self.gamma+4.16))*10.**(-11.84)/self.Tsun**self.gamma)/(self.Teff_lims[1]-self.Teff_lims[0]))
        self.Scenario = Scenario
        
        # import matplotlib
        # matplotlib.rcParams.update({'font.size': 14})
        # def func(Teff):
        #     return (Teff/self.Tsun)**self.gamma*self.g(Teff)*self.CT
        # xx = np.linspace(3900., 6300., 10000)
        # yy = [func(x) for x in xx]
        # plt.close()
        # plt.plot(xx, yy)
        # plt.xlabel('Host star effective temperature [K]', size=14)
        # plt.ylabel('Fractional occurrence rate', size=14)
        # plt.show()
        
        # import pdb; pdb.set_trace()
        
        pass
    
    def redraw_params(self):
        """
        """
        
        ww = np.random.randint(0, self.mcpost.shape[0])
        
        self.F0 = self.mcpost[ww, 0]
        self.alpha = self.mcpost[ww, 2]
        self.beta = self.mcpost[ww, 1]
        self.gamma = self.mcpost[ww, 3]
        
        self.CR = 1./(self.Rp_lims[1]**(self.alpha+1.)/(self.alpha+1.)-self.Rp_lims[0]**(self.alpha+1.)/(self.alpha+1.))
        self.CF = 1./(self.Fp_lims[1]**(self.beta+1.)/(self.beta+1.)-self.Fp_lims[0]**(self.beta+1.)/(self.beta+1.))
        self.CT = 1./(((self.Teff_lims[1]**(self.gamma+5.49)/(self.gamma+5.49)-5117.**(self.gamma+5.49)/(self.gamma+5.49))*10.**(-16.77)/self.Tsun**self.gamma+(5117.**(self.gamma+4.16)/(self.gamma+4.16)-self.Teff_lims[0]**(self.gamma+4.16)/(self.gamma+4.16))*10.**(-11.84)/self.Tsun**self.gamma)/(self.Teff_lims[1]-self.Teff_lims[0]))
        
        pass
    
    def iCDF_R(self,
               x):
        """
        Parameters
        ----------
        x: float, 0 <= x <= 1
            Uniformly distributed random number.
        
        Returns
        -------
        Rp: float
            Planet radius (Rearth) distributed according to Bryson et al. 2021
            planet distribution.
        """
        
        return ((self.alpha+1.)*x/self.CR+self.Rp_lims[0]**(self.alpha+1.))**(1./(self.alpha+1.))
    
    def iCDF_F(self,
               x):
        """
        Parameters
        ----------
        x: float, 0 <= x <= 1
            Uniformly distributed random number.
        
        Returns
        -------
        Fp: float
            Planet incident host star flux (Searth) distributed according to
            Bryson et al. 2021 planet distribution.
        """
        
        return ((self.beta+1.)*x/self.CF+self.Fp_lims[0]**(self.beta+1.))**(1./(self.beta+1.))
    
    def g(self,
          Teff): # K
        """
        Parameters
        ----------
        Teff: float
            Host star effective temperature (K).
        
        Returns
        -------
        g(Teff): float
            Equation 4 of Bryson et al. 2021.
        """
        
        if (Teff < 5117.):
            return 10**(-11.84)*Teff**3.16
        else:
            return 10**(-16.77)*Teff**4.49
    
    def draw(self,
             Rp_range=[0.5, 2.5], # Rearth
             Fp_range=[0.2, 2.2], # Searth
             Nplanets=None,
             Scale=1.,
             Star=None):
        """
        Parameters
        ----------
        Rp_range: list
            Requested planet radius range (Rearth).
        Fp_range: list
            Requested planet incident host star flux (Searth).
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
        tempF0 = np.array(self.F0).copy()*(Star.Teff/self.Tsun)**self.gamma*self.g(Star.Teff)*self.CT*Scale
        
        # If the number of planets is not given, draw it from a Poisson
        # distribution. Note that the final number of drawn planets might be
        # smaller than the drawn number because of clipping to the requested Rp
        # and Porb range.
        if (Nplanets is None):
            Nplanets = np.random.poisson(tempF0)
            for i in range(Nplanets):
                tempRp = self.iCDF_R(np.random.rand()) # Rearth
                tempFp = self.iCDF_F(np.random.rand()) # Searth
                if (Rp_range[0] <= tempRp <= Rp_range[1] and Fp_range[0] <= tempFp <= Fp_range[1]):
                    tempap = np.sqrt(Star.Rad**2*(Star.Teff/self.Tsun)**4/tempFp) # au
                    tempPorb = np.sqrt(4.*np.pi**2*(tempap*self.au)**3/self.G/(Star.Mass*self.Msun))/86400. # d
                    Rp += [tempRp] # Rearth
                    Porb += [tempPorb] # d
        
        # If the number of planets is given, draw exactly this number of
        # planets in the requested Rp and Porb range.
        else:
            while (len(Rp) < Nplanets):
                tempRp = self.iCDF_R(np.random.rand()) # Rearth
                tempFp = self.iCDF_F(np.random.rand()) # Searth
                if (Rp_range[0] <= tempRp <= Rp_range[1] and Fp_range[0] <= tempFp <= Fp_range[1]):
                    tempap = np.sqrt(Star.Rad**2*(Star.Teff/self.Tsun)**4/tempFp) # au
                    tempPorb = np.sqrt(4.*np.pi**2*(tempap*self.au)**3/self.G/(Star.Mass*self.Msun))/86400. # d
                    Rp += [tempRp] # Rearth
                    Porb += [tempPorb] # d
        
        return np.array(Rp), np.array(Porb)
    
    def SummaryPlot(self,
                    Ntest=100000,
                    Rp_range=[0.5, 2.5], # Rearth
                    Fp_range=[0.2, 2.2], # Searth
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
        
        Sun = Star.Star('Sun',
                        10., # pc
                        'G',
                        1., # Rsun
                        5778., # K
                        1., # Msun
                        RA=None, # deg
                        Dec=None, # deg
                        Vmag=None, # mag
                        Jmag=None, # mag
                        Hmag=None, # mag
                        WDSsep=None, # au
                        WDSdmag=None, # mag
                        lGal=None, # deg
                        bGal=None) # deg
        
        Rp = []
        Porb = []
        for i in range(Ntest):
            tempRp, tempPorb = self.draw(Rp_range,
                                         Fp_range,
                                         Star=Sun)
            Rp += [tempRp]
            Porb += [tempPorb]
        Rp = np.concatenate(Rp)
        Porb = np.concatenate(Porb)
        
        print('--> Bryson2021Model1Hab2High:\n%.2f planets/star in [%.1f, %.1f] Rearth and [%.1f, %.1f] Searth' % (len(Rp)/float(Ntest), Rp_range[0], Rp_range[1], Fp_range[0], Fp_range[1]))
        
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
        plt.suptitle('Bryson2021Model1Hab2High')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'PlanetDistribution_Bryson2021Model1Hab2High.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
