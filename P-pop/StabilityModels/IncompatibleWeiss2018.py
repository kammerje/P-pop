from __future__ import division


# =============================================================================
# IMPORTS
# =============================================================================

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


# =============================================================================
# WEISS2018
# =============================================================================

class StabilityModel():
    """
    https://ui.adsabs.harvard.edu/abs/2018AJ....155...48W/abstract
    """
    
    def __init__(self):
        
        # Max number of re-draws if correlation is too extreme or system is
        # unstable
        self.MaxTrials = np.inf
        
        # Mutual Hill radius above which planets are considered to be stable
        self.Threshold = 8.
        
        # Constants
        self.yr2d = 365.24
        self.Mearth = 5.972e24 # kg
        self.Msun = 1.989e30 # kg
        
        pass
    
    def draw(self,
             Ms, # Msun
             PlanetDistribution,
             MassModel,
             EccentricityModel,
             returnTrials=False):
        """
        Parameters
        ----------
        Ms: float
            Mass (Msun) of the star.
        PlanetDistribution: instance
            Instance of class PlanetDistribution.
        MassModel: instance
            Instance of class MassModel.
        EccentricityModel: instance
            Instance of class EccentricityModel.
        returnTrials: bool
            If True, returns number of trials required for drawing a stable
            system.
        
        Returns
        -------
        Rp: array
            Radius (Rearth) of drawn planets of stable system.
        Porb: array
            Orbital period (d) of drawn planets of stable system.
        Mp: array
            Mass (Mearth) of drawn planets of stable system.
        ep: array
            Eccentricity of drawn planets of stable system.
        _Trials: int
            Number of trials required for drawing a stable system.
        """
        
        # Draw planet radius and planet orbital period
        Rp, Porb = PlanetDistribution.draw() # Rearth, d
        
        # Sort planets by orbital period
        ww = np.argsort(Porb)
        Rp = Rp[ww] # Rearth
        Porb = Porb[ww] # d
        
        # Number of planets
        Nplanets = len(Rp)
        Trials = 1
        
        # Multi-planet system --> check stability
        if (Nplanets > 1):
            
            # Draw planet radius and planet orbital period according to the
            # correlations from Weiss et al. 2018 for multi-planet systems
            Rp = self.RadiusCorrSkew(Rp) # Rearth
            Porb = self.PeriodCorrSkew(Porb) # d
            
            # Sort planets by orbital period
            ww = np.argsort(Porb)
            Rp = Rp[ww] # Rearth
            Porb = Porb[ww] # d
            
            # Draw planet mass and planet eccentricity
            Mp = MassModel.RadiusToMass(Rp) # Mearth
            ep = EccentricityModel.getEccentricity(Porb)
            
            # Check system stability
            Stability = self.CheckSystemStability(Ms,
                                                  Porb,
                                                  Mp,
                                                  ep)
            
            # If the system is unstable, draw a new system
            while (Trials <= self.MaxTrials and Stability == False):
                
                # Draw planet radius and planet orbital period; keep the same
                # number of planets
                Rp, Porb = PlanetDistribution.draw(Nplanets=Nplanets) # Rearth, d
                
                # Sort planets by orbital period
                ww = np.argsort(Porb)
                Rp = Rp[ww] # Rearth
                Porb = Porb[ww] # d
                
                # Draw planet radius and planet orbital period according to
                # the correlations from Weiss et al. 2018 for multi-planet
                # systems
                Rp = self.RadiusCorrSkew(Rp) # Rearth
                Porb = self.PeriodCorrSkew(Porb) # d
                
                # Sort planets by orbital period
                ww = np.argsort(Porb)
                Rp = Rp[ww] # Rearth
                Porb = Porb[ww] # d
                
                # Draw planet mass and planet eccentricity
                Mp = MassModel.RadiusToMass(Rp) # Mearth
                ep = EccentricityModel.getEccentricity(Porb)
                
                # Check system stability
                Stability = self.CheckSystemStability(Ms,
                                                      Porb,
                                                      Mp,
                                                      ep)
                Trials += 1
        
        # Single-planet system
        else:
            
            # Draw planet mass and planet eccentricity
            Mp = MassModel.RadiusToMass(Rp) # Mearth
            ep = EccentricityModel.getEccentricity(Porb)
        
        if (returnTrials == False):
            return Rp, Porb, Mp, ep
        else:
            return Rp, Porb, Mp, ep, Trials
    
    def RadiusCorrSkew(self,
                       Rp, # Rearth
                       Ratio_lims=[0., np.inf],
                       Rp_lims=[0.1, 100.]): # Rearth
        """
        Parameters
        ----------
        Rp: array
            Radius (Rearth) of drawn planets sorted by orbital period.
        Ratio_lims: list
            Requested planet radius ratio range.
        Rp_lims: list
            Requested planet radius range (Rearth).
        
        Returns
        -------
        Rp: array
            Radius (Rearth) of drawn planets according to Weiss et al. 2018.
        """
        
        # Parameters from Emile
        mu, sigma, skew = 0.69007816, 0.70029624, 3.52948576
        
        for i in range(len(Rp)-1):
            
            # Draw radius ratio inside limits
            RadiusRatio = -1.
            Rp[i+1] = -1.
            Trials = 1
            while (Trials <= self.MaxTrials and (RadiusRatio <= Ratio_lims[0] or RadiusRatio >= Ratio_lims[1] or Rp[i+1] <= Rp_lims[0] or Rp[i+1] >= Rp_lims[1])):
                RadiusRatio = stats.skewnorm(skew, loc=mu, scale=sigma).rvs()
                Rp[i+1] = Rp[i]*RadiusRatio # Rearth
                Trials += 1
        
        return Rp
    
    def PeriodCorrSkew(self,
                       Porb, # d
                       Ratio_lims=[0., np.inf]):
        """
        Parameters
        ----------
        Porb: array
            Orbital period (d) of drawn planets sorted by orbital period.
        Ratio_lims: list
            Requested planet orbital period ratio range.
        
        Returns
        -------
        Porb: array
            Orbital period (d) of drawn planets according to Weiss et al. 2018.
        """
        
        # Parameters from Emile
        mu, sigma, skew = 1.28170899, 0.91055566, 9.35402502
        
        for i in range(len(Porb)-1):
            
            # Draw period ratio inside limits
            PeriodRatio = -1.
            Trials = 1
            while (Trials <= self.MaxTrials and (PeriodRatio <= Ratio_lims[0] or PeriodRatio >= Ratio_lims[1])):
                PeriodRatio = stats.skewnorm(skew, loc=mu, scale=sigma).rvs()
                Porb[i+1] = Porb[i]*PeriodRatio # d
                Trials += 1
        
        return Porb
    
    def CheckSystemStability(self,
                             Ms, # Msun
                             Porb, # d
                             Mp, # Mearth
                             ep):
        """
        Parameters
        ----------
        Ms: float
            Mass of the star (Msun).
        Porb: array
            Orbital period (d) of drawn planets sorted by orbital period.
        Mp: array
            Mass (Mearth) of drawn planets sorted by orbital period.
        ep: array
            Eccentricity of drawn planets sorted by orbital period.
        
        Returns
        -------
        Stability: bool
            If True, the system is stable.
        """
        
        # Compute semi-major axes in au
        tempPorb = Porb/self.yr2d # yr
        ap = tempPorb**(2./3.)*Ms**(1./3.) # au
        
        # Check stability between each pair of planets
        for i in range(len(tempPorb)-1):
            Stability = self.CheckStability(Ms,
                                            ap[i],
                                            ap[i+1],
                                            Mp[i],
                                            Mp[i+1],
                                            ep[i],
                                            ep[i+1])
            if (Stability == False):
                return False
        return True
    
    def CheckStability(self,
                       Ms, # Msun
                       ap_in, # au
                       ap_out, # au
                       Mp_in, # Mearth
                       Mp_out, # Mearth
                       ep_in,
                       ep_out):
        """
        https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.4575H/abstract
        
        Parameters
        ----------
        Ms: float
            Mass of the star (Msun).
        ap_in: float
            Semi-major axis (au) of inner planet.
        ap_out: float
            Semi-major axis (au) of outer planet.
        Mp_in: array
            Mass (Mearth) of inner planet.
        Mp_out: array
            Mass (Mearth) of outer planet.
        ep_in: array
            Eccentricity of inner planet.
        ep_out: array
            Eccentricity of outer planet.
        
        Returns
        -------
        Stability: bool
            If True, the planet pair is stable.
        """
        
        # Convert masses to kg
        tempMs = Ms*self.Msun # kg
        tempMp_in = Mp_in*self.Mearth # kg
        tempMp_out = Mp_out*self.Mearth # kg
        
        # Compute Hill radius
        HillRadius = (ap_in+ap_out)/2.*((tempMp_in+tempMp_out)/(3.*tempMs))**(1./3.)
        Delta = (ap_out*(1.-ep_out)-ap_in*(1.-ep_in))/HillRadius
        
        return Delta > self.Threshold
    
    def SummaryPlot(self,
                    PlanetDistribution,
                    MassModel,
                    EccentricityModel,
                    Ntest=100000,
                    block=True):
        """
        Parameters
        ----------
        PlanetDistribution: instance
            Instance of class PlanetDistribution.
        MassModel: instance
            Instance of class MassModel.
        EccentricityModel: instance
            Instance of class EccentricityModel.
        Ntest: int
            Number of test draws for summary plot.
        block: bool
            If True, blocks plots when showing.
        """
        
        Ntest = Ntest//10
        Rp_in = []
        Porb_in = []
        Rp_out = []
        Porb_out = []
        Trials_out = []
        for i in range(Ntest):
            Rp, Porb = PlanetDistribution.draw()
            Mp = MassModel.RadiusToMass(Rp)
            ep = EccentricityModel.getEccentricity(Porb)
            Rp_in += [Rp]
            Porb_in += [Porb]
            Rp, Porb, Mp, ep, Trials = self.draw(1.,
                                                 PlanetDistribution,
                                                 MassModel,
                                                 EccentricityModel,
                                                 returnTrials=True)
            Rp_out += [Rp]
            Porb_out += [Porb]
            Trials_out += [Trials]
        Rp_in = np.concatenate(Rp_in)
        Porb_in = np.concatenate(Porb_in)
        Rp_out = np.concatenate(Rp_out)
        Porb_out = np.concatenate(Porb_out)
        
        print('Weiss2018: %.2f (%.2f) Rearth median planet radius after (before) stabilization, %.2f (%.2f) d median planet orbital period after (before) stabilization, %.2f+-%.2f trials on average' % (np.median(Rp_out), np.median(Rp_in), np.median(Porb_out), np.median(Porb_in), np.mean(Trials_out), np.std(Trials_out)))
        
        Colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        f, ax = plt.subplots(1, 2)
        temp = np.hstack((Rp_in, Rp_out))
        ax[0].hist(Rp_in, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), color=Colors[0], label='Before')
        ax[0].hist(Rp_out, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), color=Colors[1], label='After', alpha=0.5)
        ax[0].set_xscale('log')
        ax[0].grid(axis='y')
        ax[0].set_xlabel('Planet radius [$R_\oplus$]')
        ax[0].set_ylabel('Number')
        ax[0].legend()
        temp = np.hstack((Porb_in, Porb_out))
        ax[1].hist(Porb_in, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), color=Colors[0], label='Before')
        ax[1].hist(Porb_out, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), color=Colors[1], label='After', alpha=0.5)
        ax[1].set_xscale('log')
        ax[1].grid(axis='y')
        ax[1].set_xlabel('Planet orbital period [d]')
        ax[1].set_ylabel('Number')
        ax[1].legend()
        plt.suptitle('Weiss2018')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show(block=block)
        
        pass
