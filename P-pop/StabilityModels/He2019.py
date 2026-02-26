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

import Star


# =============================================================================
# HE2019
# =============================================================================

class StabilityModel():
    """
    https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.4575H/abstract
    """
    
    def __init__(self):
        
        # Max number of re-draws if system is unstable and number of times no
        # stable system could be found after self.MaxTrials trials.
        self.MaxTrials = 100.
        self.Nfails = 0
        
        # Mutual Hill radius above which planets are considered to be stable.
        self.Threshold = 8.
        
        # Constants.
        self.yr2d = 365.24
        self.Mearth = 5.972e24 # kg
        self.Msun = 1.989e30 # kg
        
        pass
    
    def draw(self,
             Star,
             PlanetDistribution,
             MassModel,
             EccentricityModel,
             Nplanets=None,
             Scale=1.,
             returnTrials=False):
        """
        Parameters
        ----------
        Star: instance
            Instance of class Star.
        PlanetDistribution: instance
            Instance of class PlanetDistribution.
        MassModel: instance
            Instance of class MassModel.
        EccentricityModel: instance
            Instance of class EccentricityModel.
        Nplanets: None, int
            Number of planets to be drawn.
        Scale: float
            Scaling factor for the planet occurrence rates.
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
        
        # Try to draw a stable system with the correct number of planets.
        if (returnTrials == False):
            Rp, Porb, Mp, ep = self.draw_temp(Star=Star,
                                              PlanetDistribution=PlanetDistribution,
                                              MassModel=MassModel,
                                              EccentricityModel=EccentricityModel,
                                              Nplanets=Nplanets,
                                              Scale=Scale,
                                              returnTrials=returnTrials)
        else:
            Rp, Porb, Mp, ep, Trials = self.draw_temp(Star=Star,
                                                      PlanetDistribution=PlanetDistribution,
                                                      MassModel=MassModel,
                                                      EccentricityModel=EccentricityModel,
                                                      Nplanets=Nplanets,
                                                      Scale=Scale,
                                                      returnTrials=returnTrials)
        
        # If no stable system could be drawn, try to draw a stable system with
        # a new number of planets.
        while(Rp is None):
            self.Nfails += 1
            if (returnTrials == False):
                Rp, Porb, Mp, ep = self.draw_temp(Star=Star,
                                                  PlanetDistribution=PlanetDistribution,
                                                  MassModel=MassModel,
                                                  EccentricityModel=EccentricityModel,
                                                  Scale=Scale,
                                                  returnTrials=returnTrials)
            else:
                Rp, Porb, Mp, ep, Trials = self.draw_temp(Star=Star,
                                                          PlanetDistribution=PlanetDistribution,
                                                          MassModel=MassModel,
                                                          EccentricityModel=EccentricityModel,
                                                          Scale=Scale,
                                                          returnTrials=returnTrials)
        
        if (returnTrials == False):
            return Rp, Porb, Mp, ep
        else:
            return Rp, Porb, Mp, ep, Trials
    
    def draw_temp(self,
                  Star,
                  PlanetDistribution,
                  MassModel,
                  EccentricityModel,
                  Nplanets=None,
                  Scale=1.,
                  returnTrials=False):
        """
        Parameters
        ----------
        Star: instance
            Instance of class Star.
        PlanetDistribution: instance
            Instance of class PlanetDistribution.
        MassModel: instance
            Instance of class MassModel.
        EccentricityModel: instance
            Instance of class EccentricityModel.
        Nplanets: None, int
            Number of planets to be drawn.
        Scale: float
            Scaling factor for the planet occurrence rates.
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
        
        if (PlanetDistribution.returns == ['Rp', 'Porb']):
            
            # Draw planet radius and planet orbital period.
            Rp, Porb = PlanetDistribution.draw(Nplanets=Nplanets,
                                               Scale=Scale,
                                               Star=Star) # Rearth, d
            
            # Sort planets by orbital period.
            ww = np.argsort(Porb)
            Rp = Rp[ww] # Rearth
            Porb = Porb[ww] # d
            
            # Draw planet mass and planet eccentricity.
            Mp = MassModel.RadiusToMass(Rp) # Mearth
            ep = EccentricityModel.getEccentricity(Porb)
        
        elif (PlanetDistribution.returns == ['Mp', 'Porb']):
            
            # Draw planet mass and planet orbital period.
            Mp, Porb = PlanetDistribution.draw(Nplanets=Nplanets,
                                               Scale=Scale,
                                               Star=Star) # Mearth, d
            
            # Sort planets by orbital period.
            ww = np.argsort(Porb)
            Mp = Mp[ww] # Rearth
            Porb = Porb[ww] # d
            
            # Draw planet radius and planet eccentricity.
            Rp = MassModel.MassToRadius(Mp) # Rearth
            ep = EccentricityModel.getEccentricity(Porb)
        
        else:
            raise UserWarning()
        
        # The drawn system is a multi-planet system.
        Nplanets = len(Rp)
        Trials = 1
        if (Nplanets > 1):
            
            # Check system stability.
            Stability = self.CheckSystemStability(Star.Mass,
                                                  Porb,
                                                  Mp,
                                                  ep)
            
            # If the system is unstable and the max number of re-draws has not
            # been exceeded, draw a new system.
            while ((Stability == False) and (Trials < self.MaxTrials)):
                
                if (PlanetDistribution.returns == ['Rp', 'Porb']):
                    
                    # Draw planet radius and planet orbital period. Keep the same
                    # number of planets.
                    Rp, Porb = PlanetDistribution.draw(Nplanets=Nplanets,
                                                       Scale=Scale,
                                                       Star=Star) # Rearth, d
                    
                    # Sort planets by orbital period.
                    ww = np.argsort(Porb)
                    Rp = Rp[ww] # Rearth
                    Porb = Porb[ww] # d
                    
                    # Draw planet mass and planet eccentricity.
                    Mp = MassModel.RadiusToMass(Rp) # Mearth
                    ep = EccentricityModel.getEccentricity(Porb)
                
                elif (PlanetDistribution.returns == ['Mp', 'Porb']):
                    
                    # Draw planet mass and planet orbital period.
                    Mp, Porb = PlanetDistribution.draw(Nplanets=Nplanets,
                                                       Scale=Scale,
                                                       Star=Star) # Mearth, d
                    
                    # Sort planets by orbital period.
                    ww = np.argsort(Porb)
                    Mp = Mp[ww] # Rearth
                    Porb = Porb[ww] # d
                    
                    # Draw planet radius and planet eccentricity.
                    Rp = MassModel.MassToRadius(Mp) # Rearth
                    ep = EccentricityModel.getEccentricity(Porb)
                
                else:
                    raise UserWarning()
                
                # Check system stability.
                Stability = self.CheckSystemStability(Star.Mass,
                                                      Porb,
                                                      Mp,
                                                      ep)
                Trials += 1
        else:
            Stability = True
        
        if (Stability == True):
            if (returnTrials == False):
                return Rp, Porb, Mp, ep
            else:
                return Rp, Porb, Mp, ep, Trials
        else:
            if (returnTrials == False):
                return None, None, None, None
            else:
                return None, None, None, None, None
    
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
        
        # Compute semi-major axes in au.
        tempPorb = Porb/self.yr2d # yr
        ap = np.power(tempPorb, 2./3.)*Ms**(1./3.) # au
        
        # Check stability between each pair of planets.
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
        
        # Convert masses to kg.
        tempMs = Ms*self.Msun # kg
        tempMp_in = Mp_in*self.Mearth # kg
        tempMp_out = Mp_out*self.Mearth # kg
        
        # Compute Hill radius.
        HillRadius = (ap_in+ap_out)/2.*((tempMp_in+tempMp_out)/(3.*tempMs))**(1./3.)
        Delta = (ap_out*(1.-ep_out)-ap_in*(1.+ep_in))/HillRadius
        
        return Delta > self.Threshold
    
    def SummaryPlot(self,
                    PlanetDistribution,
                    MassModel,
                    EccentricityModel,
                    Ntest=100000,
                    FigDir=None,
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
        FigDir: str
            Directory to which summary plots are saved.
        block: bool
            If True, blocks plots when showing.
        """
        
        Ntest = Ntest//10
        Rp_in = []
        Porb_in = []
        Rp_out = []
        Porb_out = []
        Trials_out = []
        Sun = Star.Star('Sun',
                        10., # pc
                        'G',
                        1., # Rsun
                        5780., # K
                        1., # Msun
                        0., # deg
                        0.) # deg
        for i in range(Ntest):
            if (PlanetDistribution.returns == ['Rp', 'Porb']):
                Rp, Porb = PlanetDistribution.draw(Star=Sun)
                Mp = MassModel.RadiusToMass(Rp)
                ep = EccentricityModel.getEccentricity(Porb)
                Rp_in += [Rp]
                Porb_in += [Porb]
            elif (PlanetDistribution.returns == ['Mp', 'Porb']):
                Mp, Porb = PlanetDistribution.draw(Star=Sun)
                Rp = MassModel.MassToRadius(Mp)
                ep = EccentricityModel.getEccentricity(Porb)
                Rp_in += [Rp]
                Porb_in += [Porb]
            else:
                raise UserWarning()
            Rp, Porb, Mp, ep, Trials = self.draw(Sun,
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
        
        print('--> He2019:\n%.2f (%.2f) Rearth median planet radius after (before) stabilization\n%.2f (%.2f) d median planet orbital period after (before) stabilization\n%.2f+-%.2f trials on average' % (np.median(Rp_out), np.median(Rp_in), np.median(Porb_out), np.median(Porb_in), np.mean(Trials_out), np.std(Trials_out)))
        
        Colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        Weight = 1./len(Rp_in)
        f, ax = plt.subplots(1, 2)
        temp = np.hstack((Rp_in, Rp_out))
        ax[0].hist(Rp_in, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), weights=np.ones_like(Rp_in)*Weight, color=Colors[0], alpha=0.5, label='Before')
        ax[0].hist(Rp_out, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), weights=np.ones_like(Rp_out)*Weight, color=Colors[1], alpha=0.5, label='After')
        ax[0].set_xscale('log')
        ax[0].grid(axis='y')
        ax[0].set_xlabel('Planet radius [$R_\oplus$]')
        ax[0].set_ylabel('Fraction')
        ax[0].legend()
        temp = np.hstack((Porb_in, Porb_out))
        ax[1].hist(Porb_in, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), weights=np.ones_like(Porb_in)*Weight, color=Colors[0], alpha=0.5, label='Before')
        ax[1].hist(Porb_out, bins=np.logspace(np.log10(np.min(temp)), np.log10(np.max(temp)), 25), weights=np.ones_like(Porb_out)*Weight, color=Colors[1], alpha=0.5, label='After')
        ax[1].set_xscale('log')
        ax[1].grid(axis='y')
        ax[1].set_xlabel('Planet orbital period [d]')
        ax[1].set_ylabel('Fraction')
        ax[1].legend()
        plt.suptitle('He2019')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if (FigDir is not None):
            plt.savefig(FigDir+'StabilityModel.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
