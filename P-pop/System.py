"""
# =============================================================================
# P-POP
# A Monte-Carlo tool to simulate exoplanet populations
# =============================================================================
"""


# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np
import scipy.optimize as so


# =============================================================================
# SYSTEM
# =============================================================================

class System():
    
    def __init__(self,
                 Star,
                 PlanetDistribution,
                 MassModel,
                 EccentricityModel,
                 OrbitModel,
                 AlbedoModel,
                 ExozodiModel,
                 StabilityModel=None,
                 Nstar=0,
                 Nuniverse=0,
                 Scale=1.):
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
        OrbitModel: instance
            Instance of class OrbitModel.
        AlbedoModel: instance
            Instance of class AlbedoModel.
        ExozodiModel: instance
            Instance of class ExozodiModel.
        StabilityModel: instance
            Instance of class StabilityModel.
        Nstar: int
            Number of the host star to which the system belongs to.
        Nuniverse: int
            Number of the universe to which the system belongs to.
        Scale: float
            Scaling factor for the planet occurrence rates.
        """
        
        # Constants.
        self.G = 6.674e-11 # m^3/kg/s^2
        self.Rsun = 695700000. # m
        self.Msun = 1.989e30 # kg
        self.au = 149597870700. # m
        self.sigma = 5.670e-8 # W/m^2/K^4
        
        self.Star = Star
        self.Nstar = Nstar
        self.Nuniverse = Nuniverse
        self.Scale = Scale
        
        if (StabilityModel is None):
            
            if (PlanetDistribution.returns == ['Rp', 'Porb']):
                
                # Draw planet radius and planet orbital period.
                self.Rp, self.Porb = PlanetDistribution.draw(Scale=self.Scale,
                                                             Star=self.Star) # Rearth, d
                
                # Sort planets by orbital period.
                ww = np.argsort(self.Porb)
                self.Rp = self.Rp[ww] # Rearth
                self.Porb = self.Porb[ww] # d
                
                # Draw planet mass and planet eccentricity.
                self.Mp = MassModel.RadiusToMass(self.Rp) # Mearth
                self.ep = EccentricityModel.getEccentricity(self.Porb)
            
            elif (PlanetDistribution.returns == ['Mp', 'Porb']):
                
                # Draw planet mass and planet orbital period.
                self.Mp, self.Porb = PlanetDistribution.draw(Scale=self.Scale,
                                                             Star=self.Star) # Mearth, d
                
                # Sort planets by orbital period.
                ww = np.argsort(self.Porb)
                self.Mp = self.Mp[ww] # Rearth
                self.Porb = self.Porb[ww] # d
                
                # Draw planet radius and planet eccentricity.
                self.Rp = MassModel.MassToRadius(self.Mp) # Rearth
                self.ep = EccentricityModel.getEccentricity(self.Porb)
            
            else:
                raise UserWarning()
        
        else:
            
            # Draw a stable system.
            self.Rp, self.Porb, self.Mp, self.ep = StabilityModel.draw(Star,
                                                                       PlanetDistribution,
                                                                       MassModel,
                                                                       EccentricityModel,
                                                                       Scale=self.Scale) # Rearth, d, Mearth
        
        # Draw planet orbits.
        self.ip, self.Omegap, self.omegap, self.thetap = OrbitModel.getOrbit(self.ep) # rad, rad, rad, rad
        
        # Draw planet albedos.
        self.Abond = AlbedoModel.getAbond(self.Rp)
        self.AgeomVIS = AlbedoModel.getAgeomVIS(self.Rp)
        self.AgeomMIR = AlbedoModel.getAgeomMIR(self.Rp)
        
        # Draw system exozodiacal dust level.
        self.z = ExozodiModel.getExozodiLevel()
        
        # Compute the other properties of the drawn planets.
        self.ap = self.getap() # au
        self.rp = self.getrp() # au
        self.rpproj = self.getrpproj() # au
        self.AngSep = self.getAngSep() # arcsec
        self.maxrpproj = self.getmaxrpproj() # au
        self.maxAngSep = self.getmaxAngSep() # arcsec
        self.Fp = self.getFp() # Searth
        self.fp = self.getfp()
        self.Tp = self.getTp() # K
        
        pass
    
    def getap(self):
        """
        Returns
        -------
        ap: array
            Semi-major axis (au) of drawn planets.
        """
        
        return ((self.G*self.Star.Mass*self.Msun*(self.Porb*86400.)**2)/(4.*np.pi**2))**(1./3.)/self.au
    
    def getrp(self):
        """
        Returns
        -------
        rp: array
            Physical separation (au) of drawn planets.
        """
        
        return self.ap*(1.-self.ep**2)/(1.+self.ep*np.cos(self.thetap))
    
    def getrpproj(self):
        """
        Returns
        -------
        rpproj: array
            Projected separation (au) of drawn planets.
        """
        
        return self.rp*np.sqrt(np.cos(self.omegap+self.thetap)**2+np.cos(self.ip)**2*np.sin(self.omegap+self.thetap)**2)
    
    def getAngSep(self):
        """
        Returns
        -------
        AngSep: array
            Projected angular separation (arcsec) of drawn planets.
        """
        
        return self.rpproj/self.Star.Dist
    
    def getmaxrpproj(self):
        """
        Returns
        -------
        maxrpproj: array
            Max projected separation (au) of drawn planets.
        """
        
        maxrpproj = [] # au
        
        for i in range(len(self.Rp)):
            
            # Function for computing the projected separation of the drawn
            # planets and its derivative with respect to the true anomaly of
            # the drawn planets.
            f = lambda theta: self.ap[i]*(1.-self.ep[i]**2)/(1.+self.ep[i]*np.cos(theta))*np.sqrt(np.cos(self.omegap[i]+theta)**2+np.cos(self.ip[i])**2*np.sin(self.omegap[i]+theta)**2)
            df = lambda theta: self.ap[i]*(1.-self.ep[i]**2)*self.ep[i]*np.sin(theta)/(1.+self.ep[i]*np.cos(theta))**2*np.sqrt(np.cos(self.omegap[i]+theta)**2+np.cos(self.ip[i])**2*np.sin(self.omegap[i]+theta)**2)+self.ap[i]*(1.-self.ep[i]**2)/(1.+self.ep[i]*np.cos(theta))*0.5*(np.cos(self.omegap[i]+theta)**2+np.cos(self.ip[i])**2*np.sin(self.omegap[i]+theta)**2)**(-1./2.)*(2.*np.cos(self.ip[i])**2-2.)*np.sin(self.omegap[i]+theta)*np.cos(self.omegap[i]+theta)
            
            # Maximize the above function based on different priors for the
            # true anomaly of the drawn planets.
            x0 = np.linspace(0., 2.*np.pi, 8, endpoint=False)
            xx = so.fsolve(df, x0)
            ff = f(xx)
            
            maxrpproj += [np.max(ff)] # au
        
        return np.array(maxrpproj)
    
    def getmaxAngSep(self):
        """
        Returns
        -------
        maxAngSep: array
            Max projected angular separation (arcsec) of drawn planets.
        """
        
        return self.maxrpproj/self.Star.Dist
    
    def getFp(self):
        """
        Returns
        -------
        Fp: array
            Host star flux (Searth) incident on drawn planets.
        """
        
        return self.sigma*self.Star.Teff**4*(self.Star.Rad*self.Rsun)**2/(self.rp*self.au)**2/1361.
    
    def getfp(self):
        """
        Returns
        -------
        fp: array
            Lambertian reflectance of drawn planets.
        """
        
        # Compute phase angle.
        alphap = np.arccos(-np.sin(self.ip)*np.sin(self.omegap+self.thetap)) # rad
        
        return np.abs((np.sin(alphap)+(np.pi-alphap)*np.cos(alphap))/np.pi)
    
    def getTp(self):
        """
        Returns
        -------
        Tp: array
            Equilibrium temperature (K) of drawn planets.
        """
        
        return ((self.Star.Rad*self.Rsun)**2*(1.-self.Abond)/(4.*(self.rp*self.au)**2))**(1./4.)*self.Star.Teff
    
    def write(self,
              Name):
        """
        Parameters
        ----------
        Name: str
            Name of the output planet table.
        """
        
        Table = open(Name+'.txt', 'w')
        
        # Old header.
        Table.write('nMC\tRp\tPorb\tMp\tecc\tinc\tOmega\tomega\ttheta\tAbond\tAgeomVIS\tAgeomMIR\tzodis\ta\trp\tang_sep\tang_sep_max\tFinc\tf\tTp\tnstar\tRs\tMs\tTs\tdist\tstype\tra\tdec\tlGal\tbGal\tWDSsep\tname\t\n')
        
        # New header.
        Table.write('Nuniverse\tRp\tPorb\tMp\tep\tip\tOmegap\tomegap\tthetap\tAbond\tAgeomVIS\tAgeomMIR\tz\tap\trp\tAngSep\tmaxAngSep\tFp\tfp\tTp\tNstar\tRs\tMs\tTs\tDs\tStype\tRA\tDec\tlGal\tbGal\tWDSsep\tname\t\n')
        
        Table.close()
        
        self.append(Name)
        
        pass
    
    def append(self,
               Name):
        """
        Parameters
        ----------
        Name: str
            Name of the output planet table.
        
        ..Notes::
            Table of planets with the following columns:
            - Nuniverse: number of the universe to which the planet belongs to.
            - Rp: planet radius (Rearth).
            - Porb: planet orbital period (d).
            - Mp: planet mass (Mearth).
            - ep: planet eccentricity.
            - ip: planet inclination (rad).
            - Omegap: planet longitude of ascending node (rad).
            - omegap: planet argument of periapsis (rad).
            - thetap: planet true anomaly (rad).
            - Abond: planet Bond albedo.
            - AgeomVIS: planet geometric albedo in the visible.
            - AgeomMIR: planet geometric albedo in the mid-infrared.
            - z: exozodiacal dust level.
            - ap: planet semi-major axis (au).
            - rp: planet physical separation (au).
            - AngSep: planet projected angular separation (arcsec).
            - maxAngSep: max planet projected angular separation (arcsec).
            - Fp: planet incident host star flux (Searth).
            - fp: planet Lambertian reflectance.
            - Tp: planet equilibrium temperature (K).
            - Nstar: number of the star.
            - Rs: host star radius (Rsun).
            - Ms: host star mass (Msun).
            - Ts: host star effective temperature (K).
            - Ds: host star distance (pc).
            - Stype: host star spectral type.
            - RA: host star right ascension (deg).
            - Dec: host star declination (deg).
            - lGal: Galactic longitude (deg).
            - bGal: Galactic latitude (deg).
            - name: name of the star.
        """
        
        Table = open(Name+'.txt', 'a')
        
        # Write the simulated planets to the planet population table.
        if (self.Star.lGal is None or self.Star.bGal is None):
            for i in range(len(self.Rp)):
                Table.write('%.0f\t' % self.Nuniverse
                            +'%.5f\t' % self.Rp[i]
                            +'%.5f\t' % self.Porb[i]
                            +'%.5f\t' % self.Mp[i]
                            +'%.5f\t' % self.ep[i]
                            +'%.5f\t' % self.ip[i]
                            +'%.5f\t' % self.Omegap[i]
                            +'%.5f\t' % self.omegap[i]
                            +'%.5f\t' % self.thetap[i]
                            +'%.5f\t' % self.Abond[i]
                            +'%.5f\t' % self.AgeomVIS[i]
                            +'%.5f\t' % self.AgeomMIR[i]
                            +'%.5f\t' % self.z
                            +'%.5f\t' % self.ap[i]
                            +'%.5f\t' % self.rp[i]
                            +'%.5f\t' % self.AngSep[i]
                            +'%.5f\t' % self.maxAngSep[i]
                            +'%.5f\t' % self.Fp[i]
                            +'%.5f\t' % self.fp[i]
                            +'%.5f\t' % self.Tp[i]
                            +'%.0f\t' % self.Nstar
                            +'%.5f\t' % self.Star.Rad
                            +'%.5f\t' % self.Star.Mass
                            +'%.5f\t' % self.Star.Teff
                            +'%.5f\t' % self.Star.Dist
                            +self.Star.Stype
                            +'\t%.5f\t' % self.Star.RA
                            +'%.5f\t' % self.Star.Dec
                            +'None\t'
                            +'None\t'
                            +'None\t'
                            +''.join(self.Star.Name.split())
                            +'\n')
        else:
            for i in range(len(self.Rp)):
                Table.write('%.0f\t' % self.Nuniverse
                            +'%.5f\t' % self.Rp[i]
                            +'%.5f\t' % self.Porb[i]
                            +'%.5f\t' % self.Mp[i]
                            +'%.5f\t' % self.ep[i]
                            +'%.5f\t' % self.ip[i]
                            +'%.5f\t' % self.Omegap[i]
                            +'%.5f\t' % self.omegap[i]
                            +'%.5f\t' % self.thetap[i]
                            +'%.5f\t' % self.Abond[i]
                            +'%.5f\t' % self.AgeomVIS[i]
                            +'%.5f\t' % self.AgeomMIR[i]
                            +'%.5f\t' % self.z
                            +'%.5f\t' % self.ap[i]
                            +'%.5f\t' % self.rp[i]
                            +'%.5f\t' % self.AngSep[i]
                            +'%.5f\t' % self.maxAngSep[i]
                            +'%.5f\t' % self.Fp[i]
                            +'%.5f\t' % self.fp[i]
                            +'%.5f\t' % self.Tp[i]
                            +'%.0f\t' % self.Nstar
                            +'%.5f\t' % self.Star.Rad
                            +'%.5f\t' % self.Star.Mass
                            +'%.5f\t' % self.Star.Teff
                            +'%.5f\t' % self.Star.Dist
                            +self.Star.Stype
                            +'\t%.5f\t' % self.Star.RA
                            +'%.5f\t' % self.Star.Dec
                            +'%.5f\t' % self.Star.lGal
                            +'%.5f\t' % self.Star.bGal
                            +'%.5f\t' % self.Star.WDSsep
                            +''.join(self.Star.Name.split())
                            +'\n')
        
        Table.close()
        
        pass
