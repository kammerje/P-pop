"""
# =============================================================================
# P-POP
# A Monte-Carlo tool to simulate exoplanet populations
# =============================================================================
"""


# =============================================================================
# IMPORTS
# =============================================================================

import os
import sys
import time

import Star, System


# =============================================================================
# SYSTEMGENERATOR
# =============================================================================

class SystemGenerator():
    
    def __init__(self,
                 StarCatalog,
                 StypeToModel,
                 ScalingModel,
                 MassModel,
                 EccentricityModel,
                 StabilityModel,
                 OrbitModel,
                 AlbedoModel,
                 ExozodiModel,
                 Stypes,
                 Dist_range, # pc
                 Dec_range, # deg
                 Teff_range, # K
                 Scenario,
                 SummaryPlots,
                 Ntest,
                 FigDir,
                 block):
        """
        Parameters
        ----------
        StarCatalog: module
            Module of type StarCatalog.
        StypeToModel: dict
            Dictionary mapping the spectral types to planet distribution
            modules.
        ScalingModel: module
            Module of type ScalingModel.
        MassModel: module
            Module of type MassModel.
        EccentricityModel: module
            Module of type EccentricityModel.
        StabilityModel: None, module
            Module of type StabilityModel.
        OrbitModel: module
            Module of type OrbitModel.
        AlbedoModel: module
            Module of type AlbedoModel.
        ExozodiModel: module
            Module of type ExozodiModel.
        Stypes: list
            Spectral types to be included.
        Dist_range: list
            Distance range (pc) to be included.
        Dec_range: list
            Declination range (deg) to be included.
        Teff_range: list
            Effective temperature range (K) to be included.
        Scenario: 'baseline', 'pessimistic', 'optimistic', 'mc'
            Scenario for planet occurrence rates.
        SummaryPlots: bool
            If True, makes summary plots after importing a module.
        Ntest: int
            Number of test draws for summary plot.
        FigDir: str, None
            Directory to which summary plots are saved.
        block: bool
            If True, blocks plots when showing.
        """
        
        # Get scenario.
        if ((Scenario == 'baseline') or (Scenario == 'optimistic') or (Scenario == 'pessimistic') or (Scenario == 'mc')):
            self.Scenario = Scenario
        else:
            print('--> WARNING: '+str(Scenario)+' is an unknown scenario')
            self.Scenario = 'baseline'
        print('--> Using scenario '+str(self.Scenario))
        
        # Create directory to which summary plots are saved.
        if (FigDir is not None and not os.path.exists(FigDir)):
            os.makedirs(FigDir)
        
        # Flag indicating whether the output planet population table has
        # already been created.
        self.TableFlag = False
        
        # Get star catalog.
        self.StarCatalog = self.getStarCatalog(StarCatalog,
                                               Stypes,
                                               Dist_range,
                                               Dec_range,
                                               Teff_range)
        if (SummaryPlots == True):
            self.StarCatalog.SummaryPlot(FigDir=FigDir,
                                         block=block)
        
        # Get planet distributions.
        self.PlanetDistributions = self.getPlanetDistributions(StypeToModel,
                                                               Scenario)
        if (SummaryPlots == True):
            # temp = sorted(set(val for val in self.PlanetDistributions.values()))
            temp = set(val for val in self.PlanetDistributions.values())
            for val in temp:
                val.SummaryPlot(Ntest=Ntest,
                                FigDir=FigDir,
                                block=block)
        
        # Get scaling model.
        self.ScalingModel = self.getScalingModel(ScalingModel)
        if (SummaryPlots == True and self.ScalingModel is not None):
            self.ScalingModel.SummaryPlot(Ntest=Ntest,
                                          FigDir=FigDir,
                                          block=block)
        
        # Get mass model.
        self.MassModel = self.getMassModel(MassModel)
        if (SummaryPlots == True):
            self.MassModel.SummaryPlot(Ntest=Ntest,
                                       FigDir=FigDir,
                                       block=block)
        
        # Get eccentricity model.
        self.EccentricityModel = self.getEccentricityModel(EccentricityModel)
        if (SummaryPlots == True):
            self.EccentricityModel.SummaryPlot(Ntest=Ntest,
                                               FigDir=FigDir,
                                               block=block)
        
        # Get stability model.
        self.StabilityModel = self.getStabilityModel(StabilityModel)
        if (SummaryPlots == True and self.StabilityModel is not None):
            try:
                self.StabilityModel.SummaryPlot(self.PlanetDistributions['G'],
                                                self.MassModel,
                                                self.EccentricityModel,
                                                Ntest=Ntest,
                                                FigDir=FigDir,
                                                block=block)
            except:
                self.StabilityModel.SummaryPlot(self.PlanetDistributions[sorted(self.PlanetDistributions.keys())[0]],
                                                self.MassModel,
                                                self.EccentricityModel,
                                                Ntest=Ntest,
                                                FigDir=FigDir,
                                                block=block)
        
        # Get orbit model.
        self.OrbitModel = self.getOrbitModel(OrbitModel)
        if (SummaryPlots == True):
            self.OrbitModel.SummaryPlot(self.EccentricityModel,
                                        Ntest=Ntest,
                                        FigDir=FigDir,
                                        block=block)
        
        # Get albedo model.
        self.AlbedoModel = self.getAlbedoModel(AlbedoModel)
        if (SummaryPlots == True):
            self.AlbedoModel.SummaryPlot(Ntest=Ntest,
                                         FigDir=FigDir,
                                         block=block)
        
        # Get exozodiacal dust model.
        self.ExozodiModel = self.getExozodiModel(ExozodiModel,
                                                 Scenario)
        if (SummaryPlots == True):
            self.ExozodiModel.SummaryPlot(Ntest=Ntest,
                                          FigDir=FigDir,
                                          block=block)
        
        pass
    
    def getStarCatalog(self,
                       StarCatalog,
                       Stypes,
                       Dist_range, # pc
                       Dec_range, # deg
                       Teff_range): # K
        
        return StarCatalog.StarCatalog(Stypes,
                                       Dist_range,
                                       Dec_range,
                                       Teff_range)
    
    def getPlanetDistributions(self,
                               StypeToModel,
                               Scenario):
        """
        Parameters
        ----------
        StypeToModel: dict
            Dictionary mapping the spectral types to planet distribution
            modules.
        Scenario: 'baseline', 'pessimistic', 'optimistic', 'mc'
            Scenario for planet occurrence rates.
        
        Returns
        -------
        PlanetDistributions: dict
            Dictionary mapping the spectral types to planet distribution
            instances.
        """
        
        # Get one planet distribution instance for each spectral type. If
        # multiple spectral types have the same planet distribution don't
        # create a new instance for each of them, but simply map it to the
        # already existing one.
        PlanetDistributions = {}
        for i, key_i in enumerate(StypeToModel):
            if (i == 0):
                PlanetDistributions[key_i] = StypeToModel[key_i].PlanetDistribution(Scenario)
            else:
                isnew = True
                for j, key_j in enumerate(StypeToModel):
                    if (j < i and StypeToModel[key_i] == StypeToModel[key_j]):
                        PlanetDistributions[key_i] = PlanetDistributions[key_j]
                        isnew = False
                if (isnew == True):
                    PlanetDistributions[key_i] = StypeToModel[key_i].PlanetDistribution(Scenario)
        
        return PlanetDistributions
    
    def getScalingModel(self,
                        ScalingModel):
        """
        Parameters
        ----------
        ScalingModel: None, module
            Module of type ScalingModel.
        
        Returns
        -------
        ScalingModel: None, instance
            Instance of class ScalingModel.
        """
        
        if (ScalingModel is None):
            return None
        else:
            return ScalingModel.ScalingModel()
    
    def getMassModel(self,
                     MassModel):
        """
        Parameters
        ----------
        MassModel: module
            Module of type MassModel.
        
        Returns
        -------
        MassModel: instance
            Instance of class MassModel.
        """
        
        return MassModel.MassModel()
    
    def getEccentricityModel(self,
                             EccentricityModel):
        """
        Parameters
        ----------
        EccentricityModel: module
            Module of type EccentricityModel.
        
        Returns
        -------
        EccentricityModel: instance
            Instance of class EccentricityModel.
        """
        
        return EccentricityModel.EccentricityModel()
    
    def getStabilityModel(self,
                          StabilityModel):
        """
        Parameters
        ----------
        StabilityModel: None, module
            Module of type StabilityModel.
        
        Returns
        -------
        StabilityModel: None, instance
            Instance of class StabilityModel.
        """
        
        if (StabilityModel is None):
            return None
        else:
            return StabilityModel.StabilityModel()
    
    def getOrbitModel(self,
                      OrbitModel):
        """
        Parameters
        ----------
        OrbitModel: module
            Module of type OrbitModel.
        
        Returns
        -------
        OrbitModel: instance
            Instance of class OrbitModel.
        """
        
        return OrbitModel.OrbitModel()
    
    def getAlbedoModel(self,
                       AlbedoModel):
        """
        Parameters
        ----------
        AlbedoModel: module
            Module of type AlbedoModel.
        
        Returns
        -------
        AlbedoModel: instance
            Instance of class AlbedoModel.
        """
        
        return AlbedoModel.AlbedoModel()
    
    def getExozodiModel(self,
                       ExozodiModel,
                       Scenario):
        """
        Parameters
        ----------
        ExozodiModel: module
            Module of type ExozodiModel.
        Scenario: 'baseline', 'pessimistic', 'optimistic', 'mc'
            Scenario for exozodi level.
        
        Returns
        -------
        ExozodiModel: instance
            Instance of class ExozodiModel.
        """
        
        return ExozodiModel.ExozodiModel(Scenario)
    
    def SimulateUniverses(self,
                          Name,
                          Nuniverses=1):
        """
        Parameters
        ----------
        Name: str
            Name of the output planet table.
        Nuniverses: int
            Number of universes to be simulated.
        """
        
        # Print.
        print('Simulating %.0f universe(s)' % Nuniverses)
        
        t0 = time.time()
        
        # Go through the number of universes.
        Nstars = len(self.StarCatalog.SC)
        PDs = list(set(self.PlanetDistributions.values()))
        for i in range(Nuniverses):
            
            # Re-draw occurrence rate parameters for 'mc' scenario.
            if (self.Scenario == 'mc'):
                for j in range(len(PDs)):
                    PDs[j].redraw_params()
            
            for j in range(Nstars):
                
                # Get star.
                self.Star = self.getStar(j)
                
                # Apply scaling for the planet occurrence rates.
                if (self.ScalingModel is None):
                    Scale = 1.
                else:
                    Scale = self.ScalingModel.getScale(self.Star)
                
                # Get system.
                self.System = self.getSystem(j,
                                             i,
                                             Scale)
                
                # If there is a planet distribution for the star's spectral
                # type (i.e. if the system is not None), create a new
                # planet population table (if it hasn't already been created)
                # and write the simulated planets to it.
                if (self.System is not None):
                    if (self.TableFlag == False):
                        self.System.write(Name)
                        self.TableFlag = True
                    else:
                        self.System.append(Name)
            
            sys.stdout.write('\r--> Universe %.0f of %.0f' % ((i+1), Nuniverses))
            sys.stdout.flush()
        print('')
        
        # # Go through the star catalog.
        # Nstars = len(self.StarCatalog.SC)
        # for i in range(Nstars):
            
        #     # Get star.
        #     self.Star = self.getStar(i)
            
        #     # Apply scaling for the planet occurrence rates.
        #     if (self.ScalingModel is None):
        #         Scale = 1.
        #     else:
        #         Scale = self.ScalingModel.getScale(self.Star)
            
        #     # Go through the number of universes to be simulated.
        #     for j in range(Nuniverses):
                
        #         # Get system.
        #         self.System = self.getSystem(i,
        #                                      j,
        #                                      Scale)
                
        #         # If there is a planet distribution for the star's spectral
        #         # type (i.e. if the system is not None), create a new
        #         # planet population table (if it hasn't already been created)
        #         # and write the simulated planets to it.
        #         if (self.System is not None):
        #             if (self.TableFlag == False):
        #                 self.System.write(Name)
        #                 self.TableFlag = True
        #             else:
        #                 self.System.append(Name)
                
        #     sys.stdout.write('\r--> Star %.0f of %.0f, scaling = %.1f' % ((i+1), Nstars, Scale))
        #     sys.stdout.flush()
        # print('')
        
        t1 = time.time()
        
        # Print.
        if (self.StabilityModel is None):
            print('--> Finished after %.0f s' % (t1-t0))
        else:
            print('--> Finished after %.0f s, drawing stable system failed %.0f times' % (t1-t0, self.StabilityModel.Nfails))
        
        pass
    
    def getStar(self,
                index):
        """
        Parameters
        ----------
        index: int
            Index of requested star in the star catalog.
        
        Returns
        -------
        Star: instance
            Instance of class Star.
        """
        
        try:
            tempVmag = self.StarCatalog.SC[index]['Vmag']
        except:
            tempVmag = None
        try:
            tempJmag = self.StarCatalog.SC[index]['Jmag']
        except:
            tempJmag = None
        try:
            tempHmag = self.StarCatalog.SC[index]['Hmag']
        except:
            tempHmag = None
        try:
            tempWDSsep = self.StarCatalog.SC[index]['WDSsep']
        except:
            tempWDSsep = None
        try:
            tempWDSdmag = self.StarCatalog.SC[index]['WDSdmag']
        except:
            tempWDSdmag = None
        try:
            templGal = self.StarCatalog.SC[index]['lGal']
        except:
            templGal = None
        try:
            tempbGal = self.StarCatalog.SC[index]['bGal']
        except:
            tempbGal = None
        
        return Star.Star(self.StarCatalog.SC[index]['Name'],
                         self.StarCatalog.SC[index]['Dist'],
                         self.StarCatalog.SC[index]['Stype'],
                         self.StarCatalog.SC[index]['Rad'],
                         self.StarCatalog.SC[index]['Teff'],
                         self.StarCatalog.SC[index]['Mass'],
                         self.StarCatalog.SC[index]['RA'],
                         self.StarCatalog.SC[index]['Dec'],
                         tempVmag,
                         tempJmag,
                         tempHmag,
                         tempWDSsep,
                         tempWDSdmag,
                         templGal,
                         tempbGal)
    
    def getSystem(self,
                  Nstar,
                  Nuniverse,
                  Scale):
        """
        Parameters
        ----------
        Nstar: int
            Number of the host star to which the system belongs to.
        Nuniverse: int
            Number of the universe to which the system belongs to.
        Scale: float
            Scaling factor for the planet occurrence rates.
        
        Returns
        -------
        System: None, instance
            Instance of class System.
        """
        
        # If there is no planet distribution for the star's spectral type,
        # return None.
        if (self.Star.Stype in self.PlanetDistributions.keys()):
            return System.System(self.Star,
                                 self.PlanetDistributions[self.Star.Stype],
                                 self.MassModel,
                                 self.EccentricityModel,
                                 self.OrbitModel,
                                 self.AlbedoModel,
                                 self.ExozodiModel,
                                 self.StabilityModel,
                                 Nstar,
                                 Nuniverse,
                                 Scale)
        else:
            return None
