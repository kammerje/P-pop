"""
# =============================================================================
# P-pop
# A Monte-Carlo tool to simulate exoplanet populations
#
# Authors: Jens Kammerer, Sascha Quanz, Emile Fontanet
# Version: 1.0.0
# Last edited: 26.02.2026
# =============================================================================
#
# P-pop is introduced in Kammerer & Quanz 2018
# (https://ui.adsabs.harvard.edu/abs/2018A%26A...609A...4K/abstract). Please
# cite this paper if you use P-pop for your research.
#
# P-pop makes use of forecaster from Chen & Kipping 2017
# (https://ui.adsabs.harvard.edu/abs/2017ApJ...834...17C/abstract).
"""


# Don't print annoying warnings. Comment out if you want to see them.
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# IMPORTS
# =============================================================================

# Import catalogs, distributions, and models here.
import SystemGenerator
from StarCatalogs import alphaCenA,\
                         CrossfieldBrightSample,\
                         ExoCat1,\
                         LTC2,\
                         LTC3,\
                         LTC4,\
                         HPIC_LTC4_combined
from PlanetDistributions import Fressin2013,\
                                Burke2015,\
                                Dressing2015,\
                                SAG13,\
                                SAG13_extrap,\
                                Weiss2018,\
                                Weiss2018KDE,\
                                HabitableNominal,\
                                HabitablePessimistic,\
                                Fernandes2019Symm,\
                                Bryson2021Model1Hab2Low,\
                                Bryson2021Model1Hab2High,\
                                Bergsten2022
from ScalingModels import BinarySuppression
from MassModels import Chen2017
from EccentricityModels import Circular
from StabilityModels import He2019
from OrbitModels import Random
from AlbedoModels import Uniform,\
                         Constant
from ExozodiModels import Ertel2018,\
                          Ertel2020,\
                          Median


# =============================================================================
# SETUP
# =============================================================================

# Select the star catalog, the spectral types, the distance range, the
# declination range, and the effective temperature range which should be
# included here.
# StarCatalog = CrossfieldBrightSample # used in Kammerer & Quanz 2018
# StarCatalog = ExoCat1 # used by NASA
# StarCatalog = LTC2 # LIFE Target Catalog (version 2)
# StarCatalog = LTC3 # LIFE Target Catalog (version 3)
StarCatalog = LTC4 # LIFE Target Catalog (version 4)
Stypes = ['A', 'F', 'G', 'K', 'M'] # list of str
Dist_range = [0., 20.] # pc, list of float, [min, max]
Dec_range = [-90., 90.] # deg, list of float, [min, max]
# Teff_range = [4800., 6300.] # K, list of float, [min, max]
Teff_range = None

# Select the planet distributions, the scenario, and the scaling model which
# should be used here. A different planet distribution can be assigned to each
# spectral type.
# dict, StarCatalog.Stype as keys, PlanetDistribution as data
# StypeToModel = {'A': Fressin2013, 'F': Fressin2013, 'G': Fressin2013, 'K': Fressin2013, 'M': Fressin2013}
# StypeToModel = {'F': Burke2015, 'G': Burke2015, 'K': Burke2015, 'M': Dressing2015}
# StypeToModel = {'F': SAG13, 'G': SAG13, 'K': SAG13, 'M': Dressing2015} # used in the ESA Voyage 2050 White Paper
# StypeToModel = {'F': SAG13, 'G': SAG13, 'K': SAG13, 'M': SAG13}
# StypeToModel = {'F': Weiss2018, 'G': Weiss2018, 'K': Weiss2018, 'M': Weiss2018}
# StypeToModel = {'F': Weiss2018KDE, 'G': Weiss2018KDE, 'K': Weiss2018KDE, 'M': Weiss2018KDE}
# StypeToModel = {'A': SAG13, 'F': SAG13, 'G': SAG13, 'K': SAG13, 'M': Dressing2015}
StypeToModel = {'A': SAG13, 'F': SAG13, 'G': SAG13, 'K': SAG13, 'M': SAG13}
# StypeToModel = {'A': HabitableNominal, 'F': HabitableNominal, 'G': HabitableNominal, 'K': HabitableNominal, 'M': HabitableNominal}
# StypeToModel = {'A': HabitablePessimistic, 'F': HabitablePessimistic, 'G': HabitablePessimistic, 'K': HabitablePessimistic, 'M': HabitablePessimistic}
# StypeToModel = {'A': Fernandes2019Symm, 'F': Fernandes2019Symm, 'G': Fernandes2019Symm, 'K': Fernandes2019Symm, 'M': Fernandes2019Symm}
# StypeToModel = {'A': Bryson2021Model1Hab2Low, 'F': Bryson2021Model1Hab2Low, 'G': Bryson2021Model1Hab2Low, 'K': Bryson2021Model1Hab2Low, 'M': Bryson2021Model1Hab2Low}
# StypeToModel = {'A': Bryson2021Model1Hab2High, 'F': Bryson2021Model1Hab2High, 'G': Bryson2021Model1Hab2High, 'K': Bryson2021Model1Hab2High, 'M': Bryson2021Model1Hab2High}
# StypeToModel = {'A': Bergsten2022, 'F': Bergsten2022, 'G': Bergsten2022, 'K': Bergsten2022, 'M': Bergsten2022}
Scenario = 'baseline'
# Scenario = 'pessimistic' # for 1-sigma lower error bars on planet distribution
# Scenario = 'optimistic' # for 1-sigma upper error bars on planet distribution
# Scenario = 'mc' # draw parameters from MCMC posterior for each universe
# ScalingModel = None
ScalingModel = BinarySuppression

# Select the mass model, the eccentricity model, the stability model, the orbit
# model, the albedo model, and the exozodiacal dust model which should be used
# here.
MassModel = Chen2017 # Forecaster
EccentricityModel = Circular
# StabilityModel = None
StabilityModel = He2019
OrbitModel = Random
AlbedoModel = Uniform
# AlbedoModel = Constant
# ExozodiModel = Ertel2018
ExozodiModel = Ertel2020
# ExozodiModel = Median

# Select whether you want to generate summary plots after loading the catalogs,
# distributions, and models selected above, how many test draws should be
# done for generating these plots, and where you want to save them.
SummaryPlots = True
# SummaryPlots = False
Ntest = 100000 # int
# Ntest = 10000 # int
# FigDir = None # if you don't want to save the summary plots
FigDir = 'Figures/' # str, should end with a slash ("/")
# block = True # interrupt script and display summary plots
block = False # don't interrupt script, only save summary plots

# Select a name for the output planet population table and how many universes
# should be simulated.
Name = 'TestPlanetPopulation' # str
Nuniverses = 10 # int


# =============================================================================
# P-POP
# =============================================================================

# Don't modify the following code.
SysGen = SystemGenerator.SystemGenerator(StarCatalog,
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
                                         block)
SysGen.SimulateUniverses(Name,
                         Nuniverses)
