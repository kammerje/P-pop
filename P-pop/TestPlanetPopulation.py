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

import ReadPlanetPopulation as RPP


# =============================================================================
# SETUP
# =============================================================================

# Select the name of the planet population table to be read.
PathPlanetTable = 'TestPlanetPopulation.txt' # str

# Select a model for the computation of the habitable zone.
Model = 'MS'
#Model = 'POST-MS'


# =============================================================================
# READ PLANET POPULATION
# =============================================================================

# The next five lines read the planet population table and its photometry.
PP = RPP.PlanetPopulation(PathPlanetTable)
PP.ComputeHZ(Model)

print('Number of planets in the planet population table')
print(len(PP.Rp))

print('Planet radius (Rearth) of the first planet in the planet population table')
print(PP.Rp[0])

print('Host star radius (Rsun) of the first planet in the planet population table')
print(PP.Rs[0])

import pdb; pdb.set_trace()
