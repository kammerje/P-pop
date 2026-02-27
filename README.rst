#####
P-pop
#####
.. image:: logo.png
P-pop is a Python tool for generating synthetic planet populations. A brief documentation can be found in :code:`P-pop/Documentation.pdf`.

Installation
************

Start by cloning the Git repository:

::

	git clone https://github.com/kammerje/P-pop.git

If you would like to install a specific branch (e.g., the :code:`develop` branch):

::

	git clone https://github.com/kammerje/P-pop.git@develop

It is **highly** recommended that you create a unique Conda environment to hold all of the P-pop dependencies:

::

	conda create -n p-pop python=3.9
	conda activate p-pop

With the Conda environment created, move to the cloned repository and install the dependencies and P-pop itself:

::

	cd where/you/saved/the/git/repo
	pip install -r requirements.txt
	pip install -e .

You should now be able to run P-pop:

::

	python P-pop/P-pop.py

This will simulate an example planet population called :code:`TestPlanetPopulation.txt`. You can then run the script :code:`P-pop/TestPlanetPopulation.py` as an example of how to read in and work with a P-pop planet population:

::

	python P-pop/TestPlanetPopulation.py

To customize your generated exoplanet population, adapt the script :code:`P-pop/P-pop.py` according to your wishes. Below you can find a list of the provided star catalogs, planet distributions, and parameter models.

Package content
***************

Available star catalogs:

- :code:`alphaCenA`: From `Crossfield (2013) <https://ui.adsabs.harvard.edu/abs/2013A%26A...551A..99C/abstract>`__, but using only the star alpha Cen A from that catalog. Serves as an example for how to select custom subsets of any of the provided star catalogs.
- :code:`ExoCat1`: From `Turnbull (2015) <https://ui.adsabs.harvard.edu/abs/2015arXiv151001731T/abstract>`__.
- :code:`HPIC_LTC4_combined`: Custom combination of HPIC and LTC4.
- :code:`HPIC`: From `Tuchow et al. (2024) <https://ui.adsabs.harvard.edu/abs/2024arXiv240208038T/abstract>`__.
- :code:`LTC2`: LIFE Target Catalog version 2.
- :code:`LTC3`: LIFE Target Catalog version 3 from `Quanz et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022A%26A...664A..21Q/abstract>`__.
- :code:`LTC4`: LIFE Target Catalog version 4 from `Menti et al. (2024) <https://ui.adsabs.harvard.edu/abs/2024RNAAS...8..267M/abstract>`__.

Available planet distributions:

- :code:`Bergsten2022`: From `Bergsten et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022AJ....164..190B/abstract>`__.
- :code:`Bryson2021Model1Hab2High`: From `Bryson et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021AJ....161...36B/abstract>`__. Requires unpublished parameter posterior file. Contact Steve Bryson directly to obtain it.
- :code:`Bryson2021Model1Hab2Low`: From `Bryson et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021AJ....161...36B/abstract>`__. Requires unpublished parameter posterior file. Contact Steve Bryson directly to obtain it.
- :code:`Burke2015`: From `Burke et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015ApJ...809....8B/abstract>`__.
- :code:`Dressing2015`: From `Dressing et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015ApJ...807...45D/abstract>`__.
- :code:`EarthTwin`: Simulate an Earth twin (1 Earth radius, 1 Earth insolation) around every star.
- :code:`Fernandes2019Symm`: From `Fernandes et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019ApJ...874...81F/abstract>`__.
- :code:`Fressin2013`: From `Fressin et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013ApJ...766...81F/abstract>`__.
- :code:`HabitableNominal`: Simulate a nominal population of habitable zone planets with 0.3 planets/star for AFGK-stars and 0.5 planets/star for M-dwarfs based on a log-normal distribution of planet radius from 0.5-1.5 Earth radii and of planet insolation from 0.35-1.75 Earth insolations.
- :code:`HabitablePessimistic`: Simulate a pessimistic population of habitable zone planets with 0.15 planets/star for AFGK-stars and 0.25 planets/star for M-dwarfs based on a log-normal distribution of planet radius from 0.5-1.5 Earth radii and of planet insolation from 0.35-1.75 Earth insolations.
- :code:`SAG13_extrap`: From `Kopparapu et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018ApJ...856..122K/abstract>`__. Extrapolated out to orbital periods of 20000 days.
- :code:`SAG13_rv`: From `Kopparapu et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018ApJ...856..122K/abstract>`__. Different implementation using the continuous random variable class from SciPy. Experimental!
- :code:`SAG13`: From `Kopparapu et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018ApJ...856..122K/abstract>`__.
- :code:`Weiss2018`: From `Kopparapu et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018ApJ...856..122K/abstract>`__, but ensuring the creation of peas-in-a-pod multi-planet systems according to `Weiss et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018AJ....155...48W/abstract>`__.
- :code:`Weiss2018KDE`: From `Kopparapu et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018ApJ...856..122K/abstract>`__, but ensuring the creation of peas-in-a-pod multi-planet systems according to `Weiss et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018AJ....155...48W/abstract>`__. Different implementation using a kernel density estimation. Experimental!

Available eccentricity models:

- :code:`Circular`: Place all planets on circular orbits.

Available orbit models:

- :code:`Quadrature`: Place all planets at maximum orbital elongation.
- :code:`Random`: Distribute planets randomly on the sphere.

Available mass models:
- :code:`Chen2017`: Use :code:`Forecaster` (`Chen & Kipping (2017) <https://ui.adsabs.harvard.edu/abs/2017ApJ...834...17C/abstract>`__) to forecast planet masses from planet radii and vice versa.
- :code:`EarthMass`: Assign a mass of 1 Earth mass to every planet.

Available albedo models:
- :code:`Constant`: Assign a constant Bond albedo of 0.4, geometric visible albedo of 0.3, and geometric mid-infrared albedo of 0.05 to every planet.
- :code:`EarthAlbedo`: Assign a constant Bond albedo of 0.306, geometric visible albedo of 0.434, and geometric mid-infrared albedo of 0.05 to every planet.
- :code:`Uniform`: Distribute the albedos uniformly with a Bond albedo in [0.0, 0.8), a geometric visible albedo in [0.0, 0.6), and a geometric mid-infrared albedo in [0.0, 0.1).

Available exozodi models:
- :code:`Ertel2018`: From `Ertel et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018AJ....155..194E/abstract>`__.
- :code:`Ertel2020`: From `Ertel et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020AJ....159..177E/abstract>`__.
- :code:`Median`: Assign an exozodi level of 3 zodi to every system.

Available stability models:
- :code:`He2019`: From `He et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.4575H/abstract>`__. Re-draw planets generated in multi-planet systems until mutually stable orbits are found. Note that this alters the original planet radius and orbital period distribution slightly. A summary plot showing the impact of this can be created.

Available scaling models:
- :code:`BinarySuppression`: Suppress the planet occurrence rate around <50 au binaries to 30% of its nominal value according to `Kraus et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016AJ....152....8K/abstract>`__.
