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

Package content
***************

Available planet distributions:

- :code:`Bergsten2022`: From `Bergsten et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022AJ....164..190B/abstract>`__
