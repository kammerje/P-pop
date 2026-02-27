#####
.. image:: logo.png
#####

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
