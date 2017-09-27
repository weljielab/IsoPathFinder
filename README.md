# IsoPathFinder
Workflow for isotope mapping


# See below for how to run the program through iPython (Ensure you have ipython installed (we recommend the Anaconda suite, at https://www.anaconda.com/download/)):

To label the metabolic network, run the ReactionModeling_IsotopeTracing_Rhoades_NoNameConvert.ipynb notebook (e.g. if in the terminal, call 'ipython', then %run ReactionModeling_IsotopeTracing_Rhoades_NoNameConvert.ipynb'. NOTE, if you want to change the tracer and model organism, you will need to edit the appropriate cells in the .ipynb, which could be done by opening the notebook through Jupyter.

Once the labeled network is created (which will be stored locally), the functions for the IsoPathFinder can be loaded through iPython by '%run IsoPathFinder-Parallel_WithNameConversion_Rhoades.ipynb'. Once the functions are loaded, the BuildPaths function will allow you to search for an isotopologue at a desired path length (up to 15). You can do this by calling e.g. BuildPaths('Serine M+3',10) from the IPython interpreter. You will get a lot of warnings, thats ok, it should still run fine. Check the local directory to see that the csv result is written.

One could run the scripts purely out of a Jupyter Notebook as well, which we recommend if one wants to change the starting tracer or load their own model organism .dat from MetaCyc.
