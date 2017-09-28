# IsoPathFinder

File information and workflow for simulating isotope tracing at atomic resolution and constructing paths from the downstream metabolite of interest back to the source isotope tracer.

# File descriptions


*CpdConvertDict.csv* - Matches MetaCyc-specific metabolite names to more commonly used names of major metabolites of interest (used in the IsoPathFinder script)

*IsoPathFinder-Parallel_WithNameConversion_Rhoades* - Python scripts to take precompiled labeled metabolic networks and trace paths to labeled metabolite of interest from the initial isotope tracer used to build the labeled network. If using the Jupyter notebook, this would be where specific metabolite and path lengths would be entered (see comments within the notebook for more detail). Runs in linux and can use multiple cores with the pathos package installed.

*PubChemCpdMatchAdd.csv* - SMILES strings matched to metabolite names, courtesy of the PubChemPy module.

*ReactionModeling_IsotopeTracing_Rhoades_NoNameConvert* - Python script to build the labeled metabolite network. If you are interested in simulating any tracer, the notebook can be modified accordingly (in the paper, we built a 13C6 glucose network). Detailed comments exist in the notebook to both change the labeled metabolite(s) to simulate, as well as exlusion of any uninformative metabolites (e.g. we 13C-carbon dioxide at each step). Runs in linux and can use multiple cores with the pathos package installed.

*ReactionTrimmer_Parallel_Rhoades* - Custom example script to filter the paths based on prior information. In the paper, thermodynamics and circadian rhythm expression datasets were incorporated to trim paths. Note we leave this here as a template, and trimming reactions should be rationalized based on the experimental objectives. Runs in linux and can use multiple cores with the pathos package installed.

*Tepper2013_GibbsECs.csv* - Estimated Gibbs free energies from Tepper et al., 2013, used in the ReactionTrimmer_Parallel_Rhoades script.

*Xu2011Gill2015_ECs_Cyclers.csv* - Significantly circadian enzyme expression from Xu et al., 2011 and Gill et al., 2015, used to trim paths without any enzyme cyclers through the ReactionTrimmer_Parallel_Rhoades script.

*\*atom-mappings-smiles\** - contains the atom-mapping reactions for fly, human, and mouse from MetaCyc. Note this would be read in the ReactionModeling_IsotopeTracing... script, and allows the user to specify which species to simulate isotope labeling. Other organisms can be pulled from MetaCyc as well, or trimming of reactions by organism from the total set of atom mapping solutions (atom-mappings-smilesAdd.dat).

*reaction-linksAdd.txt* - Links the Metacyc reaction names to ECs.


# Running the program through iPython:

Ensure you have ipython installed (we recommend the Anaconda suite, at https://www.anaconda.com/download/)

To label the metabolic network, run the ReactionModeling_IsotopeTracing_Rhoades_NoNameConvert.ipynb notebook (e.g. if in the terminal, call 'ipython', then %run ReactionModeling_IsotopeTracing_Rhoades_NoNameConvert.ipynb'. NOTE, if you want to change the tracer and model organism, you will need to edit the appropriate cells in the .ipynb, which could be done by opening the notebook through Jupyter (details comments within).

Once the labeled network is created (which will be stored locally), the functions for the IsoPathFinder can be loaded through iPython by '%run IsoPathFinder-Parallel_WithNameConversion_Rhoades.ipynb'. Once the functions are loaded, the BuildPaths function will allow you to search for an isotopologue at a desired path length (up to 15). You can do this by calling e.g. BuildPaths('Serine M+3',10) from the IPython interpreter. You will get a lot of warnings, thats ok, it should still run fine. Check the local directory to see that the csv result is written.

One could run the scripts purely out of a Jupyter Notebook as well, which we recommend if one wants to change the starting tracer or load their own model organism .dat from MetaCyc.
