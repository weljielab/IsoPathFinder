# IsoPathFinder
Workflow for isotope mapping


# You can run the program one of two ways:
The first is from the terminal. Open IPython by typing 'ipython' in the terminal. You can then load the IsoPathFinder functions by typing '%run IsoPathFinder-Parallel_WithNameConversion_Rhoades.ipynb'. Once the functions are loaded, the BuildPaths function will allow you to search for an isotopologue at a desired path length (up to 15). You can do this by calling e.g. BuildPaths('Serine M+3',10) from the IPython interpreter. You will get a lot of warnings, thats ok, it should still run fine. Check the local directory to see that the csv result is written.

The other way is to run the script out of a Jupyter Notebook. Open jupyter by calling 'jupyter notebook' from the terminal, launching a Firefox page. From here you can open the IsoPathFinder-Parallel_WithNameConversoin_Rhoades.ipynb notebook. Run all the cells to load the functions, then at the very bottom you can add a cell and call e.g. BuildPaths('Serine M+3',10). The advantage here is that you can see the code and modify if problems come up. You will get a lot of warnings, thats ok, it should still run fine. Check the local directory to see that the csv result is written, which you can copy back to the lab shared drive via 'sudo cp Serine\ M+3_Paths_10Rxns.csv /media/sf_pythonVM/'
