#!/bin/bash

#Build the labeled metabolic network - note this example builds 20 reaction rounds
python3 ReactionModeling.py

#Check that the network exists, then run a BuildPaths('Serine M+3', 9) example
python3 IsoPathFinderDemo.py
