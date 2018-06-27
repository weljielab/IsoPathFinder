try:
    from IsoPathFinder import *
except:
    print('Did you build the reaction network already? If not, run ReactionModeling.py first')

#Example to run Serine m+3
BuildPaths('Serine M+3', 9)

