# TODO
- add json setting file  
- give pre setting function to the gui
- neighbor list memory structure improvement
- superquadric (maybe based on liggghts?)

# Current glitches
- gpu doesn't match with cpu results, but possibly caused by the difference in round off errors
- sort cell list doesn't match with naive implementation
- gpu calculation result changes everytime it runs at large number of particles

# WIP
- morton key sorted BVH in gpu
- maybe, need fix in BVH. doesn't match with cell list when particle collides each other. Matches when particles doesn't collide to each other.

