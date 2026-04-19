# TODO
- give pre setting function to the gui
- neighbor list memory structure improvement
- remake the gpu structure
- superquadric (maybe based on liggghts?)

# Current glitches
- nondimensionalization seem to be glitched. turned off for now
- gpu doesn't match with cpu results, but possibly caused by the difference in round off errors
- sort cell list doesn't match with naive implementation

# WIP
- morton key sorted BVH in gpu was implemented but needs check.
- maybe, need fix in BVH. doesn't match with cell list when particle collides each other. Matches when particles doesn't collide to each other.

