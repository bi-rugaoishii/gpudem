# Current glitches
- gpu doesn't match with cpu results, but possibly caused by the difference in round off errors
- sort cell list doesn't match with naive implementation

# WIP
- morton key sorted BVH in gpu
- maybe, need fix in BVH. doesn't match with cell list when particle collides each other. Matches when particles doesn't collide to each other.
# TODO
- fix the convention of normal displacement sign and tangential force for gpu
- gpu verlet
- neighbor list memory structure improvement

