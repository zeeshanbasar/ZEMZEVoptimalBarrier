# ZEMZEVoptimalBarrier

working PoC of ZEM/ZEV-based, fuel optimal, collision avoidance and planetary landing algorithm

### UPDATE: 29-09-2022 ###

- testable code with multiple barrier shapes
- the constant l4 can be adjusted on case-by-case basis to avoid going through the barrier
- negative l4 may allow the trajectory to come back inside
- works against all possible cases


### UPDATE: 30-09-2022 ###
- it turns out the algo has good (except for a few corner cases) fuel optimality
- the odd bounce off issue doesnt happend anymore, idk how; dont ask, dont care
- uploading the working code, plus a few test cases and test case generation file


### UPDATE: 13-10-2022 ###
- fixed the issue where trajectory was bouncing off the higher level, instead of allowable lower level
-- replaced 2-norm with inf-norm
- moved the barrier 5 down by 500m; did it to make the barrier similar to the one in Gong et al. (2022)
- so for smaller tf the trjaectories are much, much smoother, but the acceleration demand is obviously very high
