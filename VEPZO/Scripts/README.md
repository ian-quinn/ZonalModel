Sample code for model serialization, set.py as main entry.

mesh.py tesselates one not self-intersected polygon into rough mesh and assign solar energy to it. The input should be a nested list of 2D point sequence starting with the exterior boundary that followed by other vertice loops representing the inner holes.

```python
# input
loops = [[[3, 0], [10, 0], [10, 7], [7, 10], [0, 10], [0, 3], [3, 0]], 
		 [[3, 3], [7, 3], [7, 7], [3, 7], [3, 3]]]
mask_wwr = [0.8, 0.8, 0.2, 0.2, 0.2, 0.8]
mask_adia = [0, 0, 0, 1, 1, 0]
# output mask for meshing with scale factor 1 and 0.5
[[0 3 3 3 3 3 3 3 0 0 0 0]
 [3 1 1 1 1 1 1 1 0 0 0 0]
 [3 1 1 1 1 1 1 1 1 0 0 0]
 [3 1 1 1 1 1 1 1 1 1 0 0]				[[0 3 3 3 3 0 0]
 [3 1 1 1 1 1 1 1 1 1 1 0]				 [3 1 1 1 1 0 0]
 [3 1 1 1 1 2 2 1 1 1 1 0]				 [3 1 1 1 1 1 0]
 [3 1 1 1 1 2 2 1 1 1 1 0]				 [3 1 1 2 1 1 0]
 [3 1 1 1 1 1 1 1 1 1 1 0]				 [3 1 1 1 1 1 0]
 [0 0 1 1 1 1 1 1 1 1 1 0]				 [0 0 1 1 1 1 0]
 [0 0 0 1 1 1 1 1 1 1 1 0]				 [0 0 0 0 0 0 0]]
 [0 0 0 0 1 1 1 1 1 1 1 0]
 [0 0 0 0 0 0 0 0 0 0 0 0]]
# also visualize solar load distribution in snapshots
# see gif animation below (matplotlib)
```
<img src="animation.gif?raw=true">

0 - outside
1 - floorplan
2 - interior void
3 - outside (adiabatic)

set.py translates the mask into modelica model with proper boundary condition assigned.