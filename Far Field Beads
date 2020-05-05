import numpy as np
import fmtrack
import math

cell_threshold = 1.0
X_DIM = 149.95
Y_DIM = 149.95
Z_DIM = 140.0

filenames_final_cell = './06092019 Info/Cell_CytoD/Gel 1 CytoD%s.tif'
filenames_init_cell = './06092019 Info/Cell_Normal/Gel 1 Normal2%s.tif'
filenames_final_beads = './06092019 Info/Beads_CytoD/Gel 1 CytoD%s.tif'
filenames_init_beads = './06092019 Info/Beads_Normal/Gel 1 Normal2%s.tif'

mesh_init = fmtrack.FMMesh()
mesh_init.get_cell_surface(filenames_init_cell, 0, X_DIM, Y_DIM, Z_DIM, cell_threshold)

mesh_final = fmtrack.FMMesh()
mesh_final.get_cell_surface(filenames_final_cell, 0, X_DIM, Y_DIM, Z_DIM, cell_threshold)

beads_init = fmtrack.FMBeads()
beads_init.get_bead_centers(filenames_init_beads, 1, X_DIM, Y_DIM, Z_DIM)

beads_final = fmtrack.FMBeads()
beads_final.get_bead_centers(filenames_final_beads, 1, X_DIM, Y_DIM, Z_DIM)


points_init = mesh_init.points
points_final = mesh_final.points 

# Find shortest distance  between each initial bead and the cell boundary
kk = 0
X,Y,Z = beads_init.get_xyz()
c1,c2,c3 = points_init.get_xyz()
distances = np.zeros(len(X))
temp = np.zeros(len(c1))
for bead in beads_init:
    for point in points_init: 
        temp[point] = math.sqrt((X[bead]-c1[point])**2 + (Y[bead]-c2[point])**2 + (Z[bead]-c3[point])**2)
    
    min_dist = np.amin(temp)
    if min_dist > 50:
        distances[kk] = [X[bead], Y[bead], Z[bead]]
        kk = kk + 1

# Trim off the end zeros that were never filled
np.trim_zeros(distances)
    
# the d array should be all the far field beads

