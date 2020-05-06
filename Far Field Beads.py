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


points = mesh_init.points
p1,p2,p3 = (points[:,0],points[:,1],points[:,2]) 

beads = beads_init.points
b1,b2,b3 = (beads[:,0],beads[:,1],beads[:,2])
# Find shortest distance  between each initial bead and all the cell mesh vertices
kk = 0
new_beads = np.zeros((len(beads), 3))
temp = np.zeros(len(points))
for bead in range(0,len(beads)):
    jj = 0
    for point in range(0,len(points)): 
        temp[jj] = math.sqrt((b1[bead]-p1[point])**2 + (b2[bead]-p2[point])**2 + (b3[bead]-p3[point])**2)
        jj = jj + 1
    
    min_dist = np.amin(temp)
    print(min_dist)
    if min_dist > 50:
        new_beads[kk] = [b1[bead],b2[bead],b3[bead]]
        kk = kk + 1


# Trim off the end zeros that were never filled
new_beads = new_beads[:kk, :]
b1,b2,b3 = (new_beads[:,0],new_beads[:,1],new_beads[:,2])   
Far_beads = fmtrack.FMBeads()
Far_beads.load_from_positions(b1, b2, b3)
Far_beads.save_as_txt('./06092019 Info/Far Field Beads')

