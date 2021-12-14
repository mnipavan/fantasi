from dolfin import *
from mshr import *
import numpy as np

global_normal = Expression(("x[0]", "x[1]", "x[2]"), degree=3)
meshid = "comsol_extreme_fine_sphere"
meshfile = "%s.xml" % meshid
mesh = Mesh(meshfile)
mesh.init_cell_orientations(global_normal)
OutLine = "Mesh dimensions: %d" % mesh.topology().dim()
print(OutLine)

simName = "check_mesh"

q_degree = 1
V = FunctionSpace(mesh,'CG',q_degree)

rho_D = Expression('1/(4*pi)', degree=q_degree)
x_Val = Expression('x[0]', degree=q_degree)
y_Val = Expression('x[1]', degree=q_degree)
z_Val = Expression('x[2]', degree=q_degree)
x_Arr = interpolate(x_Val,V)
y_Arr = interpolate(y_Val,V)
z_Arr = interpolate(z_Val,V)
rho_curr = interpolate(rho_D,V)

outDirName = simName+"_results"
vtkfile = File(outDirName+"/solution.pvd")
print('VTK File saved')
vtkfile << (rho_curr, 0)

print('Initial probability:')
print(assemble(rho_curr*dx))
print('Average X:')
print(assemble(x_Arr*rho_curr*dx))
print('Average Y:')
print(assemble(y_Arr*rho_curr*dx))
print('Average Z:')
print(assemble(z_Arr*rho_curr*dx))
