from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
from gryphon import *

'''
Interactive Output Setting
'''
outInteractive = 0

'''
Mesh
'''
nx=10001
mesh=UnitIntervalMesh(nx)
OutLine = "Mesh dimensions: %d" % mesh.topology().dim()
print(OutLine)

'''
Constant parameters for relaxation problem
'''
midPoint = 0.5

'''
Temperature parameters
'''
alpha = 0.1
Temperature = 300
D = Constant(alpha * Temperature)   # per unit time based on that for gamFac

'''
Definitions for quadratic well
'''
quad_coeff = 100.0

'''
Parameters for FEM algorithm
'''
q_degree = 3

'''
The velocity equation
'''
drdt=Expression("-G*(x[0]-mp)",G=quad_coeff, mp=midPoint, degree=3)

'''
Set up variational form of Fokker-Planck equation...
'''

#### Basis space
V=FunctionSpace(mesh, 'CG', q_degree)
V_vec=FunctionSpace(mesh, 'CG', q_degree)

#### define initial values
rho_D=Expression('6.0 * x[0] * (1.0 - x[0])', degree = q_degree)
rho_curr=interpolate(rho_D,V)

#### Set up drift equation that will enter variational form
velocity_n=interpolate(drdt,V_vec)
velocity_p=velocity_n

#### Set up variational form for RHS of FPE ####
rho_=TrialFunction(V)
v0=TestFunction(V)
fpe_rhs  = velocity_n*rho_*v0.dx(0)*dx - D*rho_.dx(0)*v0.dx(0)*dx

'''
Time interval to computer result
'''
T = [0, 100]

print('Initial probability:')
print(assemble(rho_curr*dx))

#### Initiate and start time stepper
if (outInteractive != 0):
	plt.figure()
	plot(rho_curr.root_node(), title="Initial solution")
	plt.xlabel("Position")
	plt.ylabel("PDF")
	plt.show()

obj = ESDIRK(T, rho_curr, fpe_rhs, bcs=[], tdfBC=[], tdf=[])

# Set up time-stepping control
obj.parameters["timestepping"]["dtmin"] = 1e-18
obj.parameters["timestepping"]["dtmax"] = 1e-1
obj.parameters["timestepping"]["stepsizeselector"] = "gustafsson"
obj.parameters["timestepping"]["convergence_criterion"] = "relative"

# Set up solver verbosity
obj.parameters["verbose"] = False

# Save plot of each time step in VTK format.
obj.parameters["output"]["plot"] = True

# Set that the plot of selected step sizes should be saved in jpg.
# Available choices are jpg, png and eps.
obj.parameters["output"]["imgformat"] = "jpg"

# Call the solver which will do the actual calculation.
obj.solve()

