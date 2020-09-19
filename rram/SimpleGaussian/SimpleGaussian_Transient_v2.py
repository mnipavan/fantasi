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
rho_D=Expression('1.0', degree = q_degree)
rho_curr=interpolate(rho_D,V)

#### Set up LLG equation that will enter variational form
velocity_n=interpolate(drdt,V_vec)
velocity_p=velocity_n

#### Set up variational form
#rho_next_CN=TrialFunction(V)
#rho_next_alt=TrialFunction(V)
#v0=TestFunction(V)
#v1=TestFunction(V)
rho_=TrialFunction(V)
v0=TestFunction(V)

#F_CN  = (rho_next_CN  - rho_curr)*v0*dx - 0.5*dt*velocity_n*rho_next_CN *v0.dx(0)*dx + 0.5*dt*D*rho_next_CN.dx(0) *v0.dx(0)*dx - 0.5*dt*velocity_p*rho_curr*v0.dx(0)*dx + 0.5*dt*D*rho_curr.dx(0)*v0.dx(0)*dx
#F_alt = (rho_next_alt - rho_curr)*v1*dx - 0.6*dt*velocity_n*rho_next_alt*v1.dx(0)*dx + 0.6*dt*D*rho_next_alt.dx(0)*v1.dx(0)*dx - 0.4*dt*velocity_p*rho_curr*v1.dx(0)*dx + 0.4*dt*D*rho_curr.dx(0)*v1.dx(0)*dx
#a_CN, L_CN  = lhs(F_CN), rhs(F_CN)
#a_alt, L_alt = lhs(F_alt), rhs(F_alt)
fpe_rhs  = velocity_n*rho_*v0.dx(0)*dx - D*rho_.dx(0)*v0.dx(0)*dx
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

