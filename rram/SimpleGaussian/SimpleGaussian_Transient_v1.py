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
alpha = 0.01
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

#### Create VTK file for saving solution
vtkfile = File('GaussianTransient/solution.pvd')

'''
Set up variational form of Fokker-Planck equation...
'''

#### Basis space
V=FunctionSpace(mesh, 'CG', q_degree)
V_vec=FunctionSpace(mesh, 'CG', q_degree)

#### define initial values
sgma_val = 0.05
rho_D=Expression('aFac * exp(-0.5 * (x[0] - 0.5) * (x[0] - 0.5) / (sig * sig))', aFac = 1.0 / (sgma_val*np.sqrt(2.0*np.pi)), sig = sgma_val, dnfac = 2.0*np.pi,degree = q_degree)
rho_curr=interpolate(rho_D,V)
print('Initial probability (pre-adjustment):')
err_init=1.0 - assemble(rho_curr*dx)
print(assemble(rho_curr*dx))
rho_D=Expression('shift + aFac * exp(-0.5 * (x[0] - 0.5) * (x[0] - 0.5) / (sig * sig))', shift=err_init, aFac = 1.0 / (sgma_val*np.sqrt(2.0*np.pi)), sig = sgma_val, dnfac = 2.0*np.pi,degree = q_degree)
rho_curr=interpolate(rho_D,V)
print('Initial probability (post-adjustment):')
print(assemble(rho_curr*dx))

#### Set up LLG equation that will enter variational form
velocity_n=interpolate(drdt,V_vec)
velocity_p=velocity_n

#### Set up variational form
rho_=TrialFunction(V)
v0=TestFunction(V)
fpe_rhs  = velocity_n*rho_*v0.dx(0)*dx - D*rho_.dx(0)*v0.dx(0)*dx

#### Initiate and start time stepper
if (outInteractive != 0):
	plt.figure()
	plot(rho_curr.root_node(), title="Initial solution")
	plt.xlabel("Position")
	plt.ylabel("PDF")
	plt.show()

# Set up time-stepping control
T = [0, 1]
obj = ESDIRK(T, rho_curr, fpe_rhs, bcs=[], tdfBC=[], tdf=[], method="mumps")
obj.parameters["timestepping"]["dtmin"] = 1e-18
obj.parameters["timestepping"]["dtmax"] = 1e-1
obj.parameters["timestepping"]["dt"] = 1e-5
obj.parameters["timestepping"]["stepsizeselector"] = "gustafsson"
obj.parameters["timestepping"]["convergence_criterion"] = "relative"

# Set up solver verbosity
obj.parameters["verbose"] = True

# Save plot of each time step in VTK format.
obj.parameters["output"]["plot"] = False

# Set that the plot of selected step sizes should be saved in jpg.
# Available choices are jpg, png and eps.
#obj.parameters["output"]["imgformat"] = "jpg"
n = 1;

while True:
	# Call the solver which will do the actual calculation.
	obj.solve()
	rho_curr = obj.u
	print('VTK File saved')
	vtkfile << (rho_curr, n*T[1])
	print('Updated probability:')
	print(assemble(rho_curr*dx))
	n += 1
	if n > 2:
		break
	obj.t = 0.0
	obj.dt = 1e-5

if (outInteractive != 0):
        plt.figure()
        plot(rho_curr.root_node(), title="Final solution")
        plt.xlabel("Position")
        plt.ylabel("PDF")
        plt.show()

