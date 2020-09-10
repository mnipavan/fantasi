from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt

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
quad_coeff = 1.0

'''
Parameters for FEM algorithm
'''
q_degree = 3

'''
The velocity equation
'''
drdt=Expression("-G*(x[0]-mp)",G=quad_coeff, mp=midPoint, degree=3)

'''
Parameters for time-stepping
'''
dt=1e-15
num_steps=31
absTol=1e-11
relTol=1e-6

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
rho_next_CN=TrialFunction(V)
rho_next_alt=TrialFunction(V)
v0=TestFunction(V)
v1=TestFunction(V)

F_CN  = (rho_next_CN  - rho_curr)*v0*dx - 0.5*dt*velocity_n*rho_next_CN *v0.dx(0)*dx + 0.5*dt*D*rho_next_CN.dx(0) *v0.dx(0)*dx - 0.5*dt*velocity_p*rho_curr*v0.dx(0)*dx + 0.5*dt*D*rho_curr.dx(0)*v0.dx(0)*dx
F_alt = (rho_next_alt - rho_curr)*v1*dx - 0.6*dt*velocity_n*rho_next_alt*v1.dx(0)*dx + 0.6*dt*D*rho_next_alt.dx(0)*v1.dx(0)*dx - 0.4*dt*velocity_p*rho_curr*v1.dx(0)*dx + 0.4*dt*D*rho_curr.dx(0)*v1.dx(0)*dx
a_CN, L_CN  = lhs(F_CN), rhs(F_CN)
a_alt, L_alt = lhs(F_alt), rhs(F_alt)

#### Create VTK file for saving solution
vtkfile = File('GaussianTransient/solution.pvd')
bc=[]

#### Generate internal progess output to screen
progress = Progress('Time-stepping', num_steps)
set_log_level(LogLevel.PROGRESS)

print('Initial probability:')
print(assemble(rho_curr*dx))

#### Define the problem and solver
rho_next_CN  = Function(V)
rho_next_alt = Function(V)
probCN  = LinearVariationalProblem(a_CN , L_CN , rho_next_CN )
print("Done creating Crank-Nicolson problem...")
probAlt = LinearVariationalProblem(a_alt, L_alt, rho_next_alt)
print("Done creating lower order problem...")
solvCN  = LinearVariationalSolver(probCN)
print("Done creating solver for Crank-Nicolson problem...")
solvAlt = LinearVariationalSolver(probAlt)
print("Done creating solver for lower order problem...")

#### Initiate and start time stepper
t=0

for n in range(num_steps):
	t_next = t + dt
	# Perform step using both methods
	solvCN.solve()
	solvAlt.solve()

	print('Calculating next time from t = ' + str(t) + ' sec, using dt = ' + str(dt))

	# Calculate error
	tempErr = abs(np.array(rho_next_CN.vector()) - np.array(rho_next_alt.vector()))
	# Get the max relative error
	maxRelErr = abs(np.array(rho_next_CN.vector()))*relTol
	# Get the error tolerance value to compare errors with
	errTolArr = maxRelErr + absTol
	# Calculate ratio of error to tolerance
	errRatio = tempErr / errTolArr
	# Find the maximum ratio and adapt time-step
	maxErrRatio = max(errRatio)
	if max(errRatio) == 0.0:
		dt = 0.8*4.0*dt
	else:
		factor = 1.0/sqrt(maxErrRatio)
		if factor > 4.0:
			factor = 4.0
		dt = 0.8*dt*factor
	while maxErrRatio > 1.0:
		print('    using dt = ' + str(dt))
		# Update temporary time step (this will be iterated until error meets tolerance)
		t_next = t + dt

		# Solve using updated dt
		solvCN.solve()
		solvAlt.solve()

		# Calculate error
		tempErr = abs(np.array(rho_next_CN.vector()) - np.array(rho_next_alt.vector()))
		# Get the max relative error
		maxRelErr = abs(np.array(rho_next_CN.vector()))*relTol
		# Get the error tolerance value to compare errors with
		errTolArr = maxRelErr + absTol
		# Calculate ratio of error to tolerance
		errRatio = tempErr / errTolArr
		# Find the maximum ratio and adapt time-step
		maxErrRatio = max(errRatio)
		if max(errRatio) == 0.0:
			dt = 0.8*4.0*dt
		else:
			factor = 1.0/sqrt(maxErrRatio)
			if factor > 4.0:
				factor = 4.0
			dt = 0.8*dt*factor

	# Update current data point once error meets tolerance
	print('    ...accepted dt = ' + str(dt))
	rho_curr.assign(rho_next_CN)
	t = t_next
	if (n%10)==0:
		print('VTK File saved')
		vtkfile << (rho_curr, t)
		print('Updated probability:')
		print(assemble(rho_curr*dx))
	progress += 1
