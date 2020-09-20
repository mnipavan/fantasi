from dolfin import *
from mshr import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''
Plot control
'''
outPlotInteractive=0

'''
Mesh
'''
global_normal = Expression(("x[0]", "x[1]", "x[2]"), degree=3)
meshid = "sphere_ico6"
meshfile = "meshes/%s.xml.gz" % meshid
mesh=Mesh(meshfile)
mesh.init_cell_orientations(global_normal)
OutLine = "Mesh dimensions: %d" % mesh.topology().dim()
print(OutLine)

'''
Time stepping and time integration
'''
tscale=1e-9                      # Scale factor to convert simulation time to real-time
num_steps = 31                   # total number of steps
#dt = T/float(num_steps)         # time step size
dt = 1e-18 / tscale              # time step size (in simulation time)
dtmax = 1e-2			 # maximum time step
tfactor = 2.0                    # maximum scale factor to expand dt
tfinal = 10e-9 / tscale          # end-time of simulation
#dt = 1e-9                        # time step size (in ns)
q_degree = 3
#dx=dx(metadata={'quadrature_degree': q_degree})
absTol=1e-15
relTol=1e-7

'''
Constant parameters for LLG problem
'''
kBoltzmann = 1.38064852e-23    # in J/K
mu0 = 4*np.pi * 1.0e-7         # in N/A^2
gamFac = 1.7595e11              # in rad/(s.T)
#gamFac = 1.7595e2              # in rad/(ns.T)

'''
Free layer description
'''
alpha = 0.0135                                   # Unitless damping factor
Ms =450e3                                        # in A/m
t_FL=1.0e-9                                      # thickness of FL
diameter=40.0e-9                                 # diameter of FL with circular cross-section
#length=50e-9
#width=3*length
magVolume = t_FL * (diameter**2) * (np.pi/4.0)   # in m^3
#magVolume=(np.pi/4)*length*width*thickness
G=Constant((gamFac*tscale*mu0)/(1+alpha**2))            # Scale factor for LLG equation

'''
Temperature parameters
'''
Temperature = 300                                                                         # in K
D = Constant(alpha * gamFac * tscale * kBoltzmann * Temperature / ((1+alpha**2)*Ms * magVolume))   # per unit time based on that for gamFac

'''
Definitions for uniaxial anisotropy
'''
delta=44.0
Eb=delta*kBoltzmann*Temperature
Ku2=Eb/magVolume
H_uni=Constant((2*Ku2)/(mu0*Ms))

'''
The LLG equation
'''
dmdt=Expression(("-G*(a*x[0]*x[2]+x[1])*H*x[2]","G*(x[0]-a*x[1]*x[2])*H*x[2]","a*G*(1-x[2]*x[2])*H*x[2]"),G=G, a=alpha, H=H_uni, degree=1)

'''
Set up variational form of Fokker-Planck equation...
'''

#### Basis space
#V=FunctionSpace(mesh,'P',3)
#V_vec=VectorFunctionSpace(mesh,'P',degree=3, dim=3)
V=FunctionSpace(mesh,'CG',q_degree)
V_vec=VectorFunctionSpace(mesh,'CG',degree=q_degree,dim=3)

#### define initial value
rho_D=Expression('1/(4*pi)', degree = 1)
rho_curr=interpolate(rho_D,V)

#### Set up LLG equation that will enter variational form
velocity_n=interpolate(dmdt,V_vec)
velocity_p=velocity_n

#### Set up variational form
rho_next_CN=TrialFunction(V)
rho_next_alt=TrialFunction(V)
v0=TestFunction(V)
v1=TestFunction(V)

F_CN  = (rho_next_CN  - rho_curr)*v0*dx - 0.5*dt*dot(velocity_n*rho_next_CN , grad(v0))*dx + 0.5*dt*D*dot(grad(rho_next_CN ),grad(v0))*dx - 0.5*dt*dot(velocity_p*rho_curr,grad(v0))*dx + 0.5*dt*D*dot(grad(rho_curr),grad(v0))*dx
F_alt = (rho_next_alt - rho_curr)*v1*dx - 0.6*dt*dot(velocity_n*rho_next_alt, grad(v1))*dx + 0.6*dt*D*dot(grad(rho_next_alt),grad(v1))*dx - 0.4*dt*dot(velocity_p*rho_curr,grad(v1))*dx + 0.4*dt*D*dot(grad(rho_curr),grad(v1))*dx
a_CN, L_CN  = lhs(F_CN), rhs(F_CN)
a_alt, L_alt = lhs(F_alt), rhs(F_alt)

#### Create VTK file for saving solution
vtkfile = File('result_files/solution.pvd')
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
n=0

if (outPlotInteractive != 0):
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
	x = np.cos(u)*np.sin(v)
	y = np.sin(u)*np.sin(v)
	z = np.cos(v)
	rho_var = rho_curr(Point(x,y,z))
	ax.plot_surface(x,y,z, color='g', edgecolor='none')
	plt.show()

while (t < tfinal):
	n = n+1

	if (dt > dtmax):
		dt = dtmax

	t_next = t + dt

	if (t_next > tfinal):
		dt = tfinal - t
		t_next = tfinal

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
		dt = 0.8*tfactor*dt
	else:
		factor = 1.0/sqrt(maxErrRatio)
		if factor > tfactor:
			factor = tfactor
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
			dt = 0.8*tfactor*dt
		else:
			factor = 1.0/sqrt(maxErrRatio)
			if factor > tfactor:
				factor = tfactor
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
		if (outPlotInteractive != 0):
			plt.figure()
			plot(rho_curr)
			plt.show()

	progress += 1
