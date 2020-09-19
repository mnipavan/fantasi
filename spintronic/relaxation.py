from dolfin import *
from mshr import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from gryphon import *

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
#dt = 1e-9                        # time step size (in ns)
q_degree = 3
#dx=dx(metadata={'quadrature_degree': q_degree})
absTol=1e-15
relTol=1e-5

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
D = Constant(alpha * gamFac * kBoltzmann * Temperature / ((1+alpha**2)*Ms * magVolume))   # per unit time based on that for gamFac

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
#rho_next_N=TrialFunction(V)
#rho_next_alt=TrialFunction(V)
#v0=TestFunction(V)
#v1=TestFunction(V)
rho_=TrialFunction(V)
v0=TestFunction(V)

#F_CN  = (rho_next_CN  - rho_curr)*v0*dx - 0.5*dt*dot(velocity_n*rho_next_CN , grad(v0))*dx + 0.5*dt*D*dot(grad(rho_next_CN ),grad(v0))*dx - 0.5*dt*dot(velocity_p*rho_curr,grad(v0))*dx + 0.5*dt*D*dot(grad(rho_curr),grad(v0))*dx
#F_alt = (rho_next_alt - rho_curr)*v1*dx - 0.6*dt*dot(velocity_n*rho_next_alt, grad(v1))*dx + 0.6*dt*D*dot(grad(rho_next_alt),grad(v1))*dx - 0.4*dt*dot(velocity_p*rho_curr,grad(v1))*dx + 0.4*dt*D*dot(grad(rho_curr),grad(v1))*dx
#a_CN, L_CN  = lhs(F_CN), rhs(F_CN)
#a_alt, L_alt = lhs(F_alt), rhs(F_alt)
fpe_rhs  = dot(velocity_n*rho_, grad(v0))*dx - D*dot(grad(rho_),grad(v0))*dx
T = [0, 10]

#### Create VTK file for saving solution
vtkfile = File('result_files/solution.pvd')
bc=[]

print('Initial probability:')
print(assemble(rho_curr*dx))

obj = ESDIRK(T, rho_curr, fpe_rhs, bcs=[], tdfBC=[], tdf=[])

# Set up time-stepping control
obj.parameters["timestepping"]["dtmin"] = 1e-18
obj.parameters["timestepping"]["dtmax"] = 1e-2
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
