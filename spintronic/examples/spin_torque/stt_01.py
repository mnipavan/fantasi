from dolfin import *
from mshr import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from gryphon import *
from HFields import *

'''
FEniCS controls
'''
parameters["linear_algebra_backend"] = "PETSc"
parameters["mesh_partitioner"] = "ParMETIS"

'''
Plot control
'''
outputStats = False

'''
FANTASI simulation name
'''
simName = "stt_01"

'''
Mesh
'''
global_normal = Expression(("x[0]", "x[1]", "x[2]"), degree=3)
meshid = "sphere_ico6"
meshfile = "meshes/%s.xml.gz" % meshid
mesh = Mesh(meshfile)
mesh.init_cell_orientations(global_normal)
OutLine = "Mesh dimensions: %d" % mesh.topology().dim()
print(OutLine)

'''
Control parameters for FEM solver
'''
q_degree = 3
absTol = 1e-15
relTol = 1e-5

'''
Control parameters for time stepping and time integration
'''
tscale = 1e-9                    # Scale factor to convert simulation time to real-time
num_steps = 31                   # total number of steps
dt = 1e-18 / tscale              # time step size (in simulation time)
startDt = 1e-4                   # First time step to try (time-scaled)
startMaxDt = 1e-2                # Maximum Dt in time-stepper (time-scaled)
startMinDt = 1e-18               # Minimum Dt in time-stepper (time-scaled)
T = [0, 1]                       # Time range of simulation per stage (time-scaled)
stageCountOuter = 1              # Number of stages to simulate
stageCountInner = 20             # Number of stages to simulate

'''
Constant parameters for LLG problem
'''
kBoltzmann = 1.38064852e-23    # in J/K
mu0 = 4*np.pi * 1.0e-7         # in N/A^2
gamFac = 1.7595e11             # in rad/(s.T)

'''
Free layer description
'''
alpha = 0.0135                                   # Unitless damping factor
Ms = 450e3                                       # in A/m
t_FL = 1.0e-9                                    # thickness of FL
diameter = 40.0e-9                               # diameter of FL with circular cross-section
magVolume = t_FL * (diameter**2) * (np.pi/4.0)   # in m^3
G = Constant((gamFac*tscale*mu0)/(1+alpha**2))   # Scale factor for LLG equation (time-scaled)

'''
Temperature parameters
'''
Temperature = 300                                                                                  # in K
D = Constant(alpha * gamFac * tscale * kBoltzmann * Temperature / ((1+alpha**2)*Ms * magVolume))   # per unit time based on that for gamFac

'''
Definitions for uniaxial anisotropy
'''
delta = 44.0
Eb = delta*kBoltzmann*Temperature
Ku2 = Eb/magVolume
H_uni = Constant((2*Ku2)/(mu0*Ms))

'''
Definitions for spin torque
'''
P_fix = 0.5
P_free = 0.5
Lambda_fix = 1.0
Lambda_free = 1.0
epsPrime = 0.0
Icurr = 25e-6
mp = np.array([0.0, 0.0, -1.0])

'''
The LLG equation
'''
dmdt = dmdt_huaz(gam_fac=G, alph_damp=alpha, Huaz=H_uni, q_degree=q_degree)
dmdt_stt = dmdt_mp(gam_fac=G, alph_damp=alpha, Pfix=P_fix, Pfree=P_free, LambFree=Lambda_free, LambFix=Lambda_fix, Ms_=Ms, Icurr=Icurr, vol=magVolume, epsPrime=epsPrime, mp=mp, q_degree=q_degree)

'''
Set up variational form of Fokker-Planck equation for initial value problem (IVP)
'''

#### Basis space
V = FunctionSpace(mesh,'CG',q_degree)
V_vec = VectorFunctionSpace(mesh,'CG',degree=q_degree,dim=3)

#### Create time series file to import nodal values
timeseries_rho = TimeSeries('rho_series')

#### Define initial value on the mesh
rho_curr = Function(V)
timeseries_rho.retrieve(rho_curr.vector(),0.0)

#### Set up LLG equation to be solved
velocity_uniaxial = interpolate(dmdt,V_vec)
velocity_stt = interpolate(dmdt_stt,V_vec)

#### Set up variational form
rho_ = TrialFunction(V)
v0 = TestFunction(V)
fpe_rhs = dot((velocity_uniaxial+velocity_stt)*rho_, grad(v0))*dx - D*dot(grad(rho_),grad(v0))*dx

#### Create VTK file for saving solution and save initial value
outDirName = simName+"_results"
vtkfile = File(outDirName+"/solution.pvd")
print('VTK File saved')
vtkfile << (rho_curr, 0)

#### Create time series file to save nodal values
timeseries_rho = TimeSeries(outDirName+"/rho_series")

#### Perform initial integration to get estimated error in the beginning
print('Initial probability:')
print(assemble(rho_curr*dx))

'''
Using Gryphon toolbox to perform time-stepping
'''
#### Start the solver to calculate transient solution
for idx1 in range(0, stageCountOuter):
    for idx in range(0, stageCountInner):
        #### Get initial Gryphon object
        obj = ESDIRK(T, rho_curr, fpe_rhs, bcs=[], tdfBC=[], tdf=[], method="mumps")

        #### Set up Gryphon time-stepping control
        obj.parameters["timestepping"]["dtmin"] = startMinDt
        obj.parameters["timestepping"]["dtmax"] = startMaxDt
        obj.parameters["timestepping"]["stepsizeselector"] = "gustafsson"      # Time step adaptation scheme
        obj.parameters["timestepping"]["convergence_criterion"] = "relative"   # Error check scheme
        obj.parameters["timestepping"]["dt"] = startDt                         # First time step to try

        #### Set up solver verbosity
        obj.parameters["verbose"] = True

        #### Save plot of solution at every internal time step
        obj.parameters["output"]["plot"] = False

        #### Set that the plot of selected step sizes should be saved in jpg.
        #### Available choices are jpg, png and eps.
        if outputStats:
            obj.parameters["output"]["imgformat"] = "jpg"

        obj.solve()

        rho_curr = obj.u
        t = obj.t
        print('VTK File saved')
        vtkfile << (rho_curr, (t*(1+idx) + idx1*T[1]))
        print('Updated probability:')
        print(assemble(rho_curr*dx))

    startDt = startDt * 10
    startMaxDt = startMaxDt * 10
    startMinDt = startMinDt * 10

    T[1] = 10 * T[1]

#### Save final solution in HDF5 file for easy import by other codes
timeseries_rho.store(rho_curr.vector(), 0.0)
