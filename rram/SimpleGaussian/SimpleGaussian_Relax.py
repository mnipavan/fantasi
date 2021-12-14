from fenics import *

'''
Mesh
'''
nx=1001
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
quad_coeff = 20.0

'''
Parameters for FEM algorithm
'''
q_degree = 1

'''
Set up variational form of Fokker-Planck equation...
'''

'''
The velocity equation
'''
drdt=Expression("-G*(x[0]-mp)", G=quad_coeff, mp=midPoint, degree=q_degree)

#### Basis space
V=FunctionSpace(mesh, 'CG', q_degree)
V_vec=FunctionSpace(mesh, 'CG', q_degree)

#### define initial value
rho_D=Expression('1.0', degree = q_degree)
rho_curr=interpolate(rho_D,V)

#### velocity function
velocity=interpolate(drdt,V_vec)

#### Set up variational form
rho_ss=TrialFunction(V)
v0=TestFunction(V)

F_CN = -1.0*velocity*rho_ss*v0.dx(0)*dx + D*rho_ss.dx(0)*v0.dx(0)*dx
a, L = lhs(F_CN), rhs(F_CN)

#### Create VTK file for saving solution
vtkfile = File('GaussianRelax/solution.pvd')

print('Initial probability:')
print(assemble(rho_curr*dx))

#### Define the problem and solve
rho_ss  = Function(V)
solve(a == L, rho_ss)

#### Save solution to VTK file
vtkfile << rho_ss

plot(rho_ss)
plot(mesh)

print('Final probability:')
print(assemble(rho_ss*dx))
