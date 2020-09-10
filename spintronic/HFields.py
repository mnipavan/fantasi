from fenics import *
from mshr import *
import numpy as np

'''
The LLG equation for Zeeman fields defined by (Hx, Hy and Hz)

The functions expect gam_fac to be "gamma / (1+alpha*alpha)"
'''
def dmdt_happx(gam_fac, alph_damp, Hx, q_degree):
	dmdt=Expression(("G*a*(1-x[0]*x[0])*H","-1.0*G*(m[2]+a*m[0]*m[1])*H","G*(m[1]-a*m[0]*m[2])*H"),G=gam_fac, a=alph_damp, H=Hx, degree=q_degree)
	return dmdt

def dmdt_happy(gam_fac, alph_damp, Hy, q_degree):
	dmdt=Expression(("G*(m[2]-a*m[0]*m[1])*H","G*a*(1-m[1]*m[1])*H","-1.0*G*(m[0]+a*m[1]*m[2])*H"),G=gam_fac, a=alph_damp, H=Hy, degree=q_degree)
	return dmdt

def dmdt_happz(gam_fac, alph_damp, Hz, q_degree):
	dmdt=Expression(("-1.0*G*(m[1]+a*m[0]*m[2])*H","G*(m[0]-a*m[1]*m[2])*H","G*a*(1-m[2]*m[2])*H"),G=gam_fac, a=alph_damp, H=Hz, degree=q_degree)
	return dmdt

'''
The LLG equation for uniaxial anisotropy fields defined by (Huax, Huay and Huaz)

The functions expect gam_fac to be "gamma / (1+alpha*alpha)"
These functions can be combined to calculate the shape anisotropy
effective field
'''
def dmdt_huax(gam_fac, alph_damp, Huax, q_degree):
	dmdt=Expression(("G*a*(1-x[0]*x[0])*H*m[0]","-1.0*G*(m[2]+a*m[0]*m[1])*H*m[0]","G*(m[1]-a*m[0]*m[2])*H*m[0]"),G=gam_fac, a=alph_damp, H=Huax, degree=q_degree)
	return dmdt

def dmdt_huapy(gam_fac, alph_damp, Huay, q_degree):
	dmdt=Expression(("G*(m[2]-a*m[0]*m[1])*H*m[1]","G*a*(1-m[1]*m[1])*H*m[1]","-1.0*G*(m[0]+a*m[1]*m[2])*H*m[1]"),G=gam_fac, a=alph_damp, H=Huay, degree=q_degree)
	return dmdt

def dmdt_huaz(gam_fac, alph_damp, Huaz, q_degree):
	dmdt=Expression(("-1.0*G*(m[1]+a*m[0]*m[2])*H*m[2]","G*(m[0]-a*m[1]*m[2])*H*m[2]","G*a*(1-m[2]*m[2])*H*m[2]"),G=gam_fac, a=alph_damp, H=Huaz, degree=q_degree)
	return dmdt

