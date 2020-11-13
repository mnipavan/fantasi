The directory contains examples to demonstrate the use of FANTASI in simulating the
relaxation of monodomain ferromagnets to thermal equilibrium.

relaxation_01.py: initial magnetization is randomly distributed. The monodomain has
                  uniaxial anisotropy in the z-axis (equivalent to 44kT energy 
				  barrier @ T = 300K). The simulation shows the evolution of the
				  magnetization from uniform distribution on the unit sphere to
				  to modes along the z-axis.
				  
relaxation_02.py: Same as relaxation_01.py except initial condition is imported
                  from a HDF5 file.

relaxation_03.py: Same as relaxation_02.py except number of degrees of freedom is 2
                  instead of 3.

relaxation_04.py: Same as relaxation_02.py except uniaxial anisotropy is along
                  y-direction.

relaxation_05.py: Same as relaxation_02.py except for P finite elements instead for
                  Lagrange finite elements.
