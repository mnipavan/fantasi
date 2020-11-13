The directory contains examples to demonstrate the use of FANTASI in simulating the
relaxation of monodomain ferromagnets to thermal equilibrium in the presence of a
Zeeman field.

large_field_01.py: initial magnetization is randomly distributed. The monodomain has
                   uniaxial anisotropy in the z-axis (equivalent to 44kT energy 
				   barrier @ T = 300K) while the applied magnetic field along the
				   +y-direction is 5 times its coercive field. The simulation shows
				   the evolution (over 10 seconds) of the magnetization from uniform
				   distribution on the unit sphere to a single modes in the
				   +y-direction.

large_field_02.py: same as large_field_01 except the applied field is along the
                   +z-direction.

large_field_03.py: same as large_field_02 except for 'P' finite elements instead of
                   Lagrange finite elements.

large_field_04.py: same as large_field_02 except the number of degrees of freedom is
                   2 instead of 3.
