The directory contains examples to demonstrate the use of FANTASI in simulating the
relaxation of monodomain ferromagnets to thermal equilibrium in the presence of a
spin current injection (spin-transfer torque).

stt_01.py: initial magnetization is read in from a user-defined HDF5 file. The
                   monodomain has uniaxial anisotropy in the z-axis (equivalent to
				   44kT energy barrier @ T = 300K) while the spin current injected
				   is polarized in the +z-direction. The simulation shows the
				   evolution (over 20 nanoseconds) of the magnetization when the
				   charge current flowing is 25 uA.
