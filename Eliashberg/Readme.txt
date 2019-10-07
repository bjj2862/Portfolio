- many functions(ex, self-energy, Green function, etc.)are struct. see struct.h

Eliashber.c:
- main program
1. put initial values and calculation initial Green function
- see para.h, initialize.h and Cal_Green_function in cal_function.h
2. iteration loop
- see iter.h
2.1 the function iter-eff-self in iter.h
 - Fourier transformation the spectral function A from momentum to real (see some_function.h )
 - see cal_eff_self_energy in cal_function.h
 - inverse Fourier transformation the self-energy from real to momentum (see some_function.h )
2.2  the function iter-imp-self in iter.h
 - Fourier transformation the spectral function G from momentum to real (see some_function.h )
 - see cal_imp_self_energy in cal_function.h
 - inverse Fourier transformation the self-energy from real to momentum (see some_function.h )

3. print output
(not yet)



	se.c : old version...
