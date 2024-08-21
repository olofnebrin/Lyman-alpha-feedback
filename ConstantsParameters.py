""" Some useful physical constants, parameters, and conversion factors. """

import numpy as np

# Fundamental constants:

c       =  2.9979245e10   # Speed of ligth (cm/s)
kB      =  1.3806488e-16  # Boltzmann constant (erg/K)
mH      =  1.6735327e-24  # Mass of hydrogen atom (g)

# Einstein-A coefficients:

A_Lya   =  6.265e8        # Lya Einstein-A coefficient (s^-1)
A_2ph   =  8.2206         # 2s->1s Einstein-A coefficient(s^-1)
A_2p2s  =  8.78e-7        # 2p->2s Einstein-A coefficient(s^-1)
A_2s2p  =  1.597e-9       # 2s->2p Einstein-A coefficient(s^-1)

# Hydrogen & helium abundance:

Y       =  0.247          # Helium mass fraction
X       =  1.0-Y          # Hydrogen mass fraction

# Conversion factors:

kms     =  1e5            # 1 km/s in cm/s
