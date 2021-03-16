"""
Constants for use in calculation of constitutive equations (Passchier 
and Truow, 1996)
"""
K = 1.38062 * 10 ** -23 # Boltzmann, J/(mol K)
R = 8.3143 # Gas constant, J/(mol K)
W = 10 ** -9 # Grain boundary thickness, m
V = 2.6 * 10 ** -5 # Molar volume of solid, m^3/mol
B = 5 * 10 ** -10 # Burgers vector, m
D = 10 ** -5 # Grain size, m

MU = 42 * 10 ** 9 # Shear Modulus, N/m^2

# Diffusion constant for self diffusion, m^2/s
D_L = 2.9 * 10 ** -5 # Liquid
D_G = 3 * 10 ** -8 # Gas

# Molar activation enthalpy for self diffusion, J/mol
H_L = 243 * 10 ** 3 # Liquid
H_G = 113 * 10 ** 3 # Gas

T_M = 1550 # Melting temperature of quartz, K
A_C = 141 # Numerical Constant, Coble Creep
A_NH = 16 # Numerical Constant, Nabarro-Herring
