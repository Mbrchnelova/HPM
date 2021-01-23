from HPM_numerics import *

def get_Sutherland_viscosity(T):
	
	mu_ref = 1.716e-5
	T_ref = 273.15
	S = 110.4

	mu = mu_ref * (T/ T_ref)**(3./2.) * (T_ref + S) / (T + S)

	return mu


def ISA_from_height(height):

	R = 287.058

	# Heights in m as measured in geopotential altitude above MSL
	heights = [-610., 11000., 20000., 32000., 47000., 51000., 71000., 84852.]
	
	# Base temperatures in C
	temperatures_C = [19.0, -56.5, -56.5, -44.5, -2.5, -2.5, -58.5, -86.28]

	temperatures_K = []
	for i in range(0, len(temperatures_C)):
		temperature_K = temperatures_C[i] + 273.15
		temperatures_K.append(temperature_K)

	# Base atmospheric pressure in Pa
	pressures = [108900.0, 22632., 5474.9, 868.02, 110.91, 66.939, 3.9564, 0.3734]

 	# Calculate densities from the ideal gas law
	densities = []
	for i in range(0, len(pressures)):
		density = pressures[i] / R / temperatures_K[i]
		densities.append(density)

	p = linear_interpolate(height, heights, pressures)
        T = linear_interpolate(height, heights, temperatures_K)
        rho = linear_interpolate(height, heights, densities)
	mu = get_Sutherland_viscosity(T)
	
	return p, T, rho, mu
