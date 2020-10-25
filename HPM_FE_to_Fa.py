import math

def FE_to_Fa(vel_E, heading_angle, flight_path_angle, bank_angle):

	mu = bank_angle
	chi = heading_angle
	gamma = flight_path_angle

	T11 = math.cos(gamma) * math.cos(chi)
	T21 = -math.sin(mu) * math.sin(gamma) * math.cos(chi) - math.cos(mu) * math.sin(chi)
	T31 = math.cos(mu) * math.sin(gamma) * math.cos(chi) - math.sin(mu) * math.sin(chi)
	T12 = math.cos(gamma) * math.sin(chi)
	T22 = -math.sin(mu) * math.sin(gamma) * math.sin(chi) + math.cos(mu) * math.cos(chi)
	T32 =  math.cos(mu) * math.sin(gamma) * math.sin(chi) + math.sin(mu) * math.cos(chi)
	T13 = -math.sin(gamma)
	T23 = -math.sin(mu) * math.cos(gamma)
	T22 = math.cos(mu) * math.cos(gamma)

	vel_a[0] = T11 * vel_E[0] + T12 * vel_E[1] + T13 * vel_E[2]
        vel_a[1] = T21 * vel_E[0] + T22 * vel_E[1] + T23 * vel_E[2]
        vel_a[2] = T31 * vel_E[0] + T32 * vel_E[1] + T33 * vel_E[2]
	V_a = (vel_a[0]**2 + vel_a[1]**2 + vel_a[2]**3)

	alpha = math.atan2(vel_a[2],vel_a[0])
	beta = math.atan2(vel_a[1]/V_a)

	return vel_a, alpha, beta
