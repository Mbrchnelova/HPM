import sys


verysmall = 1e-20
# For this interpolation, the basearray must be ordered in an increasing fashion!
def linear_interpolate(basevalue, basearray, targetarray):

	if len(basearray) > 1:
		if basearray[1] < basearray[0]:
			sys.exit("Incorrect base array ordering for interpolation. Aborting.")

	if len(basearray) == 1:
		sys.exit("Hard to interpolate an array of length of one. Get your shit together. Aborting.")

	idx_lower = 0
	idx_higher = len(basearray)

	for i in range(0, len(basearray)):
		if basevalue >= basearray[i]:
			idx_lower = max(0, i)
			idx_higher = min(len(basearray), i + 1)
			break

	if idx_lower == idx_higher:
		print "Warning - interpolation cannot be performed as the array is outside of the bounds."
		print "Performing extrapolation instead."

		targetvalue = linear_extrapolate(basevalue, basearray, targetarray)
	else:
		d1 = max(verysmall, abs(basevalue - basearray[idx_lower]))
        	d2 = max(verysmall, abs(basevalue - basearray[idx_higher]))

		w1 = 1./d1
		w2 = 1./d2

		targetvalue = (targetarray[idx_lower] * w1 + targetarray[idx_higher] * w2) / (w1 + w2) 

	return targetvalue


def linear_extrapolate(basevalue, basearray, targetarray):


	if basevalue > basearray[-1]:
		slope = (targetarray[-1] - targetarray[-2]) / (basearray[-1] - basearray[-2])

		targetvalue = targetarray[-2] + slope * (basevalue - basiearray[-2])

	else:
		slope = (targetarray[1] - targetarray[0]) / (basearray[1] - basearray[0])

		targetvalue = targetvalue[1] - slope * (basiearray[1] - basevalue)


	return targetvalue
