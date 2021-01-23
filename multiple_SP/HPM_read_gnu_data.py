import sys

def read_gnu_data(path):

	f = open(path, 'r')
	lines = f.readlines()
	line_length = len((lines[0]).split())
	no_lines = len(lines)
	for line in lines:
		if (line_length != len(line.split())):
			#print line_length
			#print  len(line.split())
			sys.exit("Corrupted data. Aborting.")

	data = []
	for j in range(0, no_lines):
                data.append([])

		for i in range(0, line_length):
			data[j].append([])

	for j in range(0, no_lines):
		line = (lines[j]).split()
		for i in range(0, line_length):
			#print line[i]
			data[j][i] = float(line[i])
			
	return data

#path = "./Output/2/out01.gnu"
#read_gnu_data(path) 
