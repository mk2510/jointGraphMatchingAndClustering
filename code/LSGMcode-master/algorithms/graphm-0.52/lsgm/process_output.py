import sys 

# number of lines to delay be for copying
delay = 1
# copy matrix of the output file
with open(sys.argv[1]) as infile, open(sys.argv[2], 'w') as outfile:
	copy = False
	for line in infile:
		if line.strip() == "Permutations:":
			copy = True
			cnt = 0
		elif line.strip() == "Experiment parameters:":
			copy = False
		elif copy:
			cnt += 1
			if cnt > delay:
				outfile.write(line)

outfile.close()
