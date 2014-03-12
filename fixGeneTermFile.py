import sys
with open(sys.argv[1]) as f:
	data = f.readlines()
for line in data:
	if "\t" in line.strip():
		print line.strip().replace("\t"," ")
	else:
		print line.strip()
