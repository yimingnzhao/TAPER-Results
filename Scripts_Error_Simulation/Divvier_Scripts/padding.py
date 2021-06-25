import sys

USAGE = "python getErrorRates.py [error file] [Divvier output file]"

if not len(sys.argv) == 3:
	print()
	print("\tIncorrect Parameters")
	print("\tUSAGE: " + USAGE)
	sys.exit()
	
error = open(sys.argv[1], "r").read().split("\n")
input = open(sys.argv[2], "r").read().split("\n")
output = [[c for c in error[i]] if i % 2 == 0 else ["-" for c in error[i]] for i in range(len(error))]

def subcolumn(i, j):
	for k in range(int(len(error) / 2)):
		if error[2*k+1][i] != input[2*k+1][j] and input[2*k+1][j] != "-":
			return False
	return True

def replacecolumn(i, j):
	for k in range(int(len(error) / 2)):
		output[2*k+1][i] = input[2*k+1][j]
	
y = 0
for x in range(len(error[1])):
	if y < len(input[1]) and subcolumn(x, y):
		replacecolumn(x, y)
		y += 1

print("\n".join(["".join(line) for line in output]))