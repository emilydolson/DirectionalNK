import sys
import itertools

N = 6

name = sys.argv[1]
pars = name.split("-")
k = int(pars[2].strip("k"))
peak = []

landscape = [[0 for i in range(N)] for j in range(2**(k+1))]

with open(name) as infile:
    for line in infile:
        if line.startswith("# Fitness"):
            break
    for pos in range(2**(k+1)):
        line = infile.readline()
        sline = line.strip().split()
        for i in range(N):
            landscape[pos][i] = float(sline[i])

    infile.readline()
    peak = infile.readline().strip().split()

peak[0] = int(peak[0])
peak[1] = float(peak[1])

print("genotype, ones, fitness, prop_max")
for genotype in itertools.product("01", repeat=N):
    fitness = 0
    extended = genotype + genotype
    ones = genotype.count("1")
    for i in range(N):
        # print(i, i+k, "0b0" + "".join(extended[i:i+k+1]))
        num = int("0b0" + "".join(extended[i:i+k+1]), base=2)
        fitness += landscape[num][i]
    fitness /= N
    print(",".join([str(int("0b0" + "".join(genotype), base=2)), str(ones), str(fitness), str(fitness/peak[1])]))


