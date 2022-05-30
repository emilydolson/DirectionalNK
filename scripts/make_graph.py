import itertools

N = 15

print("from, to")
for genotype in itertools.product("01", repeat=N):
    gene_num = int("0b0" + "".join(genotype), base=2)
    for i in range(N):
        new_genotype = list(genotype[:])
        new_genotype[i] = "1" if genotype[i] == "0" else "0"
        num = int("0b0" + "".join(new_genotype), base=2)
        # if num < gene_num:
        print(",".join([str(gene_num), str(int("0b0" + "".join(new_genotype), base=2))]))
