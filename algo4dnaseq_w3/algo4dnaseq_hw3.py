# Reading genome files
def rGenome(filename):
    g = ""
    with open(filename, "r") as f:
        for l in f:
            if not l[0] == ">":
                g += l.rstrip()
    return g


# Edit Distance algorithm
def editDistance(p, t):
    D = []
    x = len(p)
    y = len(t)
    for i in range(x+1):
        D.append([0]*(y+1))
    for i in range(x+1):
        D[i][0] = i
    for i in range(y + 1):
        D[0][i] = i
    for i in range(1, x + 1):
        for j in range(1, y + 1):
            distHor = D[i][j - 1] + 1
            distVer = D[i - 1][j] + 1
            if p[i - 1] == t[j - 1]:
                distDiag = D[i - 1][j - 1]
            else:
                distDiag = D[i - 1][j - 1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    return D[-1][-1]

chr1_genome = rGenome("/Users/. . ./Downloads/chr1.GRCh38.excerpt.fasta")
x = "GCTGATCGATCGTACG"
print("Edit distance value:", editDistance(x, chr1_genome))


# 1. What is the edit distance of the best match between pattern GCTGATCGATCGTACG and the excerpt of
# human chromosome 1?  (Don't consider reverse complements.)
# 2. What is the edit distance of the best match between pattern GATTTACCAGATTGAG and the excerpt of
# human chromosome 1?  (Don't consider reverse complements.)

# Edit Distance algorithm + Best Approximate Match
def editDistBestApproximateMatch(p, t):
    D = []
    x = len(p)
    y = len(t)
    for i in range(x+1):
        D.append([0]*(y+1))
    for i in range(x+1):
        D[i][0] = i
    for i in range(1, x + 1):
        for j in range(1, y + 1):
            distHor = D[i][j - 1] + 1
            distVer = D[i - 1][j] + 1
            if p[i - 1] == t[j - 1]:
                distDiag = D[i - 1][j - 1]
            else:
                distDiag = D[i - 1][j - 1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    return min(D[-1])

chr1_genome = rGenome("/Users/. . ./Downloads/chr1.GRCh38.excerpt.fasta")
p = "GCTGATCGATCGTACG"
print("Q1. Best approximate edit distance value in chr1 genome:")
print("Edit distance value:", editDistBestApproximateMatch(p, chr1_genome))

chr1_genome = rGenome("/Users/. . ./Downloads/chr1.GRCh38.excerpt.fasta")
p = "GATTTACCAGATTGAG"
print("Q2. Best approximate edit distance value in chr1 genome:")
print("Edit distance value:", editDistBestApproximateMatch(p, chr1_genome))
