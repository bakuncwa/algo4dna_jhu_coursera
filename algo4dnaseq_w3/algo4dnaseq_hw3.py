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


# 3. - Concerned only with overlaps that are (a) exact matches (no differences allowed), (b) are at least k bases long
# - Calling current overlap function would be too slow
# - Where k = 6, we have a read 'a' whose length 6-suffix is GTCCTA
# - GTCCTA does NOT occur in any other read in the dataset
# - 6-mer GTCCTA occurs at the END of read 'a' and nowhere else
# - It follows that 'a''s suffix cannot possibly overlap the prefix of any other read by 6 or more characters
# - Alt: if we want to find the overlaps involving a suffix of read 'a' and a prefix of some other read,
# we can ignore any reads that don't contain the length-k suffix [HENCE THIS IS THE BOOLEAN FUNCTION, if len(y) > 6, mismatch = true)
# * Let every k-mer in the dataset have an associated Python set object, which starts out empty
# * Use a python dictionary to associate each k-mer with its corresponding set
# * (1) For every k-mer in the read, add the read to the set object corresponding to that k-mer
# If our read is GATTA and k=3, we would add GATTA to the set objects for GAT, ATT, TTA (codons/amino acids)
# We do this for every read so that, at the end, each set contains all reads containing the corresponding k-mer
# * (2) For each read 'a', we find all overlaps involving suffix of 'a'.
# Take 'a''s length-k suffix, find all reads containing that k-mer (obtained from the corresponding set) and call
# overlap(a, b, min_length=k)
# DO NOT call overlap function if 'b' does not contain the length-k suffix of 'a'
# Download and parse the read sequences from the provided Phi-X FASTQ file. Just use the base sequences
# ignore read names and base qualities. Also, no two reads in the FASTQ have the same sequence of bases

# Next, find all pairs of reads with an exact suffix/prefix match of length at least 30. Don't overlap a read
# with itself; if a read has a suffix/prefix match to itself, ignore that match.  Ignore reverse complements.

# Hint 1: Your function should not take much more than 15 seconds to run on this 10,000-read dataset, and maybe much
# less than that.  (Our solution takes about 3 seconds.) If your function is much slower, there is a problem somewhere.
# Hint 2: Remember not to overlap a read with itself. If you do, your answers will be too high.
# Hint 3: You can test your implementation by making up small examples, then checking that (a) your implementation runs
# quickly, and (b) you get the same answer as if you had simply called overlap(a, b, min_length=k) on every pair of reads.
# We also have provided a couple examples you can check against.

# Picture the overlap graph corresponding to the overlaps just calculated.  How many edges are in the graph?
# In other words, how many distinct pairs of reads overlap?

# Parses the read and quality strings from FASTQ file containing sequencing reads
def rFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fx:
        while True:
            fx.readline()
            seq = fx.readline().rstrip()
            fx.readline()
            qual = fx.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities
# sequencess, qualities = rFastq("/Users/. . ./Downloads/lambda_virus.fa")

# Overlap graph function
def overlap(a, b, min_length=6):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1

# Overlap all pairs function
def overlap_graph(reads, k):
    olaps = {}
    res = {}
    for r in reads:
        q = len(r)-k+1
        for i in range(q):
            if r[i:i+k] not in olaps:
                olaps[r[i:i+k]] = [r]
            elif r[i:i+k] in olaps:
                olaps[r[i:i+k]].append(r)
    c = 0
    for r in reads:
        r_suffix = r[-k:]
        for probable_r in olaps[r_suffix]:
            if probable_r != r:
                olen = overlap(r, probable_r, k)
                if olen > 0:
                    c += 1
                    res[(r, probable_r)] = olen
    return res, c

ERR266411_1_genome = ("/Users/. . ./Downloads/ERR266411_1.for_asm.fastq")
reads, _ = rFastq(ERR266411_1_genome)
pairs_olen, pairs_c = overlap_graph(reads, 30)
print("Q3. Number of distinct pairs of reads overlap in ERR266411_1 genome:")
print("Distinct overlap pairs:", pairs_c)


# 4. Picture the overlap graph corresponding to the overlaps computed for the previous question.
# How many nodes in this graph have at least one outgoing edge?
# (In other words, how many reads have a suffix involved in an overlap?)

reads_overlapped = []
for k, v in pairs_olen:
    reads_overlapped.append(k)
print("Q4. Number of reads that have a suffix involved in an overlap of ERR266411_1 genome:")
print("Reads involved:", len(set(reads_overlapped)))



