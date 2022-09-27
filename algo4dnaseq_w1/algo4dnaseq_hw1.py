# Reading genome files
def rGenome(filename):
    g = ""
    with open(filename, "r") as f:
        for l in f:
            if not l[0] == ">":
                g += l.rstrip()
    return g


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


# 1. How many times does AGGT or its reverse complement (ACCT) occur in the lambda virus genome?
# E.g. if AGGT occurs 10 times and ACCT occurs 12 times, you should report 22.

# Naive Exact Matching algorithm
def naive(p, t):
    o = []
    y = len(t)
    x = len(p)
    for i in range(y - x + 1):
        match = True
        for j in range(x):
            if t[i + j] != p[j]:
                match = False
                break
        if match:
            o.append(i)
    return o

# Finding reverse complement of a DNA strand
def reverse_comp(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    DNA_strand = " "
    for nuc_base in seq[::-1]:
        DNA_strand = DNA_strand + comp[nuc_base]
    return DNA_strand

# Naive exact matching algorithm that is strand aware
# Instead of looking only for occurrences of P in T, look for occurrences of the reverse complement of P in T

def naive_with_rc(p, t):
    o_naive = naive(p,t)
    o_naive_rc = naive(reverse_comp(p), t)
    rc = reverse_comp(p)
    if rc == p:
        return o_naive
    elif rc != p:
        return o_naive + o_naive_rc

lambda_virus_genome = rGenome("/Users/. . ./Downloads/lambda_virus.fa")
p = "ACCT"
o = naive_with_rc(p, lambda_virus_genome)
print("ACCT in lambda_virus_genome:")
print("# Occurrences: %d" % len(o))

lambda_virus_genome = rGenome("/Users/. . ./Downloads/lambda_virus.fa")
p = "AGGT"
o = naive_with_rc(p, lambda_virus_genome)
print("AGGT in lambda_virus_genome:")
print("# Occurrences: %d" % len(o))


# 2. How many times does TTAA or its reverse complement occur in the lambda virus genome?
# Hint: TTAA and its reverse complement are equal, do not double count
# if P and its reverse complement  are identical; given match offset should be reported only once

def naive_2mm(p, t):
    o = []
    y = len(t)
    x = len(p)
    for i in range(y-x + 1):
        match = True
        mismatch_c = 0
        for j in range(x):
            if t[i+j] != p[j]:
                mismatch_c+=1
                if mismatch_c > 2:
                    match = False
                    break
        if match:
            o.append(i)
    return o

lambda_virus_genome = rGenome("/Users/. . ./Downloads/lambda_virus.fa")
p = "TTAA"
o = naive_with_rc(p, lambda_virus_genome)
print("TTAA in lambda_virus_genome:")
print("# Occurrences: %d" % len(o))


# 3. What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement in the Lambda virus genome?
# E.g. if the leftmost occurrence of ACTAAGT is at offset 40 (0-based) and the leftmost occurrence of its reverse complement
# ACTTAGT is at offset 29, then report 29.

lambda_virus_genome = rGenome("/Users/. . ./Downloads/lambda_virus.fa")
p = "ACTAAGT"
o = naive_with_rc(p, lambda_virus_genome)
print("ACTAAGT in lambda_virus_genome:")
print("# Occurrences: %d" % min(o))

lambda_virus_genome = rGenome("/Users/. . ./Downloads/lambda_virus.fa")
p = "ACTTAGT"
o = naive_with_rc(p, lambda_virus_genome)
print("ACTTAGT in lambda_virus_genome:")
print("# Occurrences: %d" % min(o))


# 4. What is the offset of the leftmost ocurrence of AGTCGA or its reverse complement in the Lambda virus genome?

lambda_virus_genome = rGenome("/Users/. . ./Downloads/lambda_virus.fa")
p = "AGTCGA"
o = naive_with_rc(p, lambda_virus_genome)
print("AGTCGA in lambda_virus_genome:")
print("# Occurrences: %d" % min(o))

lambda_virus_genome = rGenome("/Users/. . ./Downloads/lambda_virus.fa")
p = "AGTGGA"
o = naive_with_rc(p, lambda_virus_genome)
print("AGTGGA in lambda_virus_genome:")
print("# Occurrences: %d" % min(o))


# 5. As we will discuss, sometimes we would like to find approximate matches for P in T. That is, we want to find
# occurrences with one or more differences. For Questions 5 and 6, make a new version of the naive function called
# naive_2mm that allows up to 2 mismatches per occurrence. Unlike for the previous questions, do not consider the
# reverse complement here.  We're looking for approximate matches for P itself, not its reverse complement.
# For example, ACTTTA occurs twice in ACTTACTTGATAAAGT, once at offset 0 with 2 mismatches, and once at offset 4
# with 1 mismatch. So naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT') should return the list [0, 4].
# Hint: See  this notebook for a few examples you can use to test your naive_2mm function.
# How many times does TTCAAGCC occur in the Lambda virus genome when allowing up to 2 mismatches

lambda_virus_genome = rGenome("/Users/. . ./Downloads/lambda_virus.fa")
p = "TTCAAGCC"
o = naive_2mm(p, lambda_virus_genome)
print("TTCAAGCC (up to 2 matches) in lambda_virus_genome:")
print("# Occurrences: %d" % len(o))


# 6. What is the offset of the leftmost occurrence of the AGGAGGTT in the Lambda virus genome when allowing
# up to 2 matches?

lambda_virus_genome = rGenome("/Users/. . ./Downloads/lambda_virus.fa")
p = "AGGAGGTT"
o = naive_2mm(p, lambda_virus_genome)
print("AGGAGGTT (up to 2 matches) in lambda_virus_genome:")
print("# Occurrences: %d" % min(o))


# 7. Finally, download and parse the provided FASTQ file containing real DNA sequencing reads derived from a human:
# https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR037900_1.first1000.fastq
# Note that the file has many reads in it and you should examine all of them together when answering this question.
# The reads are taken from this study:
# Ajay, S. S., Parker, S. C., Abaan, H. O., Fajardo, K. V. F., & Margulies, E. H. (2011).
# Accurate and comprehensive sequencing of personal genomes. Genome research, 21(9), 1498-1505.
# This dataset has something wrong with it; one of the sequencing cycles is poor quality.
# Report which sequencing cycle has the problem.  Remember that a sequencing cycle corresponds to a particular offset
# in all the reads. For example, if the leftmost read position seems to have a problem consistently across reads, report 0.
# If the fourth position from the left has the problem, report 3. Do whatever analysis you think is needed to identify
# the bad cycle. It might help to review the "Analyzing reads by position" video.

# def createHist(qualities):
   # sequences, qualities = rFastq("/Users/. . ./Downloads/ERR037900_1.first1000.fastq")
   # phredscore = []
   # qualities = qualities[-5]
   # for phred in quals:
       # q = ord(phred)-33
       # phredscore.append(q)
   # return phredscore
# h = createHist(qualities)
# print(h)
# import matplotlib.pyplot as plt
# plt.bar(range(len(h)), h)
# plt.show()
# min_h = min(h)
# index_min_h = h.index(min_h)
# print(index_min_h)
