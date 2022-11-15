# Overlap function
def overlap(a, b, min_length=3):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1


# Shortest Common Superstring algorithm
import itertools

def scs(ss):
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]
        for i in range(len(ss)-1):
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup
    return shortest_sup


# 1. In a practical, we saw the SCS function (copied below along with Overlap) for finding the shortest common
# superstring of a set of strings. It's possible for there to be multiple different shortest common superstrings for the
# same set of input strings. Consider the input strings ABC, BCA, CAB. One shortest common superstring is
# ABCAB but another is BCABC and another is CABCA.
# What is the length of the shortest common superstring of the following strings?
# CCT, CTT, TGC, TGG, GAT, ATT

print("1. SCS length of strings CCT, CTT, TGC, TGG, GAT, ATT:")
print(len(scs(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])))


# 2. How many different shortest common superstrings are there for the input strings given in the previous question?
# Hint 1: You can modify the SCS function to keep track of this.
# Hint 2: You can look at these examples to double-check that your modified SCS is working as expected.


# SCS list function
import itertools

def scs_list(ss):
    import itertools
    shortest_sup = []
    shortest_length = 0
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]
        for i in range(len(ss)-1):
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            sup += ssperm[i+1][olen:]
        if shortest_length == 0:
            shortest_sup.append(sup)
            shortest_length = len(sup)
        elif len(sup) < shortest_length:
            shortest_length = len(sup)
            shortest_sup = [sup]
        elif len(sup) == shortest_length:
            shortest_sup.append(sup)
        else:
            None
    return sorted(shortest_sup)

shortest_list = scs_list(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])
print("2. How many different shortest common superstrings are there for the input strings given in the previous question?")
print(len(shortest_list))


# 3. Download this FASTQ file containing synthetic sequencing reads from a mystery virus:
# https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ads1_week4_reads.fq
# All the reads are the same length (100 bases) and are exact copies of substrings from the forward strand of the virus
# genome.  You don't have to worry about sequencing errors, ploidy, or reads coming from the reverse strand.
# Assemble these reads using one of the approaches discussed, such as greedy shortest common superstring.  Since there
# are many reads, you might consider ways to make the algorithm faster, such as the one discussed in the programming
# assignment in the previous module.

# How many As are there in the full, assembled genome?
# Hint: the virus genome you are assembling is exactly 15,894 bases long

import operator

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

# Maximal overlap function
def pick_maximal_overlap(reads, k):
    reada, readb = None, None
    best_olen = 0
    for a, b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen

# Greedy Shortest Common Superstring algorithm
def greedy_scs(reads, k):
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)

# Revised Greedy SCS algorithm
def revised_greedy_scs(reads, k):
    pairs_olen, pairs_c = overlap_graph(reads, k)
    pairs_olen_sort = sorted(pairs_olen.items(), key=operator.itemgetter(1), reverse=True)
    read_a, read_b, olen = pairs_olen_sort[0][0][0], pairs_olen_sort[0][0][1], pairs_olen_sort[0][1]
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        pairs_olen, pairs_count = overlap_graph(reads, k)
        if pairs_olen != {}:
            pairs_olen_sort = sorted(pairs_olen.items(), key=operator.itemgetter(1), reverse=True)
            read_a, read_b, olen = pairs_olen_sort[0][0][0], pairs_olen_sort[0][0][1], pairs_olen_sort[0][1]
        elif pairs_olen == {}:
            read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return "".join(reads)

ADS1_reads = ("/Users/. . ./Downloads/ads1_week4_reads.fq")
reads, _ = rFastq(ADS1_reads)
assemble_genome = revised_greedy_scs(reads, 10)
print("Q3. Total of As in the full assembled genome:")
print("Base count:", assemble_genome.count("A"))

# 4. How many Ts are there in the full, assembled genome from the previous question?

ADS1_reads = ("/Users/. . ./Downloads/ads1_week4_reads.fq")
reads, _ = rFastq(ADS1_reads)
assemble_genome = revised_greedy_scs(reads, 10)
print("Q3. Total of Ts in the full assembled genome:")
print("Base count:", assemble_genome.count("T"))
