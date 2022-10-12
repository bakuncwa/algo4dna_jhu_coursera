# Reading genome files
def rGenome(filename):
    g = ""
    with open(filename, "r") as f:
        for l in f:
            if not l[0] == ">":
                g += l.rstrip()
    return g


# 1. How many alignments does the naive exact matching algorithm try when matching the string
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the
# excerpt of human chromosome 1?  (Don't consider reverse complements.)
# 2. How many character comparisons does the naive exact matching algorithm try when matching the string
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human
# chromosome 1?  (Don't consider reverse complements.)

# Naive Exact Matching algorithm with alignment reads + character comparison count
def naive_with_c(p,t):
    o = []
    y = len(t)
    x = len(p)
    alignment_c = 0
    chara_comp_c = 0
    for i in range (y-x +1):
        match = True
        alignment_c += 1
        for j in range(x):
            chara_comp_c += 1
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            o.append(i)
    return o, alignment_c, chara_comp_c

chr1_genome = rGenome("/Users/. . ./Downloads/chr1.GRCh38.excerpt.fasta")
p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
o, alignment_c, chara_comp_c = naive_with_c(p, chr1_genome)
print("Q1 + Q2. Naive Exact Matching algorithm in chr1 genome:")
print("Occurrences:", o, " | Alignment reads:", alignment_c, " | Character comparisons:", chara_comp_c)


# 3. How many alignments does Boyer-Moore try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG
# (derived from human Alu sequences) to the excerpt of human chromosome 1?  (Don't consider reverse complements.)

# Arrays
import string
def z_array(s):
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s ) -1)
    for i in range(1, len(s)):
        if s[i] == s[ i -1]:
            z[1] += 1
        else:
            break
    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1
    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[ i -k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                z[k] = zkp
            else:
                nmatch = 0
                for i in range( r +1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    lp = [0] * len(p)
    for j in range(len(p ) -1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[ i -1], lp[i])
    return l


def small_l_prime_array(n):
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+ 1: 
            small_lp[len(n) - i - 1] = i + 1
    for i in range(len(n) - 2, -1, -1):
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i + 1]
    return small_lp


# Good suffix rule
def good_suffix_table(p):
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    return len(small_l_prime) - small_l_prime[1]


# Bad character rule
def dense_bad_char_tab(p, amap):
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i + 1
    return tab


# Boyer-Moore algorithm
class BoyerMoore(object):
    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        self.amap = {}
        for i in range(len(self.alphabet)):
            self.amap[self.alphabet[i]] = i
        self.bad_char = dense_bad_char_tab(p, self.amap)
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        assert c in self.amap
        ci = self.amap[c]
        assert i > (self.bad_char[i][ci] - 1)
        return i - (self.bad_char[i][ci] - 1)

    def good_suffix_rule(self, i):
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        return len(self.small_l_prime) - self.small_l_prime[1]
    
    
# Boyer-Moore algorithm with alignment reads + character comparison count
def bm_with_c(p, p_bm, t):
    i = 0
    o = []
    y = len(t)
    x = len(p)
    bm_align_c = 0
    bm_chara_c = 0
    while i < y-x +1:
        bm_align_c += 1
        shift = 1
        mismatched = False
        for j in range(x-1, -1, -1):
            bm_chara_c += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            o.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return o, bm_align_c, bm_chara_c
    
chr1_genome = rGenome("/Users/. . ./Downloads/chr1.GRCh38.excerpt.fasta")
p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
p_bm = BoyerMoore("GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG")
o, bm_align_c, bm_chara_c = bm_with_c(p, p_bm, chr1_genome)
print("Q3. Boyer-Moore algorithm in chr1 genome:")
print("Occurrences:", o, " | Alignment reads:", bm_align_c, " | Character comparisons:", bm_chara_c)


# 4. Index-assisted approximate matching. In practicals, we built a Python class called Index
# implementing an ordered-list version of the k-mer index.  The Index class is copied below.
# We also implemented the pigeonhole principle using Boyer-Moore as our exact matching algorithm.
# Implement the pigeonhole principle using Index to find exact matches for the partitions.
# Assume P always has length 24, and that we are looking for approximate matches with up to 2 mismatches (substitutions).
# We will use an 8-mer index.
# Download the Python module for building a k-mer index.
# https://d28rh4a8wq0iu5.cloudfront.net/ads1/code/kmer_index.py
# Write a function that, given a length-24 pattern P and given an Index object built on 8-mers, finds all approximate
# occurrences of P within T with up to 2 mismatches. Insertions and deletions are not allowed. Don't consider any
# reverse complements.
# How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, which is derived from a human Alu sequence, occur with up
# to 2 substitutions in the excerpt of human chromosome 1?  (Don't consider reverse complements here.)
# Hint 1: Multiple index hits might direct you to the same match multiple times, but be careful not to count
# a match more than once.
# Hint 2: You can check your work by comparing the output of your new function to that of the naive_2mm
# function implemented in the previous module.

import bisect
class Index(object):
    def __init__(self, t, k):
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1): 
            self.index.append((t[i:i + k], i)) 
        self.index.sort()

    def query(self, p):
        kmer = p[:self.k]
        i = bisect.bisect_left(self.index, (kmer, -1))
        hits = []
        while i < len(self.index):  
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
    

def approximate_match_2mm(p, t, n):
    y = len(t)
    x = len(p)
    segment_length = int(round(x / (n+1)))
    all_matches = set()
    idx_p = Index(t, segment_length)
    idx_hits = 0
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, x)
        matches = idx_p.query(p[start:end])
        for m in matches:
            idx_hits += 1
            if m < start or m-start+x > y:
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n: 
                        break
            for j in range(end, x):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches), idx_hits

chr1_genome = rGenome("/Users/. . ./Downloads/chr1.GRCh38.excerpt.fasta")
p = "GGCGCGGTGGCTCACGCCTGTAAT"
n = 2
print("Q4. Index hits up to 2 matches:")
print("Index hits:", len(approximate_match_2mm(p, chr1_genome, n)[0]))


# Naive Exact Matching algorithm (up to 2 mismatches) output comparison
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

chr1_genome = rGenome("/Users/. . ./Downloads/chr1.GRCh38.excerpt.fasta")
p = "GGCGCGGTGGCTCACGCCTGTAAT"
o = naive_2mm(p, chr1_genome)
print("Naive_2mm up to 2 matches in chr1 genome:")
print("Naive_2mm: %d" % len(o))


# Q5.Using the instructions given in Question 4, how many total index hits are there when searching for occurrences of
# GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of human chromosome 1?
# (Don't consider reverse complements.)
# Hint: You should be able to use the boyer_moore function (or the slower naive function) to double-check your answer.

chr1_genome = rGenome("/Users/. . ./Downloads/chr1.GRCh38.excerpt.fasta")
p = "GGCGCGGTGGCTCACGCCTGTAAT"
n = 2
print("Q5. Total index hits up to 2 matches:")
print("Total index hits:", (approximate_match_2mm(p, chr1_genome, n)[1]))


# 6. Let's examine whether there is a benefit to using an index built using subsequences of T rather than substrings, as we
# discussed in the "Variations on k-mer indexes" video.  We'll consider subsequences involving every N characters.
# For example, if we split ATATAT\verb|ATATAT|ATATAT into two substring partitions, we would get partitions ATA (the first half)
# and TAT (second half).  But if we split ATATAT into two  subsequences  by taking every other character, we would get
# AAA (first, third and fifth characters) and TTT\verb|TTT|TTT (second, fourth and sixth).
# Another way to visualize this is using numbers to show how each character of P is allocated to a partition.
# Splitting a length-6 pattern into two substrings could be represented as 111222, and splitting into two subsequences
# of every other character could be represented as 121212
# The following class SubseqIndex is a more general implementation of Index that additionally handles subsequences.
# It only considers subsequences that take every Nth character:

import bisect
class SubseqIndex(object):
    
    def __init__(self, t, k, ival):
        self.k = k  
        self.ival = ival  
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):
            self.index.append((t[i:i + self.span:ival], i)) 
        self.index.sort()

    def query(self, p):
        subseq = p[:self.span:self.ival]
        i = bisect.bisect_left(self.index, (subseq, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

# Write a function that, given a length-24 pattern P and given a SubseqIndex object built
# with k = 8 and ival = 3, finds all approximate occurrences of P within T with up to 2 mismatches.
# When using this function, how many total index hits are there when searching for GGCGCGGTGGCTCACGCCTGTAAT with up to 2
# substitutions in the excerpt of human chromosome 1?  (Again, don't consider reverse complements.)
# Hint: See this notebook for a few examples you can use to test your function.

def approximate_match_subseq(p, t, n, ival):
    y = len(t)
    x = len(p)
    segment_length = int(round(x / (n+1)))
    all_matches = set()
    idx_p = SubseqIndex(t, segment_length, ival)
    idx_hits = 0 
    for i in range(n+1):
        start = i
        matches = idx_p.query(p[start:])
        for m in matches:
            idx_hits += 1
            if m < start or m-start+x > y:
                continue  
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break             
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches), idx_hits

chr1_genome = rGenome("/Users/. . ./Downloads/chr1.GRCh38.excerpt.fasta")
p = "GGCGCGGTGGCTCACGCCTGTAAT"
n = 2
ival = 3
print("Q6. Total subsequence index hits up to 2 matches:")
print("Total subsequence index hits:", (approximate_match_subseq(p, chr1_genome, n, ival)[1]))
