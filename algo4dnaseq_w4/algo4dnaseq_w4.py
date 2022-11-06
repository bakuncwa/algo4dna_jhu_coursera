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
