# Tools for experimenting with balanced vectors
# Part of the University of Rochester's 2015 mathematics REU
# Written by Charlotte Aten (caten2@u.rochester.edu)

# Balanced vectors over Z_p are vectors of length k*p which have k components equal to 0,
# k components equal to 1, ..., and k components equal to p-1.

import sage.all
from sage.rings.finite_rings.constructor import FiniteField
from itertools import permutations
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import Matrix

class BalancedVectors():
    """The collection of balanced vectors of size s over Z_p."""
    def __init__(self, size, modulus):
        assert size % modulus == 0 #Ensure the size is a multiple of the modulus so the set is nonempty
        self.size = size #The size of the balanced vectors
        self.modulus = modulus #The size of the underlying field
        self.quotient = self.size/self.modulus
        self.field = FiniteField(modulus) #The underlying field
    
    # Check if a vector is balanced
    def is_balanced(self, vec):
        V = list(vec)
        for i in range(self.modulus):
            if V.count(i) != self.quotient:
                return False
        return True

    # Create a generator for all balanced vectors of size s over Z_p.
    def elements(self):
        lis = []
        for permutation in permutations([self.field(i) for i in range(self.size)], self.size):
            if permutation not in lis:
                lis.append(permutation)
                yield (permutation)
    # Create a generator for all balanced vectors beginning with 0 belonging to the above set.
    def reducedelements(self):
        for elem in self.elements():
            if elem[0] == 0:
                yield (elem)

    # Take the difference of two vectors.
    def difference(self, vec0, vec1):
        # Use "0" to represent an unbalanced vector. The unbalanced vector absorbs everything.
        if vec0 == 0 or vec1 == 0:
            return 0
        # Compute the difference vector and check if it is balanced.
        dif = tuple([vec0[i] - vec1[i] for i in range(self.size)])
        if self.is_balanced(dif) == True:
            return dif
        else:
            return 0

    # Compute the (left-handed) action of a vector on the set of balanced vectors.
    def action(self, vec):
        for x in self.elements():
            yield [x, self.difference(vec, x)]

    # Compute the (left-handed) action of a vector on the subset of balanced vectors which begin with 0.
    def reducedaction(self, vec):
        for x in self.reducedelements():
            yield [x, self.difference(vec, x)]

    def reduced_digraph(self, vec):
        D = DiGraph()
        for pair in self.reducedaction(vec):
            D.add_edge(pair[0], pair[1])
        return D
    
    def trial_generators(self, vecset):
        D = DiGraph()
        for vec in vecset:
            for pair in self.reducedaction(vec):
                D.add_edge(pair[0], pair[1])
    
    def pairact(self, vec0, vec1):
        vec2 = self.difference(vec0, vec1)
        if vec2 == 0:
            print("This action is degenerate. It produces the unbalanced element.")
            return
        while vec2 != vec1:
            print(vec2)
            vec2 = self.difference(vec0, vec2)

# Print the list of all six balanced vectors of size 4 over Z_2.
S = BalancedVectors(4, 2)
for x in S.elements():
    print(str(x))

print("\n")

x = list(S.elements())[0]
print(x)
print("")

for y in S.reducedaction(x):
    print(y)
   
print(S.reduced_digraph(x).connected_components())
