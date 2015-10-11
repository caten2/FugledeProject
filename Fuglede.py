# Tools for testing Fuglede's conjecture
# Part of the University of Rochester's 2015 mathematics REU
# Written by Charlotte Aten (caten2@u.rochester.edu)

# This program uses sagemath extensively.
import sage.all
from sage.rings.finite_rings.constructor import FiniteField
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import Matrix

# The class containing all of the tests we need to perform.
class SubsetPairs:
    """The collection of all pairs (E,B) where E and B are subsets of Z_p^d."""
    def __init__(self, modulus, dimension):
        self.modulus = int(modulus) #The modulus p
        self.dimension = int(dimension) #The dimension d
        self.field = FiniteField(modulus) #The underlying field
        self.space = MatrixSpace(self.field, 1, dimension) #The space of vectors in Z_p^d viewed as matrices
        self.elements = list(self.space) #An actual list of those vectors
        self.basis = self.space.basis() #The standard basis for the vector space
    # Preliminary construction of the possible subsets of size "size", implemented as a python generator.
    # This is a huge speed bottleneck. Right now I just have it throwing out subsets if their vectors are not "in order" so
    # permutations of the same set are eliminated. There are faster ways to do this, but I haven't found one that is elegant yet.
    def subsets(self, size):
        for elem in MatrixSpace(self.field, self.dimension, size):
            passing = True
            for i in range(size - 1):
                if self.elements.index(Matrix(elem.transpose()[i])) >= self.elements.index(Matrix(elem.transpose()[i+1])):
                    passing = False
                    break
            if passing == True:
                yield(elem)
    # Create the sets of the first type which contain the zero vector and the standard basis vectors.
    # Note that these are not the only type of vectors which need to be checked in general, 
    # so this can only find certain types of counterexamples.
    def firstsets(self, size):
        for elem in self.subsets(size):
            if elem.columns()[0] == 0 * elem.columns()[0] and elem.rref() == elem and elem.rank() == min(size, self.dimension):
                yield(elem)
    # Create the sets of the second type which contain the zero vector.
    def secondsets(self, size):
        for elem in self.subsets(size):
            if elem.columns()[0] == 0 * elem.columns()[0]:
                yield(elem)
    # Create the log-Hadamard matrix for a given pair of subsets.
    def loghadamard(self, E, B):
        return E.transpose() * B
    # Check if the difference of two rows is balanced.
    def row_difference(self, row0, row1):
        difference_vector = list(row0 - row1)
        counters = {}
        for i in range(self.modulus):
            counters[i] = 0
        for entry in difference_vector:
            value = difference_vector[entry]
            counters[value] = counters[value] + 1
        for i in range(self.modulus - 1):
            if counters[i] != counters[i+1]:
                return False
        return True
    # Perform the appropriate tests on all subsets not already eliminated.
    # At the moment this still tests (B,E) even after (E,B) has been tested.
    def runtest(self, size):
        for elem0 in S.firstsets(size):
            for elem1 in S.secondsets(size):
                H = self.loghadamard(elem0, elem1)
                passing = True
                for i,j in [(i,j) for i in range(size) for j in range(size)]:
                    if i != j and self.row_difference(H[i], H[j]) == False:
                        print(H[i], H[j], H[i]-H[j]) #Comment if you don't want to watch pairs of vectors and their differences pour down the screen.
                        passing = False
                        break
                if passing == True:
                    print(elem0, elem1, H)
                    print("")

print("Example for sets of size 4 in Z_2^3.")
S = SubsetPairs(2,3)
print("")
print("Here are subsets of the first type.")
for elem in S.firstsets(4):
    print(elem)
    print("")
print("Here are subsets of the second type.")
for elem in S.secondsets(4):
    print(elem)
    print("")
print("And finally we look for spectral pairs of size 4.")
S.runtest(4)
print("Of course in this case we found nothing but nonspectral pairs.")
print("It is easy to modify the test function to only speak up when a pair is found.")