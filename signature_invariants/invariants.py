"""Provides the invariants of
   `Diehl/Reizenstein - 2018 - Invariants of multidimensional time series based on their iterated-integral signature`"""

#import free_lie_algebra as fla
import numpy as np
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.generators import symmetric
import sympy
import itertools, functools
import unittest
import linear_combination.linear_combination as lc
import linear_combination.words as words

from collections import defaultdict

def signature_of_path(path,upto_level):
    """Convenience method for calculating the signature of a path using iisignature.
    
    Parameters
    ----------
    path : numpy.array
       A (T,d) dimensional array, with `d` the dimnsion of the path and `T` the number of time-steps.
    upto_level: int
       Compute the signature up-to (and including) this level.

    Returns
    -------
    signature : A `LinearCombination` of `ConcatenationWord`s storing the signature as its coefficients.
    """
    import iisignature

    d = np.shape(path)[-1]

    s = iisignature.sig(path,upto_level,1)
    letters = range(1,d+1)    
    def signature_generator():
        yield (words.ConcatenationWord(), 1)
        for lev in range(1,upto_level+1):
            for word, val in zip( itertools.product(letters,repeat=lev), s[lev-1]):
                yield (words.ConcatenationWord(word), val)
    return lc.LinearCombination.from_generator( signature_generator() )
    #out=[{Word(key):float(val) for key,val in zip(itertools.product(letters,repeat=lev),vals)} for lev,vals in zip(range(1,m+1),s)]
    #return Elt([{emptyWord:1}]+out)


###########################
#### GL INVARIANTS ########
###########################
def _gl_invariant_for_identity(dim,level):
    assert level % dim == 0
    weight = level // dim
    result = dict()
    for taus in itertools.product(symmetric(dim),repeat=weight):
        index = []
        sign = 1
        for tau in taus:
            sign *= tau.signature()
            for i in range(dim):
                index.append( (i^tau) + 1)
        result[ words.ShuffleWord(index) ] = sign
    return lc.LinearCombination(result)

def _tableau_to_permutation(yt):
    """The unique permutation sigma s.t. sigma yt = t_first (see Diehl/Reizenstein)"""
    return Permutation( list(map( lambda x: x-1, flatten( transpose(yt) ))) )

def dyck_word_to_tableau(dim, word):
    n = len(word) // dim
    yt = []
    word = np.array(word)
    R = np.array( range(1, dim * n + 1) )
    for i in range(1,dim+1):
        yt.append( list(R[ word == i ]) )
    return yt

def rectangular_standard_young_tableaux(height, width):
    """All STY's with shape (width^height) = (width, width, .., width)."""
    return map( lambda w: dyck_word_to_tableau(height, w), dyck_words(height, width) )

def transpose(yt):
    """Transpose of a Young tableau."""
    return [list(i) for i in zip(*yt)]

def is_increasing(a):
    """Is the list `a` strictly increasing?"""
    return all( map( lambda i: a[i-1] < a[i], range(1,len(a))) )

def flatten(iterable):
    """One-level flatten."""
    for row in iterable:
        for col in row:
            yield col

def is_standard_young_tableau(yt):
    flat = list(flatten(yt))
    return all( map( is_increasing, yt ) ) and all( map( is_increasing, transpose(yt) ) ) \
            and set( flat ) == set( range(1,len(flat)+1) )

def gl_invariants(dim,level):
    """Basis for the GL invariants."""
    if level % dim != 0:
        return
    else:
        weight = level // dim
        for_identity = _gl_invariant_for_identity(dim,level)
        for yt in rectangular_standard_young_tableaux(dim, weight):
            sigma = _tableau_to_permutation( yt )**(-1)
            def permute(word):
                letters = word
                permuted = tuple( letters[i^sigma] for i in range(len(letters)) )
                return words.ShuffleWord(permuted)
            yield lc.LinearCombination( {permute(x):y for x,y in for_identity.items()} )

def _dyck_words(dim,n,current_word,current_tally):
    # XXX this is very slow .. though faster than via young
    if len(current_word) == dim * n - 1:
        # only one possibility left
        yield current_word + [dim]
    else:
        options = list(filter( lambda s: current_tally[s -1] < n and (s == 1 or current_tally[s -1] < current_tally[s -2]), range(1,dim+1) ))
        chainz = []
        for s in options:
            new_tally = current_tally.copy()
            new_tally[s -1] += 1
            yield from _dyck_words(dim,n,current_word + [s], new_tally) # XXX python > 3.3

def dyck_words(dim,n):
    """Dyck words of total length 2n over the alphabet 1,2,..dim."""
    return _dyck_words(dim,n,[],[0] * dim)

def d_catalan(dim,n):
    """d-dimensional Catalan numbers.
       The number of standard Young tableaux of shape (n^dim) = (n,..,n)."""
    res = sympy.factorial( dim * n )
    for ell in range(0,dim):
        res /= sympy.ff(n+ell,n)
    return res

####################################
#### PERMUTATION INVARIANTS ########
####################################
def set_partitions( S, at_most ):
    """Partitions of S into <= at_most sets."""
    return itertools.chain( *(map( lambda m: sympy.utilities.iterables.multiset_partitions( S, m ), range( 1, at_most+1))))

def possible_draws_without_repetition( S, n ):
    """All (ordered) sequences of drawing n times from S, without repetition."""
    for subset in itertools.combinations(S,n):
        for p in symmetric(n):
            yield [ subset[i^p] for i in range(n) ]

def _partition_draw_to_word( order, p, draw ):
    result = []
    p_as_map = map( lambda i_p: map( lambda x: [x,i_p], p[i_p] ) , range(len(p)) )
    p_as_map = dict( [item for sublist in p_as_map for item in sublist] )
    for i in range(len(p_as_map)):
        result.append( draw[ p_as_map[i] ] )
    return result

def permutation_invariants(dim,level):
    """Basis for the permutation invariants."""
    partitions = set_partitions( range(level), dim )
    invs = []
    for p in partitions:
        inv = dict()
        for draw in possible_draws_without_repetition( range(1,dim+1), len(p) ):
            inv[ words.ShuffleWord(_partition_draw_to_word(level,p,draw)) ] = 1
        invs.append(lc.LinearCombination(inv))
    return invs


def restricted_permutations(n):
    """Permutations of 0,..,n-1 s.t. sigma(n+1) = sigma(n) + 1 mod n. There are n of them."""
    one_step = Permutation( list(range(1,n)) + [0] )
    for i in range(n):
        yield one_step**i

def restricted_permutation_invariants_for_partition(dim,level,p):
    """Put 1 in the first box; fill other boxes; permute all."""
    for draw in possible_draws_without_repetition( range(2,dim+1), len(p)-1 ):
        #inv = fla.Elt([]) # XXX
        inv = lc.LinearCombination()
        draw = [1] + draw
        word = _partition_draw_to_word( level, p, draw )
        for perm in restricted_permutations( dim ):
            permuted = map( lambda n: ((word[n]-1)^perm) +1, range(level) ) # XXX signature / words .. is 1-based but permutations are 0-based
            inv += lc.LinearCombination.lift( words.ShuffleWord(permuted) ) # XXX slow #fla.word2Elt( tuple(permuted) )
        yield inv

def restricted_permutation_invariants(dim,level):
    """Basis for the invariants wrt permutations of 1,..,dim that preserve order (mod dim)."""
    partitions = set_partitions( range(level), dim )
    invs = []
    for p in partitions:
        for inv in restricted_permutation_invariants_for_partition(dim,level,p):
            yield inv

###########################
#### SO INVARIANTS ########
###########################

def only_n_tuples( n ):
    def f(a):
        return all( map( lambda p: len(p) == n, a ) )
    return f

def equal_size_partitions(S,size):
    """Returns all partitions of the set `S` into sets of size `size`."""
    if len(S) == 0:
        return []
    else:
        return filter( only_n_tuples(size), sympy.iterables.multiset_partitions( S )) # XXX brute force

def _so_invariants_without_det(dim,level):
    if level % 2 == 1:
        yield from []
    else:
        for partitions in equal_size_partitions(range(0,level), 2):
            inv = defaultdict(int)
            for letters in itertools.product( range(1,dim+1), repeat=level//2 ):
                current_word = np.zeros( level, dtype=np.int8 )
                for i in range(level//2):
                    current_word[ partitions[i][0] ] = letters[i]
                    current_word[ partitions[i][1] ] = letters[i]
                inv[ words.ShuffleWord(current_word) ] += 1
            yield lc.LinearCombination(inv)

def _so_invariants_with_one_det(dim,level):
    # no determinant fits or, using one determinant, the other spots cannot be used for inner products
    if level <  dim or (level - dim) % 2 == 1:
        yield from []
    elif level == dim:
        yield _gl_invariant_for_identity(dim, dim)
    else:
        for det_positions in itertools.combinations(range(level),dim): # positions for the determinant
            rest = set(range(level))-set(det_positions)
            for partitions in equal_size_partitions(rest, 2):
                inv = defaultdict(int)
                for tau in symmetric(dim):
                    for letters in itertools.product( range(1,dim+1), repeat=level//2 ):
                        current_word = np.zeros( level, dtype=np.int8 )
                        for i in range( (level-dim)//2 ):
                            current_word[ partitions[i][0] ] = letters[i]
                            current_word[ partitions[i][1] ] = letters[i]
                        det_letters = [ (i^tau) + 1 for i in range(dim) ]
                        for i in range(dim):
                            current_word[ det_positions[i] ] = det_letters[i]
                        sign = tau.signature()
                        inv[ words.ShuffleWord(current_word) ] += sign
                yield lc.LinearCombination( inv )

def so_invariants(dim,level):
    """Linear generating set for the SO invariants."""
    return itertools.chain( _so_invariants_without_det(dim,level), _so_invariants_with_one_det(dim,level) )
