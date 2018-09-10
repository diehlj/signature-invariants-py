from signature_invariants.invariants import *
import pytest

def test_dyck_word_to_tableau():
    assert [ [1,2], [3,4] ] == dyck_word_to_tableau( 2, [1,1,2,2] )
    assert [ [1,3], [2,4] ] == dyck_word_to_tableau( 2, [1,2,1,2] )

def test_rectangular_standard_young_tableaux():
    for dim in [2,3]:
        for n in range(3,6):
            yts = rectangular_standard_young_tableaux(dim,n)
            assert len(list(yts)) == d_catalan(dim,n)
            for yt in yts:
                    assert is_standard_young_tableau( yt )

def test_transpose():
    assert [[1,2],[3,4]] == transpose( [[1,3],[2,4]] )

def test_is_increasing():
    assert is_increasing( [11,33,44] )
    assert not is_increasing( [11,44,33] )
    assert not is_increasing( [33,11,44] )

def test_is_standard_young_tableau():
    assert is_standard_young_tableau( [[1,3],[2,4]] )
    assert not is_standard_young_tableau( [[1,3],[2,5]] )
    assert not is_standard_young_tableau( [[3,1],[2,4]] )
    assert not is_standard_young_tableau( [[2,3],[1,4]] )

def test_gl_invariants():
    assert [lc.LinearCombination.from_str("12 - 21",words.ShuffleWord)] == list(gl_invariants(2,2))
    assert [lc.LinearCombination.from_str("1122+2211-1221-2112",words.ShuffleWord), lc.LinearCombination.from_str("1212+2121-1221-2112",words.ShuffleWord)] == list(gl_invariants(2,4))
    for weight in range(1,5):
      assert len( list(gl_invariants(2, weight*2)) ) ==\
              sympy.functions.combinatorial.factorials.binomial( 2 * weight, weight ) / (weight+1) # the Catalan number for weight


def test_d_catalan():
    assert [1, 2, 5, 14, 42, 132, 429, 1430, 4862] == [d_catalan(2,n) for n in range(1,10)]
    assert [1, 5, 42, 462, 6006, 87516, 1385670, 23371634, 414315330] == [d_catalan(3,n) for n in range(1,10)]
    assert [1, 14, 462, 24024, 1662804, 140229804, 13672405890, 1489877926680, 177295473274920] == [d_catalan(4,n) for n in range(1,10)]

def test_dyck_words():
    for dim in [2,3,4]:
        assert [ d_catalan(dim,n) for n in range(1,7-dim) ] == [ len(list(dyck_words(dim,n))) for n in range(1,7-dim) ]

@pytest.mark.slow
def test_permutation_invariants():
    assert [1, 2, 4, 8, 16, 32, 64] == [ len(list(permutation_invariants(2,n))) for n in range(1,8) ]
    assert [1, 2, 5, 14, 41, 122, 365] == [ len(list(permutation_invariants(3,n))) for n in range(1,8) ] # https://oeis.org/A124302
    assert [1, 2, 5, 15, 51, 187, 715] == [ len(list(permutation_invariants(4,n))) for n in range(1,8) ] # https://oeis.org/A007581 
    assert [1, 2, 5, 15, 52, 203] == [ len(list(permutation_invariants(6,n))) for n in range(1,7) ] # https://oeis.org/A056273
                # THIS IS _NOT_ https://oeis.org/A148092, e.g.f.: exp( x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + x^6/720 )

def test_restricted_permutation_invariants():
    for dim in range(1,5):
        for level in range(1,6):
            # since each orbit has cardinality equal to dim:
            assert dim**level == len(list(restricted_permutation_invariants(dim,level))) * dim

def test_equal_size_partitions():
    assert 1 == len( list(equal_size_partitions(range(3),3) ) )
    assert 0 == len( list(equal_size_partitions(range(3),2) ) )
    assert 1 == len( list(equal_size_partitions(range(3),1) ) )

    assert 1 == len( list(equal_size_partitions(range(4),4) ) )
    assert 0 == len( list(equal_size_partitions(range(4),3) ) )
    assert 3 == len( list(equal_size_partitions(range(4),2) ) )
    assert 1 == len( list(equal_size_partitions(range(4),1) ) )

def test_signature_of_path():
    X = np.array([[0., 0.],
                  [1., 0.],
                  [1., 1.],
                  [0., 1.]])
    s = signature_of_path(X, 3)
    a = list(gl_invariants(2,2))[0] # On this level there is only one invariant: the area.

    # Invariants are ShuffleWords. Need to transform into ConcatenationWord first .. This is suboptimal.
    assert 2. == lc.LinearCombination.inner_product( a.apply_linear_function( words.shuffle_word_to_concatenation_word ), s )


#def test_so_invariants():
#    upto = 4
#    assert [ sympy.functions.combinatorial.factorials.binomial( 2*w, w ) for w in range(1,upto)] \
#            == [ 2 * sympy.functions.combinatorial.factorials.binomial( 2*w - 1, w - 1) for w in range(1,upto)] \
#            == list(map(lambda w: words.rank( list(map(lambda inv: words.word_lc_as_vector(inv,2,2*w,words.ShuffleWord), so_invariants(2,2*w)))) , range(1,upto) ))
#    
#    dim = 2
#    upto = 20
#    for level in range(1,upto):
#        print( words.rank( words.word_lcs_as_vectors( so_invariants(dim,level), dim,level,words.ShuffleWord, from_level=level) ) )
#        # [0 2 0 6 0 20 0 70 0 252 0 ]
#        # seems to be https://oeis.org/A000984 TODO
#        
#    #dim = 3
#    #upto = 12
#    #for level in range(1,upto):
#    #    print( lc.rank( lc.word_lcs_as_vectors( so_invariants(dim,level), dim,level,words.ShuffleWord,from_level=level) ) )
#    #    # seems to be https://oeis.org/A005043
#    #    # [0, 1, 1, 3, 6, 15, 36, 91, 232, 603]
