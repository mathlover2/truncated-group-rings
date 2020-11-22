from itertools import islice, product as iprod
from operator import add

## Code for the truncating ideal

def _reducer(x,R,G,sigma_element):


    # Gives $\sum_{g \in G} \left| \left\[ x \right\]_{g} \right|} $.

    def weight(v):
        return sum(abs(x) for x in v.coefficients(False))

    def length(v):
        return v.length()

    def val_ident(v):
        return -abs(v.coefficient(G.identity()))
    
    def weed_out(f, l):
        m = f(min(l, key=f))
        return [v for v in l if f(v) == m]
    
    n = length(x)
    coeffs = sorted(x.coefficients())
    n1 = len(coeffs)
        
    # A simple case of x being already reduced.
    if 2*n1 < n:
        return x

    else:
        coeffs = set(coeffs) | {0}
        candidate_reductions = [x - y*sigma_element for y in coeffs]
        for weeder_function in [length,  weight, val_ident]:
            candidate_reductions = weed_out(weeder_function, candidate_reductions)
            if len(candidate_reductions) == 1:
                return candidate_reductions[0]
        return x - x.coefficient(G.identity())*sigma_element

        
    
def TruncatingIdeal(group_ring):
    sigma_element = sum(group_ring.basis())
    I = group_ring.ideal([sigma_element])
    I.sigma_element = sigma_element
    def reduce(self, x):
        return _reducer(group_ring(x), self.base_ring(),
                        self.ring().group(), self.sigma_element )
    I.reduce = reduce.__get__(I)
    return I

def _reducer_2(x, R, G, sigma_element):
    ident = G.identity()
    xx = x - x.coefficient(ident) * sigma_element
    return xx

def TruncatingIdeal_quick(group_ring):
    sigma_element = sum(group_ring.basis())
    I = group_ring.ideal([sigma_element])
    I.sigma_element = sigma_element
    def reduce(self, x):
        return _reducer_2(group_ring(x), self.base_ring(),
                        self.ring().group(), self.sigma_element )
    I.reduce = reduce.__get__(I)
    return I

## Code for obtaining the list of idempotent elements in the 
## group ring

def idempotents_in_AG(G, R, method='sagemath'):
    if method == 'sagemath':
        RG = G.algebra(R)
        ides = RG.central_orthogonal_idempotents()
        return ides, RG
    elif method == 'gap':
        # The group as a permutation group.
        # Required for the character table calculation
        GG = G.permutation_group()
    
        # Dict to convert elements of GG to elements of G.
        conversion_dict = dict()
        for g, gg in zip(G,GG):
            conversion_dict[gg] = g
    
        RG = G.algebra(R)
        ides = []
        for char in GG.character_table():
            acc = 0
            for i, cc in enumerate(GG.conjugacy_classes_representatives()):
                acc += R(char[i])*RG(conversion_dict[cc])
            acc /= G.order()
            ides.append(acc)
        return ides, RG
    else:
        raise ValueError("Cannot use method \'{0}\'".format(method))

# Example use:
# G = AbelianGroup([3, 3])
# ZZG = G.algebra(ZZ)
# I = TruncatingIdeal(ZZG)      
# ZZGt = QuotientRing(ZZG, I)
# l, GG = idempotents_in_AG(G, R)

## Code for searching for non trivial units
        
def _unit_mult_sums(entries, units, first_one=True):
    n = len(entries)
    
    if first_one:
        signs = ((1,)+combo for combo in iprod(units, repeat=n-1))
    else:
        signs = iprod(units, repeat=n)
        
    for combo in signs:
        yield sum(x*y for (x,y) in zip(combo, entries))

def _is_non_trivial_integral(el):
    l = el.coefficients()
    if len(l) == 1:
        return False
    for coeff in el.coefficients():
        if coeff not in ZZ:
            return False
    return True

def units_of_tgr(G_or_idemps, R=QQ):
    if G_or_idemps in Groups:
        l, RG = idempotents_in_AG(G_or_idemps, R)
    else:
        l, RG = G_or_idemps
        R = RG.base_ring()

    II = TruncatingIdeal_quick(RG)
    lt = list(map(II.reduce,l[1:]))
    units = R.roots_of_unity() if R is not QQ else [1,-1]
    
    for i in _unit_mult_sums(lt, units):
        yield II.reduce(i)

def non_trivial_integral_units_of_tgr(G_or_idemps, R=QQ):
    for i in units_of_tgr(G_or_idemps, R):
        if _is_non_trivial_integral(i):
            yield i

## Code for calculating change-of-basis matrices.            

def RG_G_to_idemp_cobm(G, R):
    """Gives a change of basis matrix for $RG$ for changing from the basis formed by elements
    of $G$ to the basis of idempotents of $RG$"""
    n = order(G)
    Gl = list(G)
    m = matrix(R, n, n)
    ides, RG = idempotents_in_AG(G, R)
    for i in range(n):
        for j in range(n):
            m[i, j] = ((ides[j].coefficient(Gl[i]))^(-1))/n
    return m

def RGt_G_to_idemp_cobm(G, R):
    """Gives a change of basis matrix for $RG_t$ for changing from the basis formed by elements
    of $G - \left\{e\right\} to the basis of idempotents of $RG$ other than $e_\chi_0$"""
    return RG_G_to_idemp_cobm(G, R).submatrix(1,1)

