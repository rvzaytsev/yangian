from itertools import chain
from collections import defaultdict
from tabulate import tabulate
import string

EMPTY = tuple()

# These are technical parameters which set up the context
FAST = True
ALGEBRA = 'free'#'commutative'  # C^N or free or commutative
CACHE = {} # for free algebra it is very easy: we store the map (n, m, i) -> formula_i([1,2,...,n], [n+1,...,n+1+m])
ALPHABET = string.printable
# arithmetic


# negates the coefficients
def negate(formulas):
    phi, psi = formulas
    phi_new = []
    psi_new = []
    for w1, w2, val in phi:
        phi_new.append((w1, w2, -val))
    for w1, w2, val in psi:
        psi_new.append((w1, w2, -val))
    return phi_new, psi_new


# adds together all the coefficients, assuming they have the same form
# makes sure there are no zeros after addition
def add(*args):
    dict_form = defaultdict(int)
    for w1, w2, val in chain.from_iterable(args):
        dict_form[(w1, w2)] += val
    result = []
    for (w1, w2), val in dict_form.items():
        if val != 0:
            result.append((w1, w2, val))
    return result


def transpose(formula):
    new_formula = []
    for w1, w2, val in formula:
        new_formula.append((w2, w1, val))
    return new_formula


def multiply(formula, scalar):
    new_formula = []
    for w1, w2, val in formula:
        new_formula.append((w1, w2, val * scalar))
    return new_formula


def append(formula, word):
    new_formula = []
    for z1, z2, val in formula:
        new_formula.append((z1, z2 + word, val))
    return new_formula


def prepend(formula, word):
    new_formula = []
    for z1, z2, val in formula:
        new_formula.append((word + z1, z2, val))
    return new_formula


# rewrites e_ij e_kl as e_kl e_ij + [e_ij, e_kl], and the commutator is expressed with a given index
# def swap(formula, index):

# TODO: adding this way is not efficient. It is not on-line, though it could be
def il_kj_to_1(formula):
    phis = [transpose(formula)]
    psis = []
    for z1, z2, val in formula:
        psi, phi = commute(z1, z2, 4)
        phis.append(multiply(phi, val))
        psis.append(multiply(psi, val))
    phi_new = add(*phis)
    psi_new = add(*psis)
    return phi_new, psi_new


def kl_ij_to_1(formula):
    phis = []
    psis = [transpose(formula)]
    for z1, z2, val in formula:
        phi, psi = commute(z1, z2, 4)
        phis.append(multiply(phi, val))
        psis.append(multiply(psi, val))
    phi_new = add(*phis)
    psi_new = add(*psis)
    return phi_new, psi_new


def il_kj_to_2(formula):
    phis = [transpose(formula)]
    psis = []
    for z1, z2, val in formula:
        psi, phi = commute(z1, z2, 2)
        phis.append(multiply(phi, val))
        psis.append(multiply(psi, val))
    phi_new = add(*phis)
    psi_new = add(*psis)
    return phi_new, psi_new


def ij_kl_to_2(formula):
    phis = []
    psis = [transpose(formula)]
    for z1, z2, val in formula:
        phi, psi = commute(z1, z2, 2)
        phis.append(multiply(phi, val))
        psis.append(multiply(psi, val))
    phi_new = add(*phis)
    psi_new = add(*psis)
    return phi_new, psi_new


def kl_ij_to_3(formula):
    phis = []
    psis = [transpose(formula)]
    for z1, z2, val in formula:
        phi, psi = commute(z1, z2, 2)
        phis.append(multiply(phi, val))
        psis.append(multiply(psi, val))
    phi_new = add(*phis)
    psi_new = add(*psis)
    return phi_new, psi_new


def ij_kl_to_4(formula):
    phis = []
    psis = [transpose(formula)]
    for z1, z2, val in formula:
        phi, psi = commute(z1, z2, 4)
        phis.append(multiply(phi, val))
        psis.append(multiply(psi, val))
    phi_new = add(*phis)
    psi_new = add(*psis)
    return phi_new, psi_new


# convert a sum of formulas 3 and 4 to a particular formula, e.g. 1
def convert(phi3, psi3, phi4, psi4, index):
    if index == 1:
        phi_from_phi3, psi_from_phi3 = il_kj_to_1(phi3)
        phi_from_phi4, psi_from_phi4 = il_kj_to_1(phi4)
        phi_from_psi3, psi_from_psi3 = [], psi3
        phi_from_psi4, psi_from_psi4 = kl_ij_to_1(psi4)
    elif index == 2:
        phi_from_phi3, psi_from_phi3 = il_kj_to_2(phi3)
        phi_from_phi4, psi_from_phi4 = il_kj_to_2(phi4)
        phi_from_psi3, psi_from_psi3 = ij_kl_to_2(psi3)
        phi_from_psi4, psi_from_psi4 = [], psi4
    elif index == 3:
        phi_from_phi3, psi_from_phi3 = phi3, []
        phi_from_phi4, psi_from_phi4 = phi4, []
        phi_from_psi3, psi_from_psi3 = [], psi3
        phi_from_psi4, psi_from_psi4 = kl_ij_to_3(psi4)
    elif index == 4:
        phi_from_phi3, psi_from_phi3 = phi3, []
        phi_from_phi4, psi_from_phi4 = phi4, []
        phi_from_psi3, psi_from_psi3 = ij_kl_to_4(psi3)
        phi_from_psi4, psi_from_psi4 = [], psi4
    phi = add(phi_from_phi3, phi_from_phi4, phi_from_psi3, phi_from_psi4)
    psi = add(psi_from_phi3, psi_from_phi4, psi_from_psi3, psi_from_psi4)
    return phi, psi

def apply_mapping(z, mapping):
    word = []
    for x in z:
        y = ""
        for c in x:
            y += mapping[c]
        if ALGEBRA == 'commutative':
            word.append(''.join(sorted(y)))
        else:
            word.append(y)
    return tuple(word)
def substitute(w, w_tilde, formula):
    merged = w + w_tilde
    mapping = {}
    for u in range(len(merged)):
        mapping[ALPHABET[u]] = merged[u]
    substituted_formula = []
    for z1, z2, val in formula:
        w1 = apply_mapping(z1, mapping)
        w2 = apply_mapping(z2, mapping)
        substituted_formula.append((w1, w2, val))
    return substituted_formula

def commute(w, w_tilde, index):
    """
    Input:
    w, w_tilde are tuples of free algebra elements (which we consider as strings)
    index is from 1 to 4 - the index of the formula

    Output:
    phi, psi - dictionaries (z1, z2) -> int
    where z1 and z2 are tuples of algebra elements
    """
    # do lazy caching
    if (ALGEBRA == 'free' or ALGEBRA == 'commutative') and ((len(w), len(w_tilde), index) not in CACHE):
        n, m = len(w), len(w_tilde)
        w_canonical, w_tilde_canonical = canonical_words(n, m)
        if w != w_canonical or w_tilde != w_tilde_canonical:
            CACHE[(n, m, index)] = commute(w_canonical, w_tilde_canonical, index)

    # check if we have an answer in cache
    if (ALGEBRA == 'free' or ALGEBRA == 'commutative') and (len(w), len(w_tilde), index) in CACHE:
        phi, psi = CACHE[(len(w), len(w_tilde), index)]
        return substitute(w, w_tilde, phi), substitute(w, w_tilde, psi)
    # do honest computation, assuming w and w_tilde are just increasing nums

    # trivial commutator
    if len(w) == 0 or len(w_tilde) == 0:
        return [], []

    # base case
    if len(w) + len(w_tilde) == 2:
        if ALGEBRA == 'C^N':
            a, b = w[0], w_tilde[0]
            if a != b:
                return [], []
            if index == 1 or index == 2:
                return [((a,), EMPTY, -1), (EMPTY, (a,), 1)], []
            elif index == 3 or index == 4:
                return [((a,), EMPTY, 1), (EMPTY, (a,), -1)], []

        ### THE BLOCK BELOW IS FOR THE FREE ALGEBRA CASE ###
        elif ALGEBRA == 'free':
            a, b = w[0], w_tilde[0]
            ab = (a + b,)
            ba = (b + a,)
            if index == 1 or index == 2:
                return [(ba, EMPTY, -1), (EMPTY, ab, 1)], []
            elif index == 3 or index == 4:
                return [(ab, EMPTY, 1), (EMPTY, ba, -1)], []
        elif ALGEBRA == 'commutative':
            a, b = w[0], w_tilde[0]
            ab = (''.join(sorted(a + b)),)
            ba = (''.join(sorted(b + a)),)
            if index == 1 or index == 2:
                return [(ba, EMPTY, -1), (EMPTY, ab, 1)], []
            elif index == 3 or index == 4:
                return [(ab, EMPTY, 1), (EMPTY, ba, -1)], []
        else:
            raise Exception(f"Algebra {ALGEBRA} not implemented. Choose C^N or free or commutative")
        # COMMUTATIVE CASE

    # due to symmetry we can exclude the case that w has length 1
    if len(w) < len(w_tilde):
        return negate(commute(w_tilde, w, 5 - index))

    # now we can assume that w has length at least 2, so we can split it
    length = len(w)
    if FAST:
        w_begin = w[:length // 2]
        w_end = w[length // 2:]
    else:
        w_begin = w[:1]
        w_end = w[1:]
    # according to the lemma we have to compute the 3rd and 4th formulas
    phi4, psi4 = commute(w_begin, w_tilde, 4)
    phi4_new = append(phi4, w_end)
    psi4_new = append(psi4, w_end)

    phi3, psi3 = commute(w_end, w_tilde, 3)
    phi3_new = prepend(phi3, w_begin)
    psi3_new = prepend(psi3, w_begin)

    return convert(phi3_new, psi3_new, phi4_new, psi4_new, index)
def canonical_words(n, m):
    w = tuple(c for c in ALPHABET[:n])
    w_tilde = tuple(c for c in ALPHABET[n:n+m])
    return w, w_tilde
def free_algebra_commute(n, m, index):
    assert ALGEBRA == 'free' or ALGEBRA == 'commutative'
    w, w_tilde = canonical_words(n, m)
    return commute(w, w_tilde, index)


def pretty_print(formula):
    formula_list = []
    # formula.sort(key=lambda item: (-len(item[0]) - len(item[1]), -item[2]))
    for z1, z2, val in formula:
        formula_list.append([';'.join(z1), ';'.join(z2), val])
    formula_list.sort(key=lambda item: (-item[2], item[0], item[1]))
    print(tabulate(formula_list, headers=['w', 'w_tilde', 'coef']))


def pretty_print_formulas(w, w_tilde, index):
    phi, psi = commute(w, w_tilde, index)
    print("φ part")
    pretty_print(phi)
    print('\nψ part')
    pretty_print(psi)


def free_algebra_io():
    assert ALGEBRA == 'free' or ALGEBRA == 'commutative'
    print(f'Computing in {ALGEBRA} algebra')
    index = int(input('Enter which formula to compute (1-4):'))
    n = int(input('Enter the number of variables for w:'))
    m = int(input('Enter the number of variables for w_tilde:'))
    w, w_tilde = canonical_words(n, m)
    print(f'w={w}')
    print(f'w_tilde={w_tilde}.')
    pretty_print_formulas(w, w_tilde, index)

def c_io():
    assert ALGEBRA == 'C^N'
    index = int(input('Enter which formula to compute (1-4):'))
    print('Use the following format: numbers (algebra vars) separated by semi-colons ("train")')
    print('For example: 1;2;1')
    w = tuple(input("Enter w: ").split(';'))
    w_tilde = tuple(input("Enter w_tilde: ").split(';'))
    print(f'w={w}', f'w_tilde={w_tilde}')
    print(f"index={1}")
    pretty_print_formulas(w, w_tilde, 1)

def project(formula, r, s):
    projection = []
    for z1, z2, val in formula:
        if len(z1) == r and len(z2) == s:
            projection.append((z1, z2, val))
    return projection

def project_total_degree(formula, degree):
    projection = []
    for z1, z2, val in formula:
        if len(z1) + len(z2) == degree:
            projection.append((z1, z2, val))
    return projection

def project_total_degree_io():
    assert ALGEBRA == 'free' or ALGEBRA == 'commutative'
    print(f'Computing in {ALGEBRA} algebra')
    index = int(input('Enter which formula to compute (1-4):'))
    n = int(input('Enter the number of variables for w:'))
    m = int(input('Enter the number of variables for w_tilde:'))
    d = int(input('Enter which degree to project to:'))
    w, w_tilde = canonical_words(n, m)
    print(f'w={w}')
    print(f'w_tilde={w_tilde}.')
    phi, psi = free_algebra_commute(n, m, index)
    phi_proj, psi_proj = project_total_degree(phi, d), project_total_degree(psi, d)
    print("φ part of degree", d)
    pretty_print(phi_proj)
    print('\nψ part of degree', d)
    pretty_print(psi_proj)


def free_algebra_projection(n, m, r=1, s=1, index=4):
    assert ALGEBRA == 'commutative'
    w, w_tilde = canonical_words(n, m)
    print(f'w={w}', f'w_tilde={w_tilde}')
    phi, psi = commute(w, w_tilde, index)
    if (n + m - r - s) % 2:
        formula, zero = phi, psi
    else:
        formula, zero = psi, phi
    assert len(project(zero, r, s)) == 0
    projection = project(formula, r, s)
    print(len(projection))
    pretty_print(projection)

def nonzero_projections(n, m, index=4):
    phi, psi = free_algebra_commute(n, m, index=4)
    phi_proj = []
    psi_proj = []
    phi_degrees = set()
    psi_degrees = set()
    for r in range(n + m):
        for s in range(n + m - r):
            if len(project(phi, r, s)) > 0:
                phi_proj.append((r, s))
                phi_degrees.add(r + s)
            if len(project(psi, r, s)) > 0:
                psi_proj.append((r, s))
                psi_degrees.add(r+ s)
    return phi_proj, psi_proj, phi_degrees, psi_degrees

def extract_subsequence(w, sequence):
    subsequence = []
    for letter in w:
        for c in letter:
            if c in sequence:
                subsequence.append(c)
    return subsequence

def get_permutation(word, permuted):
    permutation = [0]*len(word)
    for i, c in enumerate(permuted):
        permutation[word.index(c)] = i+1
    # for i, c in enumerate(permuted):
    #     permutation[i] = word.index(c) + 1
    return permutation
def get_word_permutations(n, m, index=4):
    w, w_tilde = canonical_words(n, m)
    phi, psi = commute(w, w_tilde, index)
    seconds_phi, seconds_psi = [], []
    for z1, z2, val in phi:
        second = extract_subsequence(z1+z2, w_tilde)
        second = tuple(map(str, get_permutation(w_tilde, second)))
        seconds_phi.append((z1, z2, val, second))
    for z1, z2, val in psi:
        second = extract_subsequence(z1 + z2, w_tilde)
        second = tuple(map(str, get_permutation(w_tilde, second)))
        seconds_psi.append((z1, z2, val, second))
    return seconds_phi, seconds_psi

def print_formulas_with_permutation(n, m, index=4):
    assert ALGEBRA == 'free' or ALGEBRA == 'commutative'
    print(canonical_words(n, m))
    phi, psi = get_word_permutations(n, m, index)
    print("phi part")
    pretty_print_with_permutation(phi)
    print("psi part")
    pretty_print_with_permutation(psi)

def pretty_print_with_permutation(formula):
    formula_list = []
    perm_pos = []
    perm_neg = []
    for z1, z2, val, perm in formula:
        if val > 0:
            perm_pos.append(perm)
        else:
            perm_neg.append(perm)

        formula_list.append([';'.join(z1), ';'.join(z2), val, ' '.join(perm)])
    formula_list.sort(key=lambda item: (-item[2], item[0], item[1]))
    print(tabulate(formula_list, headers=['w', 'w_tilde', 'coef', 'permutation']))
    assert sorted(perm_pos) == sorted(perm_neg)
def print_permutations(n, m):
    seconds_phi, seconds_psi, total = get_word_permutations(n, m)
    print(total)
    print("Phi part")
    for w, val in seconds_phi:
        if val == 1:
            print("+", ' '.join(w))
        else:
            print("-", ' '.join(w))
    print("Psi part")
    for w, val in seconds_psi:
        if val == 1:
            print("+", ' '.join(w))
        else:
            print("-", ' '.join(w))