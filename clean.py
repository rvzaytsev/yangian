from itertools import chain
from collections import defaultdict
from tabulate import tabulate

EMPTY = tuple()

# These are technical parameters which set up the context
FAST = True
ALGEBRA = 'free'  # C^N or free


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


def commute(w, w_tilde, index):
    """
    Input:
    w, w_tilde are tuples of free algebra elements (which we consider as strings)
    index is from 1 to 4 - the index of the formula

    Output:
    phi, psi - dictionaries (z1, z2) -> int
    where z1 and z2 are tuples of algebra elements
    """

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
        else:
            raise Exception(f"Algebra {ALGEBRA} not implemented. Choose C^N or free")
        # COMMUTATIVE CASE
        # a, b = w[0], w_tilde[0]
        # ab = ({a,b},)
        # ba = ({a,b},)
        # if index == 1 or index == 2:
        #     return [(ba, EMPTY, -1), (EMPTY, ab, 1)], []
        # elif index == 3 or index == 4:
        #     return [(ab, EMPTY, 1), (EMPTY, ba, -1)], []

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


def free_algebra_commute(n, m, index):
    assert ALGEBRA == 'free'
    w = tuple(map(str, range(1, n + 1)))
    w_tilde = tuple(map(str, range(n + 1, n + 1 + m)))
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
    assert ALGEBRA == 'free'
    index = int(input('Enter which formula to compute (1-4):'))
    n = int(input('Enter the number of variables for w:'))
    m = int(input('Enter the number of variables for w_tilde:'))

    w = tuple(map(str, range(1, n + 1)))
    w_tilde = tuple(map(str, range(n + 1, n + 1 + m)))
    print(f'w={w}')
    print(f'w_tilde={w_tilde}.')
    pretty_print_formulas(w, w_tilde, index)
