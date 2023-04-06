from itertools import product, permutations
from clean import commute, free_algebra_commute, substitute
from clean_copy import commute as commute_copy, free_algebra_commute as free_algebra_commute_copy
from misha import Formula_I, Formula_II, Formula_III, Formula_IV
from clean import ALGEBRA as CLEAN_ALGEBRA, canonical_words, pretty_print, project
from clean_copy import ALGEBRA as COPY_ALGEBRA
from collections import defaultdict
from double_bracket import compute_first_psi_term_explicitly
from itertools import chain
def normal_form(formula):
    formula_list = []
    for z1, z2, val in formula:
        formula_list.append([';'.join(z1), ';'.join(z2), val])
    formula_list.sort(key=lambda item: (-item[2], item[0], item[1]))
    return formula_list
def list_form(dict_form):
    formula = []
    for (z1, z2), val in dict_form.items():
        if val != 0:
            formula.append((z1, z2, val))
    return formula

def sort_formulas(phi, psi):
    return sorted(phi), sorted(psi)
def compare_my_and_misha(m, n):
    assert CLEAN_ALGEBRA == 'C^N', "Only run this test for C^N algebra"
    mishas = [Formula_I, Formula_II, Formula_III, Formula_IV]
    for w in product(('1', '2'), repeat=m):
        for w_tilde in product(('1', '2'), repeat=n):
            for i in range(4):
                phi_my, psi_my = commute(w, w_tilde, i+1)
                phi_my.sort()
                psi_my.sort()
                phi_misha, psi_misha = mishas[i](w, w_tilde)
                phi_misha, psi_misha = list_form(phi_misha), list_form(psi_misha)
                phi_misha.sort()
                psi_misha.sort()
                assert(phi_my == phi_misha and psi_my == psi_misha)

def compare_word_division(m, n):
    for w in product(('1', '2'), repeat=m):
        for w_tilde in product(('1', '2'), repeat=n):
            for i in range(1,4):
                phi, psi = commute(w, w_tilde, i+1)
                phi.sort()
                psi.sort()
                phi_copy, psi_copy = commute_copy(w, w_tilde, i+1)
                phi_copy.sort()
                psi_copy.sort()
                try:
                    assert phi == phi_copy and psi == psi_copy
                except:
                    print(w, w_tilde)
                    #exit(0)

def compare_word_division_free_algebra(n, m):
    for index in range(1,5):
        assert sort_formulas(*free_algebra_commute(n, m, index)) ==\
               sort_formulas(*free_algebra_commute_copy(n, m, index))


def remove_letter(formula, letter):
    dict_form = defaultdict(int)
    for w1, w2, val in formula:
        if letter in w1 or letter in w2:
            w1_new = tuple(c for c in w1 if c != letter)
            w2_new = tuple(c for c in w2 if c != letter)
            dict_form[(w1_new, w2_new)] += val
    result = []
    for (w1, w2), val in dict_form.items():
        if val != 0:
            result.append((w1, w2, val))
    return result
def free_algebra_remove_letter_test(n, m):
    w, w_tilde = canonical_words(n, m)
    for index in range(4, 5):
        full_phi, full_psi = free_algebra_commute(n, m, index)
        for i, letter in enumerate(w):
            phi_reduced = remove_letter(full_phi, letter)
            psi_reduced = remove_letter(full_psi, letter)
            phi, psi = commute(w[:i] + w[i+1:], w_tilde, index)
            pretty_print(phi_reduced)
            pretty_print(phi)
            assert sort_formulas(phi_reduced, psi_reduced) == sort_formulas(phi, psi)


def sort_letters(w):
    return tuple(''.join(sorted(c)) for c in w)
def free_to_commutative(formula):
    dict_form = defaultdict(int)
    for w1, w2, val in formula:
        dict_form[(sort_letters(w1), sort_letters(w2))] += val
    result = []
    for (w1, w2), val in dict_form.items():
        if val != 0:
            result.append((w1, w2, val))
    return result
# actually this is more of a hypothesis
def free_algebra_reduce_to_commutative_test(n, m):
    assert CLEAN_ALGEBRA == 'free' and COPY_ALGEBRA == 'commutative'
    for i in range(4,5):
        phi_free, psi_free = free_algebra_commute(n, m, i)
        phi_com, psi_com = free_algebra_commute_copy(n, m, i)
        phi_reduced, psi_reduced = free_to_commutative(phi_free), free_to_commutative(psi_free)
        if sort_formulas(phi_reduced, psi_reduced) != sort_formulas(phi_com, psi_com):
            pretty_print(phi_free)
            pretty_print(phi_reduced)
            pretty_print(phi_com)
            return

def reduce_free_var_to_c(v):
    if len(set(v)) <= 1:
        return '' if len(v) == 0 else v[0]
    return 0
def reduce_word_to_c(w):
    c_word = []
    for v in w:
        c = reduce_free_var_to_c(v)
        if c == 0:
            return 0
        c_word.append(c)
    return tuple(c_word)
def free_reduce_to_c(formula, w, w_tilde):
    formula = substitute(w, w_tilde, formula)
    dict_form = defaultdict(int)
    for w1, w2, val in formula:
        z1 = reduce_word_to_c(w1)
        z2 = reduce_word_to_c(w2)
        if z1 == 0 or z2 == 0:
            continue
        dict_form[(z1, z2)] += val
    result = []
    for (w1, w2), val in dict_form.items():
        if val != 0:
            result.append((w1, w2, val))
    return result
def free_algebra_reduce_to_c_test(w, w_tilde):
    assert CLEAN_ALGEBRA == 'free' and COPY_ALGEBRA == 'C^N'
    for i in range(1, 5):
        phi_free, psi_free = free_algebra_commute(len(w), len(w_tilde), i)
        phi_c, psi_c = commute_copy(w, w_tilde, i)
        phi_reduced, psi_reduced = free_reduce_to_c(phi_free, w, w_tilde),\
            free_reduce_to_c(psi_free, w, w_tilde)
        if sort_formulas(phi_reduced, psi_reduced) != sort_formulas(phi_c, psi_c):
            print('NOT EQUAL!!')
            pretty_print(phi_free)
            pretty_print(phi_reduced)
            pretty_print(phi_c)
            return

def test_explicit_psi(n, m):
    phi, psi = free_algebra_commute(n, m, 1)
    psi_explicit = compute_first_psi_term_explicitly(n, m)
    psi_explicit.sort(key=lambda item: (-item[2], item[0], item[1]))
    proj = []
    for z1, z2, val in psi:
        if len(z1) + len(z2) + 2 == n + m:
            proj.append((z1, z2, val))
    proj.sort(key=lambda item: (-item[2], item[0], item[1]))

    assert proj == psi_explicit

