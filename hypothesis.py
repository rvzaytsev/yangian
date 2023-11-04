from clean import free_algebra_commute, ALGEBRA, canonical_words, commute, add, negate
from itertools import permutations, product, combinations
import double_bracket
import clean
import numpy as np
import sympy
from collections import defaultdict
from tabulate import tabulate
import tqdm

def print_coefficients_free_algebra(n, m):
    print(f"words are of length {n} and {m}")
    for i in range(4,5):
        print("Computing formula #"+str(i))
        phi, psi = free_algebra_commute(n, m, i)
        print(f"Total number of coefficients: {len(phi) + len(psi)}")
        phi_coefs = [item[2] for item in phi]
        psi_coefs = [item[2] for item in psi]
        print(f"Coefficients are in {set(phi_coefs) | set(psi_coefs)}")

def check_coefs(n, m):
    for u in range(1, n+1):
        for v in range(1, m+1):
            for i in range(1, 5):
                phi, psi = free_algebra_commute(u, v, i)
                phi_coefs = [item[2] for item in phi]
                psi_coefs = [item[2] for item in psi]
                assert (set(phi_coefs) | set(psi_coefs)).issubset({-1,1})


def get_all_formulas_free(total_length):
    alphabet, _ = canonical_words(total_length, 0)
    formulas = {}
    for x in permutations(alphabet):
        for n in range(total_length+1):
            w, w_tilde = x[:n], x[n:]
            formulas[(w, w_tilde)] = commute(w, w_tilde, 4)
    return formulas

def get_all_formulas_fixed_lengths(n, m):
    alphabet, _ = canonical_words(n+m, 0)
    formulas = {}
    for x in permutations(alphabet):
        w, w_tilde = x[:n], x[n:]
        formulas[(w, w_tilde)] = commute(w, w_tilde, 4)
    return formulas

def partition(number):
    answer = set()
    answer.add((number, ))
    for x in range(1, number):
        for y in partition(number - x):
            answer.add(tuple(sorted((x, ) + y)))
    return answer

def slice_by_indices(indices, array):
    result = []
    idx = np.cumsum((0,) + indices)
    for i in range(len(indices)):
        result.append(''.join(array[idx[i]:idx[i+1]]))
    return result

def cyclic_permutations(iter):
    for i in range(len(iter)):
        yield iter[i:] + iter[:i]
def get_all_tabloids(total_length, cyclic=False, with_permutations=True):
    assert ALGEBRA == 'free' or ALGEBRA == 'commutative'
    alphabet, _ = canonical_words(total_length, 0)
    tabloids = []
    if cyclic:
        perm = cyclic_permutations
    else:
        perm = permutations
    for letters in perm(alphabet):
        for indices in partition(total_length):
            for index in set(permutations(indices)):
                merged = slice_by_indices(index, letters)
                if ALGEBRA == 'commutative':
                    merged = [''.join(sorted(letter)) for letter in merged]
                for k in range(1, len(merged)):
                    w, w_tilde = tuple(merged[:k]), tuple(merged[k:])
                    tabloids.append((w, w_tilde))
        if not with_permutations:
            break
    return tabloids

def get_all_formulas_gluing(total_length):
    formulas = {}
    for w, w_tilde in get_all_tabloids(total_length):
        formulas[(w, w_tilde)] = commute(w, w_tilde, index=4)
    return formulas

def find_good_combinations(n, m, num_coeffs):
    minimal = []
    min_coeffs = 100
    w, w_tilde = canonical_words(n, m)
    formulas = get_all_formulas_gluing(n+m)
    phi1, psi1 = formulas[(w, w_tilde)]
    for combo in combinations(formulas.items(), num_coeffs-1):
        for signs in product((-1, 1), repeat=num_coeffs-1):
            phi = phi1.copy()
            psi = psi1.copy()
            w1s = [w]
            w2s = [w_tilde]
            charsigns = ['+']
            for j in range(num_coeffs - 1):
                (v1, v2), (phi2, psi2) = combo[j]
                if w != v1 or v2 != w_tilde:
                    sign = signs[j]
                    w1s.append(v1)
                    w2s.append(v2)
                    charsigns.append('+' if sign == 1 else '-')
                    if sign == -1:
                        phi2, psi2 = negate((phi2, psi2))
                    phi, psi = add(phi, phi2), add(psi, psi2)
            length = len(phi) + len(psi)
            if length < min_coeffs:
                minimal = set()
                min_coeffs = length
            if length <= min_coeffs:
                minimal.add((tuple(w1s), tuple(w2s), tuple(charsigns)))
    return list(minimal), min_coeffs

def find_linear_combinations(n, m):
    basis2index_phi = {}
    basis2index_psi = {}
    w, w_tilde = canonical_words(n, m)
    formulas = get_all_formulas_fixed_lengths(n,m)#get_all_formulas_gluing(n+m)
    index = 0
    for phi, psi in formulas.values():
        for z1, z2, val in phi:
            if (z1, z2) not in basis2index_phi:
                basis2index_phi[z1, z2] = index
                index += 1
        for z1, z2, val in psi:
            if (z1, z2) not in basis2index_psi:
                basis2index_psi[z1, z2] = index
                index += 1

    #phi, psi = formulas.pop((w, w_tilde))
    # now construct arrays
    length = len(basis2index_phi) + len(basis2index_psi)
    brackets_matrix = np.zeros((length, len(formulas)))
    for i, (ph, ps) in enumerate(formulas.values()):
        for z1, z2, val in ph:
            brackets_matrix[basis2index_phi[z1, z2], i] = val
        for z1, z2, val in ps:
            brackets_matrix[basis2index_psi[z1, z2], i] = val
    idx = list(formulas.keys()).index((w, w_tilde))

    # _, inds = sympy.Matrix(brackets_matrix).rref()
    # expressible_inds = sorted(list(set(list(range(len(formulas)))) - set(inds)))
    # print("expressible words")
    # for i in expressible_inds:
    #     u, v = list(formulas.keys())[i]
    #     print(u, v)
    # express_columns = []
    # #print("Expressible inds", expressible_inds)
    # for idx in expressible_inds:
    #     express_columns.append(brackets_matrix[:, idx])
    # brackets_matrix = np.delete(brackets_matrix, expressible_inds, 1)
    # t=0
    # solution, err, rank, _ = np.linalg.lstsq(brackets_matrix, express_columns[t], rcond=None)
    col = brackets_matrix[:, idx].copy()
    brackets_matrix[:, idx] = 0
    # for i, (z1, z2) in enumerate(formulas.keys()):
    #     if len(z1) + len(z2) == n + m:
    #         brackets_matrix[:, i] = 0 #np.delete(brackets_matrix, idx, 1)
    solution, err, rank, _ = np.linalg.lstsq(brackets_matrix, col, rcond=None)
    assert abs(solution[idx]) < 1e-9
    print("Matrix rank", rank)
    print("Number of tabloids", len(formulas))
    # print(err[0])
    rounded_solution = solution.round().astype('int')
    print(1, list(formulas.keys())[idx])
    for i, s in enumerate(rounded_solution):
        if abs(s) > 1e-3:
            print(-s, list(formulas.keys())[i])
    # print(-1, list(formulas.keys())[expressible_inds[t]])
    # for i, s in enumerate(rounded_solution):
    #     if s != 0:
    #         print(s, list(formulas.keys())[inds[i]])
    l2_error = np.sum((brackets_matrix @ rounded_solution - col) ** 2)
    err = np.sum((brackets_matrix @ solution - col)**2)
    print("Errors:", l2_error, err)
    # l2_error = np.sum((brackets_matrix @ rounded_solution - express_columns[t]) ** 2)
    return rounded_solution, l2_error
    rounded_solution = np.clip(solution, -1, 1).round().astype('int')
    l2_error = np.sum((brackets_matrix @ rounded_solution - formula_array) ** 2)
    terms = []
    for i, val in enumerate(rounded_solution):
        if val != 0:
            terms.append((brackets[i][0], val))
    # print(solution, rounded_solution)
    return terms, l2_error, err

# brute force algo, a linear optimization would be much faster
def find_good_linear_combinations(total_length, num_coeffs):
    formulas = get_all_formulas_free(total_length)
    minimal = []
    min_coeffs = 100
    for n in range(total_length+1):
        if n != total_length // 2 + 1:
            continue
        w1, w2 = canonical_words(n, total_length-n)
        phi1, psi1 = formulas[(w1, w2)]
        for combo in combinations(formulas.items(), num_coeffs-1):
            for signs in product((-1, 1), repeat=num_coeffs-1):
                phi = phi1.copy()
                psi = psi1.copy()
                w1s = [w1]
                w2s = [w2]
                charsigns = []
                for j in range(num_coeffs - 1):
                    (v1, v2), (phi2, psi2) = combo[j]
                    sign = signs[j]
                    w1s.append(v1)
                    w2s.append(v2)
                    charsigns.append('+' if sign == 1 else '-')
                    if w1 != v1 or v2 != w2:
                        if sign == -1:
                            phi2, psi2 = negate((phi2, psi2))
                        phi, psi = add(phi, phi2), add(psi, psi2)
                length = len(phi) + len(psi)
                if length < min_coeffs:
                    minimal = []
                    min_coeffs = length
                if length <= min_coeffs:
                    minimal.append((w1s, w2s, charsigns))
    return minimal, min_coeffs

def first_word_preserves_order(n, m):
    w, w_tilde = canonical_words(n, m)
    phi, psi = commute(w, w_tilde, 4)
    for z1, z2, val in phi + psi:
        first_permutation = []
        merge = ''.join(z1+z2)
        for c in merge:
            if c in w:
                first_permutation.append(c)
        #assert len(first_permutation) == len(w)
        assert tuple(first_permutation) == w

def first_word_no_gluing(n, m):
    w, w_tilde = canonical_words(n, m)
    phi, psi = commute(w, w_tilde, 4)
    for z1, z2, val in phi + psi:
        for x in z1 + z2:
            count = 0
            for c in x:
                if c in w:
                    count += 1
            assert count <= 1
        #assert len(first_permutation) == len(w)

def detailed_sign_equality(n, m):
    w, w_tilde = canonical_words(n, m)
    phi, psi = commute(w, w_tilde, 4)
    plus, minus = 0, 0
    for z1, z2, val in phi:
        if val == 1:
            plus += 1
        elif val == -1:
            minus += 1
        else:
            raise Exception("Non trivial coefficient found!")
    assert plus == minus
    plus, minus = 0, 0
    for z1, z2, val in psi:
        if val == 1:
            plus += 1
        elif val == -1:
            minus += 1
        else:
            raise Exception("Non trivial coefficient found!")
    assert plus == minus

from sympy.combinatorics.permutations import Permutation
def calculate_permutation_sum(n, m):
    assert ALGEBRA == 'commutative'
    w, w_tilde = canonical_words(n, m)
    alphabet = w+w_tilde
    phi, psi = [], []
    for permutation in permutations(range(n+m)):
        sign = Permutation(permutation).signature()
        word = ''.join(alphabet[p] for p in permutation)
        w= word[:n]
        w_tilde = word[n:]
        phi_cur, psi_cur = commute(w, w_tilde, 4)
        if sign == -1:
            phi_cur, psi_cur = negate((phi_cur, psi_cur))
        phi = add(phi, phi_cur)
        psi = add(psi, psi_cur)
    return phi, psi

def test_odd_phi_even_psi(n, m):
    for index in range(1,5):
        phi, psi = free_algebra_commute(n, m, index)
        for w, w_tilde, val in phi:
            assert (n+m - (len(w) + len(w_tilde))) % 2 == 1
        for w, w_tilde, val in psi:
            assert (n+m - (len(w) + len(w_tilde))) % 2 == 0

def compare_explicit_formula_psi_with_bracket(n, m):
    psi = double_bracket.compute_first_psi_term_explicitly(n, m)
    degs = set()
    for z1, z2, val in psi:
        degs.add((len(z1), len(z2)))
    ts = double_bracket.compute_iterations(double_bracket.canonical_tensor(n, m))
    if len(ts) < 2:
        assert psi == []
        return [], []
    t = clean.multiply(ts[1], -1)
    proj = []
    degs_other = set()
    for z1, z2, val in t:
        deg = (len(z1), len(z2))
        if deg in degs:
            proj.append((z1, z2, val))
        else:
            degs_other.add(deg)
    psi.sort(key=lambda item: (-item[2], item[0], item[1]))
    proj.sort(key=lambda item: (-item[2], item[0], item[1]))
    # clean.pretty_print(psi)
    # clean.pretty_print(proj)
    assert psi == proj
    degs, degs_other = list(degs), list(degs_other)
    degs.sort()
    degs_other.sort()
    return degs, degs_other

def get_glued_permutations(n, m):
    w, w_tilde = clean.canonical_words(n, m)
    phi, formula = clean.commute(w, w_tilde, 4)
    letters = ''.join(w)
    dict_form = defaultdict(int)
    for z1, z2, val in formula:
        new_z1 = []
        for word in z1:
            for l in letters:
                word = word.replace(l, '')
            new_z1.append(word)
        new_z2 = []
        for word in z2:
            for l in letters:
                word = word.replace(l, '')
            new_z2.append(word)
        new = tuple(new_z1), tuple(new_z2)
        dict_form[new] += val
    assert set(dict_form.values()) == {0}
    keys = [(tuple(c for c in key[0] if c != ''),
                    tuple(c for c in key[1] if c != '')) for key in dict_form.keys()]
    keys = sorted(list(set(keys)))
    groups = defaultdict(list)
    for z1, z2 in keys:
        z2_sorted = tuple(sorted(list(z2)))
        groups[(z1, z2_sorted)].append((z1, z2))
    groups = list(groups.values())
    groups.sort(key=lambda group: (-len(group[0][0]) - len(group[0][1]), sorted(group[0][0] + group[0][1])))
    groups = [[(';'.join(v[0]), ';'.join(v[1])) for v in group] for group in groups]
    grps = []
    for group in groups:
        grps.append(tabulate(sorted(group), headers=['w', 'w_tilde']))
    return grps

def phi_contained_in_first(n):
    t = double_bracket.canonical_tensor(n, 1)
    letter = clean.ALPHABET[n]
    s = double_bracket.double_bracket(t)
    phis = []
    for a, b, val in s:
        phi, psi = clean.commute(a, b, 1)
        #clean.pretty_print(phi)
        phis.append(clean.multiply(phi, val))
    #print("FINAL")
    phi = clean.add(*phis)
    double_bracket.print_tensor(s)
    for a, b, val in phi:
        assert letter in ''.join(b)
    #clean.pretty_print(phi)

def phi_bar(n, m, index=1):
    t = double_bracket.canonical_tensor(n, m)
    s = double_bracket.double_bracket(t)
    double_bracket.print_tensor(s, sort=False)
    phis = []
    for a, b, val in s:
        phi, psi = clean.commute(a, b, index)
        print(a, b)
        clean.pretty_print(phi)
        phis.append(clean.multiply(phi, val))
    print("FINAL")
    phi = clean.add(*phis)
    clean.pretty_print(phi)


# we will work with linear combinations, so let's implement such thing
def term_to_linear_combination(w):
    return [(w, 1)]

def partitions(n, I=1):
    yield (n,)
    for i in range(I, n//2 + 1):
        for p in partitions(n-i, i):
            yield (i,) + p
def compositions(n):
    for p in partitions(n):
        for q in set(permutations(p)):
            yield q

def split_word(w, nu):
    assert ALGEBRA != 'C^N'
    idx = 0
    x = []
    for k in nu:
        m = ''
        for letter in w[idx:idx+k]:
            m = double_bracket.multiply(m, letter)
        x.append(m)
        idx += k
    return tuple(x)
def shift(w, c):
    # here c corresponds to -c from paper
    length = len(w)
    if length == 0:
        return [(tuple(), 1)] #empty word - nothing happens to it
    result = []
    for nu in compositions(length):
        x = split_word(w, nu)
        result.append((x, c**(length - len(nu))))
    return result

def shift_quadratic(quadratic, c):
    result = []
    for x1, x2, val in quadratic:
        first = shift(x1, c)
        second = shift(x2, c)
        for y1, v1 in first:
            for y2, v2 in second:
                result.append((y1, y2, v1*v2*val))
    return clean.add(result)

def commute_shifted(w, w_tilde, c):
    # first shift, then commute
    first = shift(w, c)
    second = shift(w_tilde, c)
    phis = []
    psis = []
    for (x1, v1) in first:
        for (x2, v2) in second:
            phi, psi = commute(x1, x2, 4)
            phi = clean.multiply(phi, v1*v2)
            psi = clean.multiply(psi, v1*v2)
            phis.append(phi)
            psis.append(psi)
    return sorted(clean.add(*phis)), sorted(clean.add(*psis))

def shift_commuted(w, w_tilde, c):
    # first commute, then shift
    phi, psi = commute(w, w_tilde, 4)
    return sorted(shift_quadratic(phi, c)), sorted(shift_quadratic(psi, c))

def is_shift_automorphism(n, m, k):
    for u in range(1, m+1):
        for v in range(1,n+1):
            print(u, v)
            for l in range(-k, k+1):
                w, w_tilde = clean.canonical_words(u, v)
                phi_cs, psi_cs = commute_shifted(w, w_tilde, l)
                phi_sc, psi_sc = shift_commuted(w, w_tilde, l)
                assert phi_cs == phi_sc
                assert psi_cs == psi_sc