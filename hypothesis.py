from clean import free_algebra_commute, ALGEBRA, canonical_words, commute, add, negate
from itertools import permutations, product, combinations

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

# only support linear combinations with 2 coeffs :(
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