from clean import ALGEBRA, canonical_words, free_algebra_commute
from tabulate import tabulate
from collections import defaultdict
import numpy as np
import scipy


# we will represent elements of T(A)xT(A) similarly as lists of algebra elements, and each list has a coefficient
def multiply(a, b):
    if ALGEBRA == 'free':
        return a + b
    elif ALGEBRA == 'commutative':
        return ''.join(sorted(a + b))
    elif ALGEBRA == 'C^N':
        if a == '':
            return b
        if b == '':
            return a
        if a == b:
            return a
        return 0
    else:
        raise Exception(f"Algebra {ALGEBRA} not implemented. Choose C^N or free or commutative")


def reduce(t):
    dict_form = defaultdict(int)
    for b1, b2, val in t:
        dict_form[(b1, b2)] += val
    result = []
    for (b1, b2), val in dict_form.items():
        if val != 0:
            result.append((b1, b2, val))
    return result


def double_bracket_recursive(t):
    result = []
    for basis1, basis2, val in t:
        # handle base case
        if len(basis1) == 1 and len(basis2) == 1:
            a = basis1[0]
            b = basis2[0]
            ab = multiply(a, b)
            ba = multiply(b, a)
            if ab != 0:
                result.append(([], [ab], val))
            if ba != 0:
                result.append(([ba], [], -val))
        elif len(basis1) == 0 or len(basis2) == 0:
            # T(A) unit gives zero bracket with everything
            continue
        elif len(basis2) > 1:
            length = len(basis2)
            b1, b2 = basis2[:length // 2], basis2[length // 2:]
            first_summand = double_bracket_recursive([(basis1, b1, val)])
            second_summand = double_bracket_recursive([(basis1, b2, val)])
            result.extend((t1, t2 + b2, v) for t1, t2, v in first_summand)
            result.extend((b1 + t1, t2, v) for t1, t2, v in second_summand)
        elif len(basis1) > 1:
            length = len(basis1)
            b1, b2 = basis1[:length // 2], basis1[length // 2:]
            first_summand = double_bracket_recursive([(b1, basis2, val)])
            second_summand = double_bracket_recursive([(b2, basis2, val)])
            result.extend((t1 + b2, t2, v) for t1, t2, v in first_summand)
            result.extend((t1, b1 + t2, v) for t1, t2, v in second_summand)
    return reduce(result)


def double_bracket(t):
    result = []
    for a, b, val in t:
        for i in range(len(a)):
            for j in range(len(b)):
                a_left = a[:i]
                a_right = a[i + 1:]
                b_left = b[:j]
                b_right = b[j + 1:]
                ab = multiply(a[i], b[j])
                ba = multiply(b[j], a[i])
                if ab != 0:
                    result.append((b_left + a_right, a_left + (ab,) + b_right, val))
                if ba != 0:
                    result.append((b_left + (ba,) + a_right, a_left + b_right, -val))
    return reduce(result)


# computes iterations until the result becomes zero
def compute_iterations(t):
    iterations = []
    while True:
        t = double_bracket(t)
        if len(t) == 0:
            break
        iterations.append(t)

    return iterations


def compute_iterations_recursive(t):
    iterations = []
    while True:
        t = double_bracket_recursive(t)
        if len(t) == 0:
            break
        iterations.append(t)

    return iterations


def project_first_component_zero_degree(t):
    result = []
    for a, b, val in t:
        if len(a) == 0:
            result.append((a, b, val))
    return result


def project_second_component_zero_degree(t):
    result = []
    for a, b, val in t:
        if len(b) == 0:
            result.append((a, b, val))
    return result


def project(t, p, q):
    # projects to degree p first and degree q second
    result = []
    for a, b, val in t:
        if len(a) == p and len(b) == q:
            result.append((a, b, val))
    return result


def compute_brackets(n, m):
    t = canonical_tensor(n, m)
    all_iterations = compute_iterations(t)
    brackets = []
    for i, s in enumerate(all_iterations):
        projection_first = project_first_component_zero_degree(s)
        projection_second = project_second_component_zero_degree(s)
        if len(projection_first) + len(projection_second) < len(s):
            brackets.append((f'Iteration {i + 1}', s))
        if len(projection_first) > 0:
            brackets.append((f'Iteration {i + 1}; projection 0,-', project_first_component_zero_degree(s)))
        if len(projection_second) > 0:
            brackets.append((f'Iteration {i + 1}; -,0', project_second_component_zero_degree(s)))
    return brackets


def compute_brackets_and_all_projections(n, m):
    t = canonical_tensor(n, m)
    iterations = compute_iterations(t)
    degree = n + m - 1
    brackets = []
    for i, s in enumerate(iterations):
        if i == 0:
            brackets.append((f'Iteration {1}', s))
            degree -= 1
            continue
        # brackets.append((f'Iteration {i+1}', s))
        for p in range(degree + 1):
            q = degree - p
            proj = project(s, p, q)
            if len(proj) > 0:
                brackets.append((f'Iteration {i + 1}, projection {p}, {q}', proj))

        degree -= 1
    return brackets


def print_tensor(t, sort=True):
    string_form = [(';'.join(a), ';'.join(b), val) for a, b, val in t]
    if sort:
        string_form.sort(key=lambda item: (-item[2], item[0], item[1]))
    print(tabulate(string_form, headers=['T(A)', 'T(A)', 'coef']))


def canonical_tensor(n, m):
    w, w_tilde = canonical_words(n, m)
    return [(w, w_tilde, 1)]


def express_formula(formula, brackets):
    basis2index = {}
    index = 0
    for z1, z2, val in formula:
        if (z1, z2) not in basis2index:
            basis2index[z1, z2] = index
            index += 1
    for name, bracket in brackets:
        for a, b, val in bracket:
            if (a, b) not in basis2index:
                basis2index[a, b] = index
                index += 1
    # now construct arrays
    formula_array = np.zeros(len(basis2index), dtype='int')
    brackets_matrix = np.zeros((len(basis2index), len(brackets)))
    for z1, z2, val in formula:
        formula_array[basis2index[z1, z2]] = val
    for i, (name, bracket) in enumerate(brackets):
        for a, b, val in bracket:
            brackets_matrix[basis2index[a, b], i] = val
    solution, err, rank, _ = np.linalg.lstsq(brackets_matrix, formula_array, rcond=None)
    assert rank == len(brackets)
    # print(err[0])
    rounded_solution = np.clip(solution, -1, 1).round().astype('int')
    l2_error = np.sum((brackets_matrix @ rounded_solution - formula_array) ** 2)
    terms = []
    for i, val in enumerate(rounded_solution):
        if val != 0:
            terms.append((brackets[i][0], val))
    # print(solution, rounded_solution)
    return terms, l2_error, err


# not debugged yet!
# def express_formula_sparse(formula, brackets):
#     basis2index = {}
#     index = 0
#     for z1, z2, val in formula:
#         if (z1, z2) not in basis2index:
#             basis2index[z1, z2] = index
#             index += 1
#     for name, bracket in brackets:
#         for a, b, val in bracket:
#             if (tuple(a), tuple(b)) not in basis2index:
#                 basis2index[tuple(a), tuple(b)] = index
#                 index += 1
#     # now construct arrays
#     formula_array = scipy.sparse.dok_array((len(basis2index),1), dtype='int')
#     brackets_matrix = scipy.sparse.dok_array((len(basis2index), len(brackets)))
#     for z1, z2, val in formula:
#         formula_array[basis2index[z1, z2], 0] = val
#     formula_array = formula_array.tocsr()
#     brackets_matrix = brackets_matrix.tocsr()
#     for i, (name, bracket) in enumerate(brackets):
#         for a, b, val in bracket:
#             brackets_matrix[basis2index[tuple(a), tuple(b)], i] = val
#     solution= scipy.sparse.linalg.spsolve(brackets_matrix, formula_array)
#     err = np.sum((brackets_matrix @ solution - formula_array)**2)
#     solution = solution.todense()
#     rounded_solution = np.clip(solution, -1, 1).round().astype('int')
#     #l2_error = np.sum((brackets_matrix @ rounded_solution - formula_array) ** 2)
#     terms = []
#     for i, val in enumerate(rounded_solution):
#         if val != 0:
#             terms.append((brackets[i][0], val))
#     #print(solution, rounded_solution)
#     return terms, err#l2_error

def express_through_double_brackets(n, m, all_projections=False):
    phi, psi = free_algebra_commute(n, m, 1)
    if all_projections:
        brackets = compute_brackets_and_all_projections(n, m)
    else:
        brackets = compute_brackets(n, m)
    expressed_phi, err_phi, linear_err_phi = express_formula(phi, brackets)
    expressed_psi, err_psi, linear_err_psi = express_formula(psi, brackets)
    return expressed_phi, expressed_psi, err_phi + err_psi, round(linear_err_phi[0] + linear_err_psi[0])


def compute_and_print_expression_of_commutator_via_double_brackets(n, m, all_projections=False):
    phi, psi, round_err, err = express_through_double_brackets(n, m, all_projections)
    phi.sort(key=lambda x: -x[1])
    psi.sort(key=lambda x: -x[1])

    print("Error (linear; rounding):", err, round_err)
    print('phi')
    for term, sign in phi:
        print('+' if sign == 1 else '-', term)
    print('psi')
    for term, sign in psi:
        print('+' if sign == 1 else '-', term)


def express_commutator_via_double_brackets_io():
    print(f'Computing in {ALGEBRA} algebra')
    all_projections = input("Use all projections? (y/n):")
    if all_projections == "y":
        all_projections = True
    else:
        all_projections = False
    n = int(input('Enter the number of variables for w:'))
    m = int(input('Enter the number of variables for w_tilde:'))
    t = canonical_tensor(n, m)
    print("Computing commutator and iterations/projections of double bracket for")
    print_tensor(t)
    print("Result:")
    compute_and_print_expression_of_commutator_via_double_brackets(n, m, all_projections)


def compute_brackets_io():
    print(f'Computing in {ALGEBRA} algebra')
    all_projections = input("Use all projections? (y/n):")
    if all_projections == "y":
        all_projections = True
    else:
        all_projections = False
    n = int(input('Enter the number of variables for w:'))
    m = int(input('Enter the number of variables for w_tilde:'))
    t = canonical_tensor(n, m)
    print("Computing iterations/projections of double bracket for")
    print_tensor(t)
    print("Result:")
    if all_projections:
        brackets = compute_brackets_and_all_projections(n, m)
    else:
        brackets = compute_brackets(n, m)
    for name, s in brackets:
        print(name)
        print_tensor(s)


def compute_first_psi_term_explicitly(n, m):
    alpha, beta = canonical_words(n, m)
    result = []
    for s in range(m):
        for r in range(n):
            for p in range(r + 1, n):
                left = alpha[:r] + alpha[p + 1:]

                right = beta[:s] + \
                        (multiply(multiply(beta[s], alpha[p]), alpha[r]),) + \
                        alpha[r + 1:p] + beta[s + 1:]

                result.append((left, right, 1))

                left = alpha[:r] + alpha[p + 1:]

                right = beta[:s] + alpha[r + 1:p] + \
                        (multiply(multiply(alpha[p], alpha[r]), beta[s]),) + \
                        beta[s + 1:]

                result.append((left, right, -1))
        for q in range(n):
            for r in range(q + 1, n):
                for p in range(r + 1, n):
                    left_first = alpha[:q] + \
                                 (multiply(alpha[q], alpha[p]),) + \
                                 alpha[p + 1:]

                    left_second = alpha[:q] + alpha[p + 1:]

                    right = beta[:s] + \
                            alpha[r + 1:p] + alpha[q + 1:r] + \
                            (multiply(alpha[r], beta[s]),) + \
                            beta[s + 1:]

                    result.append((left_first, right, 1))

                    right = beta[:s] + \
                            alpha[r + 1:p] + \
                            (multiply(alpha[p], alpha[q]),) + \
                            alpha[q + 1:r] + \
                            (multiply(alpha[r], beta[s]),) + \
                            beta[s + 1:]

                    result.append((left_second, right, -1))

                    right = beta[:s] + \
                            (multiply(beta[s], alpha[r]),) + \
                            alpha[r + 1:p] + \
                            (multiply(alpha[p], alpha[q]),) + \
                            alpha[q + 1:r] + \
                            beta[s + 1:]

                    result.append((left_second, right, 1))

                    right = beta[:s] + \
                            (multiply(beta[s], alpha[r]),) + \
                            alpha[r + 1:p] + alpha[q + 1:r] + \
                            beta[s + 1:]
                    result.append((left_first, right, -1))
    return result
