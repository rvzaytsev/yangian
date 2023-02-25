from clean import ALGEBRA, canonical_words, free_algebra_commute
from tabulate import tabulate
from collections import defaultdict
import numpy as np

# we will represent elements of T(A)xT(A) similarly as lists of algebra elements, and each list has a coefficient
def multiply(a, b):
    if ALGEBRA == 'free':
        return a+b
    elif ALGEBRA == 'commutative':
        return ''.join(sorted(a + b))
    elif ALGEBRA == 'C^N':
        if a == b:
            return a
        else:
            return 0
    else:
        raise Exception(f"Algebra {ALGEBRA} not implemented. Choose C^N or free or commutative")


def reduce(t):
    dict_form = defaultdict(int)
    for b1, b2, val in t:
        dict_form[(tuple(b1), tuple(b2))] += val
    result = []
    for (b1, b2), val in dict_form.items():
        if val != 0:
            result.append((list(b1), list(b2), val))
    return result
def double_bracket_recursive(t):
    result = []
    for basis1, basis2, val in t:
        # handle base case
        if len(basis1) == 1 and len(basis2) == 1:
            a = basis1[0]
            b = basis2[0]
            ab = multiply(a,b)
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
            result.extend((t1, t2+b2, v) for t1, t2, v in first_summand)
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
                a_right = a[i+1:]
                b_left = b[:j]
                b_right = b[j+1:]
                ab = multiply(a[i], b[j])
                ba = multiply(b[j], a[i])
                if ab != 0:
                    result.append((b_left+a_right, a_left + [ab] + b_right, val))
                if ba != 0:
                    result.append((b_left + [ba] + a_right,a_left + b_right, -val))
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

def compute_brackets(n, m):
    t = canonical_tensor(n, m)
    all_iterations = [t] + compute_iterations(t)
    brackets = []
    for i, s in enumerate(all_iterations):
        projection_first = project_first_component_zero_degree(s)
        projection_second = project_second_component_zero_degree(s)
        if len(projection_first) + len(projection_second) < len(s):
            brackets.append((f'Iteration {i}', s))
        if len(projection_first) > 0:
            brackets.append((f'Iteration {i}; projection 0,-', project_first_component_zero_degree(s)))
        if len(projection_second) > 0:
            brackets.append((f'Iteration {i}; -,0', project_second_component_zero_degree(s)))
    return brackets
def print_tensor(t):
    string_form = [(';'.join(a), ';'.join(b), val) for a, b, val in t]
    string_form.sort(key=lambda item: (-item[2], item[0], item[1]))
    print(tabulate(string_form, headers=['T(A)', 'T(A)', 'coef']))

def canonical_tensor(n, m):
    w, w_tilde = canonical_words(n, m)
    return [(list(w), list(w_tilde), 1)]


def express_formula(formula, brackets):
    basis2index = {}
    index = 0
    for z1, z2, val in formula:
        if (z1, z2) not in basis2index:
            basis2index[z1, z2] = index
            index += 1
    for name, bracket in brackets:
        for a, b, val in bracket:
            if (tuple(a), tuple(b)) not in basis2index:
                basis2index[tuple(a), tuple(b)] = index
                index += 1
    # now construct arrays
    formula_array = np.zeros(len(basis2index), dtype='int')
    brackets_matrix = np.zeros((len(basis2index), len(brackets)))
    for z1, z2, val in formula:
        formula_array[basis2index[z1, z2]] = val
    for i, (name, bracket) in enumerate(brackets):
        for a, b, val in bracket:
            brackets_matrix[basis2index[tuple(a), tuple(b)], i] = val
    solution, err, rank, _ = np.linalg.lstsq(brackets_matrix, formula_array, rcond=None)
    assert rank == len(brackets)
    #print(err[0])
    rounded_solution = np.clip(solution, -1, 1).round().astype('int')
    l2_error = np.sum((brackets_matrix @ rounded_solution - formula_array) ** 2)
    terms = []
    for i, val in enumerate(rounded_solution):
        if val != 0:
            terms.append((brackets[i][0], val))
    #print(solution, rounded_solution)
    return terms, l2_error


def express_through_double_brackets(n, m):
    phi, psi = free_algebra_commute(n, m, 1)
    brackets = compute_brackets(n, m)
    expressed_phi, err_phi = express_formula(phi, brackets)
    expressed_psi, err_psi = express_formula(psi, brackets)
    return expressed_phi, expressed_psi, err_phi + err_psi

def compute_and_print_expression_of_commutator_via_double_brackets(n, m):
    phi, psi, err = express_through_double_brackets(n,m)
    phi.sort(key=lambda x:-x[1])
    psi.sort(key=lambda x:-x[1])
    print("Error:", err)
    print('phi')
    for term, sign in phi:
        print('+' if sign == 1 else '-', term)
    print('psi')
    for term, sign in psi:
        print('+' if sign == 1 else '-', term)