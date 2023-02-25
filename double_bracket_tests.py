from double_bracket import double_bracket, double_bracket_recursive, \
    canonical_tensor, compute_iterations, compute_iterations_recursive
from itertools import product
def compare_two_versions_of_double_bracket(n, m):
    t = canonical_tensor(n, m)
    first = double_bracket_recursive(t)
    first.sort()
    second = double_bracket(t)
    second.sort()
    assert first == second

def compare_two_versions_up_to(k):
    for n in range(k):
        for m in range(k):
            compare_two_versions_of_double_bracket(n, m)

def compare_iterations(n, m):
    t = canonical_tensor(n,m)
    first = compute_iterations(t)
    second = compute_iterations_recursive(t)
    assert len(first) == len(second)
    for f, s in zip(first, second):
        f.sort()
        s.sort()
        assert f == s

def compare_two_versions_for_product(n, m):
    for w in product(('1', '2'), repeat=n):
        for w_tilde in product(('1', '2'), repeat=m):
            t = [(list(w), list(w_tilde), 1)]
            first = double_bracket_recursive(t)
            second = double_bracket(t)
            first.sort()
            second.sort()
            assert first == second