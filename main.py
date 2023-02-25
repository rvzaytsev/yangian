# from hypothesis import find_good_linear_combinations
# combs, m = find_good_linear_combinations(3,6)
# print('Minimum found is', m)
# for comb in combs:
#     w1s, w2s, signs = comb
#     print(';'.join(w1s[0])+'_'+';'.join(w2s[0]), end='')
#     for j, charsign in enumerate(signs):
#         print(charsign+';'.join(w1s[j+1])+'_'+';'.join(w2s[j+1]), end='')
#     print()
from clean import free_algebra_io
free_algebra_io()

from double_bracket_tests import compare_iterations
import double_bracket
double_bracket.print_tensor(
    double_bracket.compute_iterations(double_bracket.canonical_tensor(3,1))[1])
# t = [(['1', '2', '1'], ['2','1', '2'], 1)]
# double_bracket.print_tensor(t)
# for i, s in enumerate(double_bracket.compute_iterations(t)):
#     print('iteration n.', i+1)
#     double_bracket.print_tensor(s)
#     print("projection 0, -")
#     double_bracket.print_tensor(double_bracket.project_first_component_zero_degree(s))
#     print("projection -, 0")
#     double_bracket.print_tensor(double_bracket.project_second_component_zero_degree(s))
# for x in t:
#     f=double_bracket.double_bracket_recursive([x])
#     s=double_bracket.double_bracket([x])
#     f.sort()
#     s.sort()
#     assert f == s
# double_bracket.print_tensor(double_bracket.double_bracket_recursive(t))
# double_bracket.print_tensor(double_bracket.double_bracket(t))

# print_tensor(double_bracket(t))