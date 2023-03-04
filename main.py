# from clean import free_algebra_commute, pretty_print
# phi, psi = free_algebra_commute(3, 3, 1)
# pretty_print(psi)
#

# import double_bracket
# t=double_bracket.canonical_tensor(2,1)
# for name, t in double_bracket.compute_brackets(2,1):
#     print(name)
#     double_bracket.print_tensor(t)

# from double_bracket import compute_and_print_expression_of_commutator_via_double_brackets
# compute_and_print_expression_of_commutator_via_double_brackets(3,5)
import double_bracket
double_bracket.express_commutator_via_double_brackets_io()
# from hypothesis import calculate_permutation_sum
# from clean import pretty_print
# phi, psi = calculate_permutation_sum(1, 1)
# pretty_print(phi)
# pretty_print(psi)
# from clean import free_algebra_io
# free_algebra_io()
#
# double_bracket.print_tensor(
#     double_bracket.compute_iterations(double_bracket.canonical_tensor(4,2))[1])
# t = double_bracket.canonical_tensor(4,4)
# double_bracket.print_tensor(t)
# for i, s in enumerate(double_bracket.compute_iterations(t)):
#     print('iteration n.', i+1)
#     double_bracket.print_tensor(s)
    # print("projection 0, -")
    # double_bracket.print_tensor(double_bracket.project_first_component_zero_degree(s))
    # print("projection -, 0")
    # double_bracket.print_tensor(double_bracket.project_second_component_zero_degree(s))
# for x in t:
#     f=double_bracket.double_bracket_recursive([x])
#     s=double_bracket.double_bracket([x])
#     f.sort()
#     s.sort()
#     assert f == s
# double_bracket.print_tensor(double_bracket.double_bracket_recursive(t))
# double_bracket.print_tensor(double_bracket.double_bracket(t))
# import double_bracket
# t = double_bracket.canonical_tensor(2, 1)
# double_bracket.print_tensor(t)
# s = double_bracket.double_bracket(t)
# double_bracket.print_tensor(s)