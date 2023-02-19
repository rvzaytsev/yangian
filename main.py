# print('Use the following format: numbers (free algebra vars) separated by semi-colons ("train")')
# print('For example: 1;2;1')
# w = tuple(input("Enter w: ").split(';'))
# w_tilde = tuple(input("Enter w_tilde: ").split(';'))
# pretty_print_formulas(w, w_tilde, 3)
from hypothesis import print_coefficients_free_algebra, check_coefs
#print_coefficients_free_algebra(5,5)
#check_coefs(5,5)
# print_coefficients_free_algebra(7,7)
from tests import compare_word_division_free_algebra, free_algebra_reduce_to_commutative_test
free_algebra_reduce_to_commutative_test(5,5)