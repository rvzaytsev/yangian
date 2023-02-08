from clean import pretty_print_formulas, free_algebra_io
# free_algebra_io()
# print('Use the following format: numbers (free algebra vars) separated by semi-colons ("train")')
# print('For example: 1;2;1')
# w = tuple(input("Enter w: ").split(';'))
# w_tilde = tuple(input("Enter w_tilde: ").split(';'))
# pretty_print_formulas(w, w_tilde, 4)
from tests import compare_my_and_misha, compare_word_division_free_algebra, compare_word_division
# from clean import free_algebra_io
# free_algebra_io()
from hypothesis import print_coefficients_free_algebra, check_coefs
for n in range(1, 10):
    print_coefficients_free_algebra(n,n)
#check_coefs(4,4)