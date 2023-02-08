from itertools import product
from misha import Formula_I
from clean import commute
# print('Use the following format: numbers (free algebra vars) separated by semi-colons ("train")')
# print('For example: 1;2;1')
# w = tuple(input("Enter w: ").split(';'))
# w_tilde = tuple(input("Enter w_tilde: ").split(';'))
# pretty_print_formulas(w, w_tilde)
from tests import compare_my_and_misha, compare_two_versions

compare_two_versions(4,4)
