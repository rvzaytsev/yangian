from itertools import product
from clean import commute
from clean_copy import commute as commute_copy
from misha import Formula_I, Formula_II, Formula_III, Formula_IV
def normal_form(formula):
    formula_list = []
    for z1, z2, val in formula:
        formula_list.append([';'.join(z1), ';'.join(z2), val])
    formula_list.sort(key=lambda item: (-item[2], item[0], item[1]))
    return formula_list
def list_form(dict_form):
    formula = []
    for (z1, z2), val in dict_form.items():
        if val != 0:
            formula.append((z1, z2, val))
    return formula

def compare_my_and_misha(m, n):
    mishas = [Formula_I, Formula_II, Formula_III, Formula_IV]
    for w in product(('1', '2'), repeat=m):
        for w_tilde in product(('1', '2'), repeat=n):
            for i in range(4):
                phi_my, psi_my = commute(w, w_tilde, i+1)
                phi_my.sort()
                psi_my.sort()
                phi_misha, psi_misha = mishas[i](w, w_tilde)
                phi_misha, psi_misha = list_form(phi_misha), list_form(psi_misha)
                phi_misha.sort()
                psi_misha.sort()
                assert(phi_my == phi_misha and psi_my == psi_misha)

def compare_two_versions(m, n):
    for w in product(('1', '2'), repeat=4):
        for w_tilde in product(('1', '2'), repeat=4):
            for i in range(4):
                phi, psi = commute(w, w_tilde, i+1)
                phi.sort()
                psi.sort()
                phi_copy, psi_copy = commute_copy(w, w_tilde, i+1)
                phi_copy.sort()
                psi_copy.sort()
                assert(phi == phi_copy and psi == psi_copy)

