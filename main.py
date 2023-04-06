from clean import free_algebra_projection
import clean, double_bracket
clean.project_total_degree_io()
n, m = 3,3
# phi, psi = clean.free_algebra_commute(n, m, 4)

# clean.free_algebra_io()
import hypothesis
# hypothesis.first_word_preserves_order(6,6)
# n, m = 3,3
# hypothesis.find_linear_combinations(n, m)
# print(minim)
# for i, c in enumerate(coef):
#     print("Realization", i+1)
#     for u, v, s in zip(*c):
#         print(s, ';'.join(u), ';'.join(v))
# x = ('0', '1', '2', '3')
# w1 = ('0',), ('12', '3')
# w2 = ('1',), ('2', '30')
# w3 = ('0', '1'), ('23',)
# w4 = ('01', '2'), ('3',)
# w5 = ('0',), ('1', '23')
# w6 = ('01',), ('2', '3')
# w7 = ('0', '12'), ('3',)
# w8 = ('1', '2'), ('30',)
# ws = [w1, w2, w3, w4, w5, w6, w7, w8]
# phi = []
# psi = []
# for i, (u, v) in enumerate(ws):
#     f, s = len(u), len(v)
#     uv = u + v
#     for j in range(3):
#         ab = uv[j:] + uv[:j]
#         a, b = ab[:f], ab[f:]
#         ph, ps = clean.commute(a, b, 4)
#         if i > 3:
#             ph, ps = clean.negate((ph, ps))
#         phi = clean.add(phi, ph)
#         psi = clean.add(psi, ps)
# clean.pretty_print(phi)
# clean.pretty_print(psi)
# formulas = hypothesis.get_all_tabloids(n+m, False)
# for w, w_tilde in formulas:
#     print(w, w_tilde)
# sol, err =hypothesis.find_linear_combinations(n, m)
# print(sol, err)
# phi, psi = res[('0', '1', '2', '3'), ('4','5')]
# clean.pretty_print(phi)
# phi, psi = clean.free_algebra_commute(n, m, 4)
# clean.pretty_print(phi)
# for w, w_tilde in hypothesis.get_all_tabloids(6, False):
#     print(w, w_tilde)
# print(clean.canonical_words(n, m))
# phi_1, psi_1 = clean.free_algebra_commute(n, m, 1)
# phi_2, psi_2 = clean.free_algebra_commute(n, m, 2)
# phi_3, psi_3 = clean.free_algebra_commute(n, m, 3)
# phi_4, psi_4 = clean.free_algebra_commute(n, m, 4)
# ij_kl = clean.add(psi_1, psi_3)
# kl_ij = clean.add(psi_2, psi_4)
# kj_il = clean.add(phi_1, phi_2)
# il_kj = clean.add(phi_3, phi_4)

# print("phi 1")
# clean.pretty_print(phi_1)
# print("phi 3")
# clean.pretty_print(phi_2)
# print("kj_il for phi_1 + phi_2")
# clean.pretty_print(kj_il)
# print("il_kj ")
# clean.pretty_print(il_kj)

# clean.pretty_print(psi_1)
# clean.pretty_print(psi_3)
# print(len(psi_1) + len(psi_3), len(ij_kl))
# print("ij_kl")
# clean.pretty_print(ij_kl)
# print("kl_ij")
# clean.pretty_print(kl_ij)


# for n in range(7):
#     for m in range(6):
#         degs, degs_other = hypothesis.compare_explicit_formula_psi_with_bracket(n, m)
#         print(n, m, degs, degs_other)
# from tests import test_explicit_psi
# for n in range(7):
#     for m in range(7):
#         print(n, m)
#         test_explicit_psi(n, m)
# w, w_tilde = clean.canonical_words(n, m)
# print(w, w_tilde)
# clean.pretty_print_formulas(w, w_tilde, 1)
# psi = double_bracket.compute_first_psi_term_explicitly(n, m)
# degs = set()
# for z1,z2, val in psi:
#     degs.add((len(z1), len(z2)))
# print("Degrees:", degs)
# clean.pretty_print(psi)
# second = clean.multiply(double_bracket.compute_iterations(double_bracket.canonical_tensor(n, m))[1], -1)
# proj = []
# for z1, z2, val in second:
#     if (len(z1), len(z2)) in degs:
#         proj.append((z1, z2, val))
# double_bracket.print_tensor(proj)
# clean.print_formulas_with_permutation(n, m)
# firsts, seconds = clean.get_word_permutations(n, m)
# for w, val in firsts:
#     if val == 1:
#         print("+",' '.join(w))
#     else:
#         print("-", ' '.join(w))
#hypothesis.first_word_no_gluing(3, 3)
# for n in range(7):
#     for m in range(7):
#         try:
#             hypothesis.detailed_sign_equality(n, m)
#         except:
#             print(n, m)

# for m in range(8):
#     for n in range(8):
#         if n + m < 12:
#             phi_proj, psi_proj, phi_deg, psi_deg = clean.nonzero_projections(n, m, index=4)
#             phi_deg.add(10**10)
#             psi_deg.add(10**10)
#             min_deg = min(min(phi_deg), min(psi_deg))
#             if min_deg == 10**10:
#                 min_deg = None
#             print(n, m, min_deg)
#phi_proj, psi_proj, phi_deg, psi_deg = clean.nonzero_projections(n, m, index=2)

# print(f"{n, m}")
# print(f"Degree {n+m}")
# print("Phi degrees:", sorted(phi_deg))
# print("Phi projections: " + ' '.join(map(str, phi_proj)))
# print("Psi degrees:", sorted(psi_deg))
# print("Psi projections: " + ' '.join(map(str, psi_proj)))

# w, w_tilde = clean.canonical_words(n,m)
# clean.pretty_print_formulas(w, w_tilde, 4)


#free_algebra_projection(n, m, r=(n+1)//2, s=(n+1)//2)
# d=(n+m-1)
# free_algebra_projection(n, m, r=d//2, s=d//2)