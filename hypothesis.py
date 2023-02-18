from clean import free_algebra_commute

def print_coefficients_free_algebra(n, m):
    print(f"words are of length {n} and {m}")
    for i in range(4,5):
        print("Computing formula #"+str(i))
        phi, psi = free_algebra_commute(n, m, i)
        print(f"Total number of coefficients: {len(phi) + len(psi)}")
        phi_coefs = [item[2] for item in phi]
        psi_coefs = [item[2] for item in psi]
        print(f"Coefficients are in {set(phi_coefs) | set(psi_coefs)}")

def check_coefs(n, m):
    for u in range(1, n+1):
        for v in range(1, m+1):
            for i in range(1, 5):
                phi, psi = free_algebra_commute(u, v, i)
                phi_coefs = [item[2] for item in phi]
                psi_coefs = [item[2] for item in psi]
                assert (set(phi_coefs) | set(psi_coefs)).issubset({-1,1})

