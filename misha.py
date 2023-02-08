def Formula_I(w, w_tilde):
    Phi_1 = {}
    Psi_1 = {}

    if len(w) > 1:
        w_1 = w[:1]  # w_1 --- первая буква слова w, но представленная как tuple.
        w_2 = w[1:]  # w_2 --- все оставшиеся буквы.
        Phi_4_1, Psi_4_1 = Formula_IV(w_1, w_tilde)
        Phi_3_2, Psi_3_2 = Formula_III(w_2, w_tilde)

        for (z_1, z_2) in Phi_4_1:
            Phi_1[(z_2 + w_2, z_1)] = Phi_1.get((z_2 + w_2, z_1), 0) + Phi_4_1[(z_1, z_2)]

            Phi_4_z, Psi_4_z = Formula_IV(z_1, z_2 + w_2)
            for (x_1, x_2) in Phi_4_z:
                Psi_1[(x_1, x_2)] = Psi_1.get((x_1, x_2), 0) + Phi_4_1[(z_1, z_2)] * Phi_4_z[(x_1, x_2)]
            for (x_1, x_2) in Psi_4_z:
                Phi_1[(x_1, x_2)] = Phi_1.get((x_1, x_2), 0) + Phi_4_1[(z_1, z_2)] * Psi_4_z[(x_1, x_2)]

        for (z_1, z_2) in Phi_3_2:
            Phi_1[(z_2, w_1 + z_1)] = Phi_1.get((z_2, w_1 + z_1), 0) + Phi_3_2[(z_1, z_2)]

            Phi_4_z, Psi_4_z = Formula_IV(w_1 + z_1, z_2)
            for (x_1, x_2) in Phi_4_z:
                Psi_1[(x_1, x_2)] = Psi_1.get((x_1, x_2), 0) + Phi_3_2[(z_1, z_2)] * Phi_4_z[(x_1, x_2)]
            for (x_1, x_2) in Psi_4_z:
                Phi_1[(x_1, x_2)] = Phi_1.get((x_1, x_2), 0) + Phi_3_2[(z_1, z_2)] * Psi_4_z[(x_1, x_2)]

        for (z_1, z_2) in Psi_3_2:
            Psi_1[(w_1 + z_1, z_2)] = Psi_1.get((w_1 + z_1, z_2), 0) + Psi_3_2[(z_1, z_2)]

        for (z_1, z_2) in Psi_4_1:
            Psi_1[(z_2 + w_2, z_1)] = Psi_1.get((z_2 + w_2, z_1), 0) + Psi_4_1[(z_1, z_2)]

            Phi_4_z, Psi_4_z = Formula_IV(z_1, z_2 + w_2)
            for (x_1, x_2) in Phi_4_z:
                Phi_1[(x_1, x_2)] = Phi_1.get((x_1, x_2), 0) + Psi_4_1[(z_1, z_2)] * Phi_4_z[(x_1, x_2)]
            for (x_1, x_2) in Psi_4_z:
                Psi_1[(x_1, x_2)] = Psi_1.get((x_1, x_2), 0) + Psi_4_1[(z_1, z_2)] * Psi_4_z[(x_1, x_2)]

    elif len(w) == 1 and len(w_tilde) > 1:
        Phi_4_reverse, Psi_4_reverse = Formula_IV(w_tilde, w)

        for (z_1, z_2) in Phi_4_reverse:
            Phi_1[(z_1, z_2)] = - Phi_4_reverse[(z_1, z_2)]
        for (z_1, z_2) in Psi_4_reverse:
            Psi_1[(z_1, z_2)] = - Psi_4_reverse[(z_1, z_2)]

    else:  # i.e. if len(w) == 1 and len(w_tilde) == 1
        if w == w_tilde:
            assert len(w) == 1
            Phi_1[((), w)] = 1
            Phi_1[(w, ())] = -1

    return Phi_1, Psi_1


def Formula_II(w, w_tilde):
    Phi_2 = {}
    Psi_2 = {}

    if len(w) > 1:
        w_1 = w[:1]  # w_1 --- первая буква слова w, но представленная как tuple.
        w_2 = w[1:]  # w_2 --- все оставшиеся буквы.

        Phi_4_1, Psi_4_1 = Formula_IV(w_1, w_tilde)
        Phi_3_2, Psi_3_2 = Formula_III(w_2, w_tilde)

        for (z_1, z_2) in Phi_4_1:
            Phi_2[(z_2 + w_2, z_1)] = Phi_2.get((z_2 + w_2, z_1), 0) + Phi_4_1[(z_1, z_2)]

            Phi_2_z, Psi_2_z = Formula_II(z_1, z_2 + w_2)
            for (x_1, x_2) in Phi_2_z:
                Psi_2[(x_1, x_2)] = Psi_2.get((x_1, x_2), 0) + Phi_4_1[(z_1, z_2)] * Phi_2_z[(x_1, x_2)]
            for (x_1, x_2) in Psi_2_z:
                Phi_2[(x_1, x_2)] = Phi_2.get((x_1, x_2), 0) + Phi_4_1[(z_1, z_2)] * Psi_2_z[(x_1, x_2)]

        for (z_1, z_2) in Phi_3_2:
            Phi_2[(z_2, w_1 + z_1)] = Phi_2.get((z_2, w_1 + z_1), 0) + Phi_3_2[(z_1, z_2)]

            Phi_2_z, Psi_2_z = Formula_II(w_1 + z_1, z_2)
            for (x_1, x_2) in Phi_2_z:
                Psi_2[(x_1, x_2)] = Psi_2.get((x_1, x_2), 0) + Phi_3_2[(z_1, z_2)] * Phi_2_z[(x_1, x_2)]
            for (x_1, x_2) in Psi_2_z:
                Phi_2[(x_1, x_2)] = Phi_2.get((x_1, x_2), 0) + Phi_3_2[(z_1, z_2)] * Psi_2_z[(x_1, x_2)]

        for (z_1, z_2) in Psi_4_1:
            Psi_2[(z_1, z_2 + w_2)] = Psi_2.get((z_1, z_2 + w_2), 0) + Psi_4_1[(z_1, z_2)]

        for (z_1, z_2) in Psi_3_2:
            Psi_2[(z_2, w_1 + z_1)] = Psi_2.get((z_2, w_1 + z_1), 0) + Psi_3_2[(z_1, z_2)]

            Phi_2_z, Psi_2_z = Formula_II(w_1 + z_1, z_2)  # это мы уже вычисляли btw
            for (x_1, x_2) in Phi_2_z:
                Phi_2[(x_1, x_2)] = Phi_2.get((x_1, x_2), 0) + Psi_3_2[(z_1, z_2)] * Phi_2_z[(x_1, x_2)]
            for (x_1, x_2) in Psi_2_z:
                Psi_2[(x_1, x_2)] = Psi_2.get((x_1, x_2), 0) + Psi_3_2[(z_1, z_2)] * Psi_2_z[(x_1, x_2)]

    elif len(w) == 1 and len(w_tilde) > 1:
        Phi_3_reverse, Psi_3_reverse = Formula_III(w_tilde, w)

        for (z_1, z_2) in Phi_3_reverse:
            Phi_2[(z_1, z_2)] = - Phi_3_reverse[(z_1, z_2)]
        for (z_1, z_2) in Psi_3_reverse:
            Psi_2[(z_1, z_2)] = - Psi_3_reverse[(z_1, z_2)]

    else:  # i.e. if len(w) == 1 and len(w_tilde) == 1
        if w == w_tilde:
            assert len(w) == 1
            Phi_2[((), w)] = 1
            Phi_2[(w, ())] = -1

    return Phi_2, Psi_2


def Formula_III(w, w_tilde):
    Phi_3 = {}
    Psi_3 = {}

    if len(w) > 1:
        w_1 = w[:1]  # w_1 --- первая буква слова w, но представленная как tuple.
        w_2 = w[1:]  # w_2 --- все оставшиеся буквы.
        Phi_4_1, Psi_4_1 = Formula_IV(w_1, w_tilde)
        Phi_3_2, Psi_3_2 = Formula_III(w_2, w_tilde)

        for (z_1, z_2) in Phi_4_1:
            Phi_3[(z_1, z_2 + w_2)] = Phi_3.get((z_1, z_2 + w_2), 0) + Phi_4_1[(z_1, z_2)]

        for (z_1, z_2) in Phi_3_2:
            Phi_3[(w_1 + z_1, z_2)] = Phi_3.get((w_1 + z_1, z_2), 0) + Phi_3_2[(z_1, z_2)]

        for (z_1, z_2) in Psi_3_2:
            Psi_3[(w_1 + z_1, z_2)] = Psi_3.get((w_1 + z_1, z_2), 0) + Psi_3_2[(z_1, z_2)]

        for (z_1, z_2) in Psi_4_1:
            Psi_3[(z_2 + w_2, z_1)] = Psi_3.get((z_2 + w_2, z_1), 0) + Psi_4_1[(z_1, z_2)]

            Phi_2_z, Psi_2_z = Formula_II(z_1, z_2 + w_2)
            for (x_1, x_2) in Phi_2_z:
                Phi_3[(x_1, x_2)] = Phi_3.get((x_1, x_2), 0) + Psi_4_1[(z_1, z_2)] * Phi_2_z[(x_1, x_2)]
            for (x_1, x_2) in Psi_2_z:
                Psi_3[(x_1, x_2)] = Psi_3.get((x_1, x_2), 0) + Psi_4_1[(z_1, z_2)] * Psi_2_z[(x_1, x_2)]

    elif len(w) == 1 and len(w_tilde) > 1:
        Phi_2_reverse, Psi_2_reverse = Formula_II(w_tilde, w)

        for (z_1, z_2) in Phi_2_reverse:
            Phi_3[(z_1, z_2)] = - Phi_2_reverse[(z_1, z_2)]
        for (z_1, z_2) in Psi_2_reverse:
            Psi_3[(z_1, z_2)] = - Psi_2_reverse[(z_1, z_2)]

    else:  # i.e. if len(w) == 1 and len(w_tilde) == 1
        if w == w_tilde:
            assert len(w) == 1
            Phi_3[(w, ())] = 1
            Phi_3[((), w)] = -1

    return Phi_3, Psi_3


def Formula_IV(w, w_tilde):
    Phi_4 = {}
    Psi_4 = {}

    if len(w) > 1:
        w_1 = w[:1]  # w_1 --- первая буква слова w, но представленная как tuple.
        w_2 = w[1:]  # w_2 --- все оставшиеся буквы.
        Phi_4_1, Psi_4_1 = Formula_IV(w_1, w_tilde)
        Phi_3_2, Psi_3_2 = Formula_III(w_2, w_tilde)

        for (z_1, z_2) in Phi_4_1:
            Phi_4[(z_1, z_2 + w_2)] = Phi_4.get((z_1, z_2 + w_2), 0) + Phi_4_1[(z_1, z_2)]

        for (z_1, z_2) in Phi_3_2:
            Phi_4[(w_1 + z_1, z_2)] = Phi_4.get((w_1 + z_1, z_2), 0) + Phi_3_2[(z_1, z_2)]

        for (z_1, z_2) in Psi_4_1:
            Psi_4[(z_1, z_2 + w_2)] = Psi_4.get((z_1, z_2 + w_2), 0) + Psi_4_1[(z_1, z_2)]

        for (z_1, z_2) in Psi_3_2:
            Psi_4[(z_2, w_1 + z_1)] = Psi_4.get((z_2, w_1 + z_1), 0) + Psi_3_2[(z_1, z_2)]

            Phi_4_z, Psi_4_z = Formula_IV(w_1 + z_1, z_2)
            for (x_1, x_2) in Phi_4_z:
                Phi_4[(x_1, x_2)] = Phi_4.get((x_1, x_2), 0) + Psi_3_2[(z_1, z_2)] * Phi_4_z[(x_1, x_2)]
            for (x_1, x_2) in Psi_4_z:
                Psi_4[(x_1, x_2)] = Psi_4.get((x_1, x_2), 0) + Psi_3_2[(z_1, z_2)] * Psi_4_z[(x_1, x_2)]

    elif len(w) == 1 and len(w_tilde) > 1:
        Phi_1_reverse, Psi_1_reverse = Formula_I(w_tilde, w)

        for (z_1, z_2) in Phi_1_reverse:
            Phi_4[(z_1, z_2)] = - Phi_1_reverse[(z_1, z_2)]
        for (z_1, z_2) in Psi_1_reverse:
            Psi_4[(z_1, z_2)] = - Psi_1_reverse[(z_1, z_2)]

    else:  # i.e. if len(w) == 1 and len(w_tilde) == 1
        if w == w_tilde:
            assert len(w) == 1
            Phi_4[(w, ())] = 1
            Phi_4[((), w)] = -1

    return Phi_4, Psi_4