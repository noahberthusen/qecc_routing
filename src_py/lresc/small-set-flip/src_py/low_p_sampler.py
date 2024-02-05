import os
import sys
import pickle
import random
import math
import time

import matplotlib.pyplot as plt
import numpy as np
import decoder


class SplittingMethod():
    """
    This module contains the code required to conduct sampling of error
    configurations from the XZZX code under a biased noise model toward
    dephasing. The notation of the functions follows that used by the article
    'The role of entropy in topological quantum error correction'
    (arXiv:1812.05117) Appendix D as well as the original article where the
    method was conceived: 'Simulation of rare events in quantum error correction'
    (PhysRevA.88.062308).

    This code also utilises the sampled configurations together with the splitting
    method (PhysRevA.88.062308) to calculate logical error probabilities at very
    low physical error proabilities where Monte carlo simulations would be
    untractable.
    """

    def __init__(self, ccode, logical2):
        self.ccode = ccode
        self.logical2 = logical2

    """
    In the following functions we define the tools required to conduct the
    metropolis sampling algorithm and splitting method of PhysRevA.88.062308.
    """

    def g(self, x):
        """
        Detailed balance function as introduced by PhysRevA.88.062308.
        """
        return 1 / (1 + x)


    def prob_err(self, n, p, E):
        """
        Calculates the probability of observing a particular error configuration

        :param n: Number of qubits.
        :param p: Probability of an error on any one physical qubit.
        :param E: Weight of the observed error E.
        :return: Probability of observing error configuration
        """
        print(p, (p ** E) * ((1 - p) ** (n - E)))
        return (p ** E) * ((1 - p) ** (n - E))


    def prob_err1_div_prob_err2(self, p, E1, E2):
        """
        Calculates the probability of observing a particular error configuration
        on divided by the probability of observing a second error configuration
        on the coprime XZZX lattice given the number of high rate and low rate
        errors for each observed error.

        :param p: Probability of an error on any one physical qubit.
        :type p: float
        :param E1: Weight of the observed error E1.
        :param E2: Weight of the observed error E2.
        :return: Probability of observing error configuration 1 divided by
        probability of observing error configuration 2
        """
        return (p ** (E1 - E2)) * ((1 - p) ** (E2 - E1))


    def A_j(self, n, p1, p2, E):
        """
        Probability of an error configuration with H high rate errors and L
        low rate errors under the distribution (p2, p2_h, p2_l) divided by
        the probability distribution of the same error configuration under the
        distribution (p1, p1_h, p1_l).

        :param p1: Probability of an error on any one physical qubit under the
        first distribution.
        :param p2: Probability of an error on any one physical qubit under the
        second distribution.
        :return: Probability of error under first distribution over probability
        of error under the second distribution.
        :rtype: float
        """
        assert p1 < p2
        print('aj')
        # print((self.prob_err(n, p2, E) / self.prob_err(n, p1, E)), np.exp(np.log(self.prob_err(n, p2, E)) - np.log(self.prob_err(n, p1, E))))
        return np.exp(np.log(self.prob_err(n, p2, E)) - np.log(self.prob_err(n, p1, E)))


    def expec_G_CAj(self, n, p1, p2, C, Es):
        """
        Calculates the expectation value of the constant C times Aj, as defined in
        :func:`A_j`, numerically given an array of error configurations (see
        PhysRevA.88.062308 equation (8)).

        :param p1: Probability of an error on any one physical qubit under the
        first distribution.
        :type p1: float
        :param p2: Probability of an error on any one physical qubit under the
        second distribution.
        :type p2: float
        :param C: Constant as defined in PhysRevA.88.062308 equation (6).
        :type C: float
        :param E: Dictionary of error configurations in the form
        E = {(H1, L1): n1, ...., (Hk, Lk): nk} where each entry takes the form
        (# of high rate errs, # of low rate errs) : # of times such configuration
        appears.
        :type E: dict
        :return: Un-normalised expectation value of C times Aj calculated using
        the provided error configurations in E.
        :rtype: float
        """
        assert p1 < p2
        total = 0
        for E, i in Es.items():
            total += i * self.g(C * self.A_j(n, p1, p2, E))
        return total


    def expec_G_CAj1(self, n, p1, p2, C, Es):
        """
        Calculates the expectation value of one over the constant C times Aj,
        as defined in :func:`A_j`, numerically given an array of error
        configurations (see PhysRevA.88.062308 equation (8)).

        :param p1: Probability of an error on any one physical qubit under the
        first distribution.
        :type p1: float
        :param p2: Probability of an error on any one physical qubit under the
        second distribution.
        :type p2: float
        :param C: Constant as defined in PhysRevA.88.062308 equation (6).
        :type C: float
        :param E: Dictionary of error configurations in the form
        E = {(H1, L1): n1, ...., (Hk, Lk): nk} where each entry takes the form
        (# of high rate errs, # of low rate errs) : # of times such configuration
        appears.
        :type E: dict
        :return: Un-normalised expectation value of 1 / (C times Aj) calculated
        using the provided error configurations in E.
        :rtype: float
        """
        assert p1 < p2
        total = 0
        for E, i in Es.items():
            total += i * self.g(1 / (C * self.A_j(n, p1, p2, E)))
        return total


    def sample_E(self, p, num_iter, num_iter_cut):
        """
        Sample error configurations using the metropolis hastings algorithm
        presented in PhysRevA.88.062308.

        :param p: Probability of an error on any one physical qubit.
        :param N_cut: The number of simulations ran before we start to sample
        error configurations.
        :type N_cut: int
        """
        # Initialise in the required error state. <---------------- UNSURE WHAT WE WANT HERE
        mask = decoder.random_mask(self.ccode, 0) # <---------------- CHANGE THIS
#         E0 = ([(1, 17), (1, 34), (2, 38), (2, 53), (3, 26), (4, 37), (7, 54), (8, 3), (8, 27), (8, 47), (8, 53), (11, 2), (11, 15), (12, 4), (12, 29), (13, 55), (14, 6), (15, 27), (16, 45), (17, 27), (18, 0), (18, 17), (19, 34), (20, 8), (20, 37), (21, 1), (23, 26), (24, 30), (24, 37), (24, 48), (25, 59), (28, 16), (29, 40), (30, 8), (31, 13), (31, 22), (33, 12), (33, 49), (33, 50), (34, 40), (36, 49), (36, 54), (37, 23), (40, 3), (40, 13), (40, 50), (40, 52), (41, 35), (42, 20), (43, 42), (44, 10), (45, 12), (46, 46), (46, 55), (47, 39), (47, 50), (48, 31), 
# (48, 35), (48, 37), (49, 18), (49, 30), (49, 39), (49, 43), (49, 55), (51, 21), (51, 30), (52, 57), (53, 18), (57, 2), (58, 3), (58, 38)], [(1, 43), (2, 10), (3, 41), (3, 45), (4, 3), (5, 19), (9, 13), (9, 42), (13, 37), (14, 46), (15, 43), (16, 5), (16, 15), (19, 1), (19, 14), (21, 14), (21, 48), (25, 24), (27, 4), (28, 10), (28, 34), (29, 3), (31, 4), (31, 41), 
# (33, 3), (33, 40), (33, 49), (35, 13), (37, 2), (38, 6), (38, 10), (38, 36), (42, 24), (44, 35), (45, 35), (46, 48), (47, 
# 2), (47, 6), (48, 23), (48, 44)])
        E0 = ([(2, 0), (5, 8), (6, 2), (7, 7), (10, 2), (10, 8), (15, 15)], [(5, 4), (8, 0), (10, 7)])
        # E0 = ([(2, 5), (3, 6), (5, 5)], [(0, 2), (2, 2)])

        vv_xerror_array = np.array([[False for v2 in range(self.ccode.n)] for v1 in range(self.ccode.n)])
        cc_xerror_array = np.array([[False for c2 in range(self.ccode.m)] for c1 in range(self.ccode.m)])

        for (v1, v2) in E0[0]:
            vv_xerror_array[v1][v2] = not vv_xerror_array[v1][v2]
        for (c1, c2) in E0[1]:
            cc_xerror_array[c1][c2] = not cc_xerror_array[c1][c2]


        # Initialise the dictionary of sampled error configurations in the form:
        # E = {E1: n1, ...., Ek: nk}
        #   = { weight of error : # of times such an error configuration appears }
        Es = {}

        # Run simulations.
        for counter in range(num_iter):
            # Propose new error configuration by introducing an error onto a qubit.
            proposed_vv_xerror_array = np.copy(vv_xerror_array)
            proposed_cc_xerror_array = np.copy(cc_xerror_array)
            e = decoder.random_sized_error(self.ccode, 1)

            # Apply error to system
            for (v1, v2) in e[0]:
                proposed_vv_xerror_array[v1][v2] = not vv_xerror_array[v1][v2]
            for (c1, c2) in e[1]:
                proposed_cc_xerror_array[c1][c2] = not cc_xerror_array[c1][c2]

            # Total number of errors in current system
            currentE = np.count_nonzero(vv_xerror_array) + np.count_nonzero(cc_xerror_array)

            # Total number of high rate and low rate errors in Proposed lattice
            proposedE = np.count_nonzero(proposed_vv_xerror_array) + np.count_nonzero(proposed_cc_xerror_array)
            E = (list(zip(*np.where(proposed_vv_xerror_array))), list(zip(*np.where(proposed_cc_xerror_array))))

            # Metropolis step.
            r = random.random()
            if r <= min([1, self.prob_err1_div_prob_err2(p, proposedE, currentE)]):
                # Non-trivial metropolis step. Attemp to correct the new error
                # configuration. If logical error achieved then the proposed error
                # configuration becomes the new error configuration.

                res = decoder.run_algo_qcode(self.ccode, E, mask, self.logical2)
                if (res == 2 or res == 0): # <----------------- WATCH THIS
                    # update new error
                    vv_xerror_array = np.copy(proposed_vv_xerror_array)
                    cc_xerror_array = np.copy(proposed_cc_xerror_array)
                    currentE = proposedE

            # Start sampling the error configuration after sufficiently long when
            # equilibration has been achieved. This varies on a case to case
            # basis. For the simulations in the article 'The XZZX surface code' we
            # found that after 50000 iterations equilibration had been achieved.
            if counter >= num_iter_cut:
                E = currentE
                if E not in Es:
                    Es[E] = 1
                else:
                    Es[E] += 1
        return Es


    def R_j_given_E(self, n, p1, p2, Es1, Es2):
        """
        Given two probability distributions, P1 and P2, and the error
        configurations sampled under each, it calculates the ratio of logical
        failure rate under the first distribution over the logical failure rate
        under the second distribution. We call this ration R_j.

        Note:

        * R_j is the backbone of the splitting method (PhysRevA.88.062308). R_j
        can be calculated as Rj = C <g(C*Aj)>_P1 / <g(1/(C*Aj))>_P2= C top/bottom,
        where <.>_k is the expectation value of . calculated using the error
        configurations derived from the probability distribution characterised by
        k. When top = bottom, then Rj = C is the required value.

        :param p1: Probability of an error on any one physical qubit under the
        first distribution.
        :type p1: float
        :param p2: Probability of an error on any one physical qubit under the
        second distribution.
        :type p2: float
        :param Es1: Error configurations sampled under the first distribution.
        :param Es2: Error configurations sampled under the second distribution.

        :return: Rj, the ratio of logical failure under distribution characterized
        by p_1 to logical failure under distribution characterized by p_2.
        :rtype: float
        """
        assert p1 < p2
        y1, y2 = [], []  # Where R_j = C y1/y2

        C_vals = list(np.linspace(0.0001, 1, 100))

        for C in C_vals:
            val1 = self.expec_G_CAj(n, p1, p2, C, Es1) / sum(list(Es1.values()))
            val2 = self.expec_G_CAj1(n, p1, p2, C, Es2) / sum(list(Es2.values()))
            y1.append(val1)
            y2.append(val2)

        indx = np.argmin(np.abs(np.array(y1) - np.array(y2)))

        plt.plot(C_vals, y1, label='g_caj')
        plt.plot(C_vals, y2, label='g_caj1')
        plt.plot(C_vals, [1-y for y in y2], label='test')
        plt.legend(loc='lower right')
        C_star = C_vals[indx]  # This is R_j
        return C_star


    def generate_probs(p_start, p_end, d, n):
        """
        Generates an array of physical error probabilities which will be used as
        the pivot points to calculate the ratios of logical failure rate according
        to the splitting method (PhysRevA.88.062308). The sequence of
        probabilities is calculated according to equation (17) in
        PhysRevA.88.062308.

        :param p_start: The small physical error rate where we want to calculate
        the logical error rate at.
        :type p_start: float
        :param p_end: The larger physical error rate where calculating the logical
        error rate is tractable, say, by Monte carlo simulations.
        :type p_end: float
        :param d: Size of dx(d+1) lattice.
        :type d: int
        :param n: Number of physical qubits in the lattice. For the coprime code
        it is d*(d+1).
        :type n: int
        :return: A list of the physical error probabilities to conduct the
        metropolis sampling and splitting method at.
        :rtype: lst of float
        """
        assert p_start < p_end
        probs = []
        p_current = p_start
        probs.append(p_current)
        while True:
            w = max([d / 2, p_current * n])
            p_next = p_current * (2 ** (1 / math.sqrt(w)))
            if p_next > p_end:
                del probs[-1]
                probs.append(p_end)
                break
            probs.append(p_next)
            p_current = p_next
        return probs


    def calculate_and_save_samples(self, probs, N, N_cut):
        """
        Given a physical error rate where we desire to calculate the logical error
        rate at and a physical error rate where we know the logical error rate, it
        partitions the range of probabilities according to eqn. (17) of
        PhysRevA.88.062308, samples error configurations at every partition
        element using the metropolis sampling technique and saves these error
        configurations into files in the working directory.

        :param p_to_find: The small physical error rate where we want to calculate
        the logical error rate at.
        :type p_to_find: float
        :param bias: Bias coefficient
        :type bias: float
        :param d: Size of dx(d+1) lattice.
        :type d: int
        :param known_physical_err: The physical error rate where calculating the
        logical error rate is tractable, say, by Monte carlo simulations.
        :type known_physical_err: float
        :param N: Number of metropolis sampling steps.
        :type N: int
        :param N_cut: The number of simulations ran before we start to sample
        error configurations.
        :type N_cut: int
        """
        # probs = generate_probs(p_to_find, known_physical_err, d, d * (d + 1))
        # probs = [0.01, 0.02]
        for i in range(len(probs)):
            sampling_p = probs[i]
            smp_E = self.sample_E(sampling_p, N, N_cut)
            f = open("hgp_" + "_" + str(probs[0]) + "_" + str(i) + ".txt", "w")
            for E, i in smp_E.items():
                f.write(str(E) + " " + str(i) + "\n")
            f.close()


    def calculate_and_save_single_sample(p_to_find, bias, d, ix_of_sample, known_physical_err, N, N_cut):
        """
        Given a physical error rate where we desire to calculate the logical error
        rate at and a physical error rate where we know the logical error rate, it
        partitions the range of probabilities according to eqn. (17) of
        PhysRevA.88.062308, samples error configurations at a specified partition
        element using the metropolis sampling technique and saves the error
        configurations at the specified physical error rate into files in the
        working directory.

        :param p_to_find: The small physical error rate where we want to calculate
        the logical error rate at.
        :type p_to_find: float
        :param bias: Bias coefficient
        :type bias: float
        :param d: Size of dx(d+1) lattice.
        :type d: int
        :param ix_of_sample: Index of required physical probability where the
        sampling will be conducted on, from the list of the generated partition
        of probabilities.
        :type ix_of_sample: int
        :param known_physical_err: The physical error rate where calculating the
        logical error rate is tractable, say, by Monte carlo simulations.
        :type known_physical_err: float
        :param N: Number of metropolis sampling steps.
        :type N: int
        :param N_cut: The number of simulations ran before we start to sample
        error configurations.
        :type N_cut: int
        """
        probs = generate_probs(p_to_find, known_physical_err, d, d * (d + 1))
        i = ix_of_sample
        sampling_p = probs[i]
        p_h, p_l = sampling_p * bias / (1 + bias), sampling_p / (2 * (1 + bias))
        smp_E = sample_E(N, d, p_h, p_l, sampling_p, bias, N_cut)
        f = open("coprime_" + str(d) + "_" + str(p_to_find) + "_" + str(bias) + "_" + str(i) + ".txt", "w")
        for w, n in smp_E.items():
            H, L = w
            f.write(str(H) + " " + str(L) + " " + str(n) + "\n")
        f.close()


    def read_samples_and_calculate_ratio(self, p_to_find, bias, d, known_physical_err, known_logical_err, N):
        """
        Reads saved error configurations sampled using the Metropolis algorithm
        and uses the splitting method to calculate the logical error rate at the
        physical error rate where a calculation by Monte carlo simulations would
        be intractable.

        :param p_to_find: The small physical error rate where we want to calculate
        the logical error rate at.
        :type p_to_find: float
        :param bias: Bias coefficient
        :type bias: float
        :param d: Size of dx(d+1) lattice.
        :type d: int
        :param known_physical_err: The physical error rate where calculating the
        logical error rate is tractable, say, by Monte carlo simulations.
        :type known_physical_err: float
        :param known_logical_err: The logical error rate at the known physical
        error rate.
        :type known_logical_err: float
        :param N: Number of metropolis sampling steps.
        :type N: int
        :return: The logical error rate at a low physical error rate.
        :rtype: float
        """
        # Generate partition of physical error rates.
        probs = generate_probs(p_to_find, known_physical_err, d, d * (d + 1))

        # Load in saved error configurations.
        sampled_E = []  # List of dictionaries.
        for i in range(len(probs)):
            E_dic = {}
            f = open("coprime_" + str(d) + "_" + str(p_to_find) + "_" + str(bias) + "_" + str(i) + ".txt", "r")
            lines = [line.rstrip().split(" ") for line in f]
            for H, L, n in lines:
                E_dic[(int(H), int(L))] = int(n)
            sampled_E.append(E_dic)

        # Calculate the logical error rate at low physical error rate using the
        # splitting method.
        P = known_logical_err
        for j in range(len(probs) - 1):
            p1, p2 = probs[j], probs[j + 1]
            p1_h, p1_l = p1 * bias / (1 + bias), p1 / (2 * (1 + bias))
            p2_h, p2_l = p2 * bias / (1 + bias), p2 / (2 * (1 + bias))
            P *= self.R_j_given_E(p1_h, p1_l, p1, p2_h, p2_l, p2, sampled_E[j + 1], sampled_E[j], d)

            # R_j_given_E(n, p1, p2, sampled_E[j], sampled_E[j + 1])

        return P


    def test(p_to_find, bias, d, p_known, P_known, N, N_cut, keep_sampled_errs=False):
        """
        Calculates the logical error rate of the coprime XZZX code under biased
        dephasing noise at low physical error rate, potentially where simulations
        via other methods, such as Monte carlo simulations, are not tractable.

        :param p_to_find: The small physical error rate where we want to calculate
        the logical error rate at.
        :type p_to_find: float
        :param bias: Bias coefficient
        :type bias: float
        :param d: Size of dx(d+1) lattice.
        :type d: int
        :param known_physical_err: The physical error rate where calculating the
        logical error rate is tractable, say, by Monte carlo simulations.
        :type known_physical_err: float
        :param known_logical_err: The logical error rate at the known physical
        error rate.
        :type known_logical_err: float
        :param N: Number of metropolis sampling steps.
        :type N: int
        :param N_cut: The number of simulations ran before we start to sample
        error configurations.
        :type N_cut: int
        :param keep_sampled_errs: Whether to keep the sampled error files in the
        current folder or delete them.
        :type keep_sampled_errs: boolean
        """
        start_time = time.time()
        calculate_and_save_samples(p_to_find, bias, d, p_known, N, N_cut)
        calculated_logical_err = read_samples_and_calculate_ratio(p_to_find, bias, d, p_known, P_known, N)

        print([{"code: Rotated coprime XZ " + str(d) + "x" + str(d + 1),
                "decoder: Rotated coprime XZ MWPM",
                "error_model: Biased noise toward dephasing",
                "bias: " + str(bias),
                "error_probability: " + str(p),
                "logical_failure_rate: " + str(calculated_logical_err),
                "known_p: " + str(p_known),
                "known_P: " + str(P_known),
                "n_Metropolis_runs: " + str(N),
                "n_discarded: " + str(N_cut),
                "wall_time: " + str(time.time() - start_time)}])

        if not keep_sampled_errs:
            probs = generate_probs(p_to_find, p_known, d, d * (d + 1))
            for i in range(len(probs)):
                os.remove("coprime_" + str(d) + "_" + str(p_to_find) + "_" + str(bias) + "_" + str(i) + ".txt")


    if __name__ == "__main__":

        try:
            assert len(sys.argv) >= 7
            all_arguments = []
            for i in range(len(sys.argv)):
                if i > 0:
                    argument = sys.argv[i]
                    key, value = argument.split('=')
                    if key == "p":
                        all_arguments.append("p")
                        p = float(value)
                    elif key == "bias":
                        all_arguments.append("bias")
                        bias = float(value)
                    elif key == "L":
                        all_arguments.append("L")
                        d = int(value)
                    elif key == "p_known":
                        all_arguments.append("p_known")
                        p_known = float(value)
                    elif key == "P_known":
                        all_arguments.append("P_known")
                        P_known = float(value)
                    elif key == "N":
                        all_arguments.append("N")
                        N = int(value)
                    elif key == "N_cut":
                        N_cut = int(value)
                    else:
                        assert key == "keep_sampled_errs"
                        keep_sampled_errs = bool(value)

            assert all(item in all_arguments for item in ["p", "bias", "L", "p_known", "P_known", "N"])

        except Exception:
            print("#######################")
            print()
            print("Incorrect input syntax.")
            print(".......................")
            print("Run this program as: ")
            print()
            print("python coprime_low_p_sampler.py p=<target_physical_error_rate>"
                + " bias=<bias_coefficient> L=<size_of_coprime_Lx(L+1)_lattice>"
                + " p_known=<physical_error_rate_where_logical_rate_is_known>"
                + " P_known=<logical_rate_at_p_known>"
                + " N=<tot_number_of_metropolis_steps>")
            print()
            print("<target_physical_error_rate>: float between 0 and 1 (inclusive).")
            print("<bias_coefficient>: float greater than 0.")
            print("<size_of_coprime_Lx(L+1)_lattice>: int positive and odd.")
            print("<physical_error_rate_where_logical_rate_is_known>: float between 0 and 1 (inclusive).")
            print("<logical_rate_at_p_known>: float between 0 and 1 (inclusive).")
            print("<tot_number_of_metropolis_steps>: int greater than zero.")
            print(".......................")
            print("Optional arguments:")
            print("N_cut=<int, greater than 0>: The number of simulations ran"
                + " before starting to sample error configurations. The"
                + " default is 0.")
            print("keep_sampled_errs=<Boolean (True or False)>: Whether to keep "
                + "the sampled error configurations for use in the future or "
                + "not. The default is False.")
            print(".......................")
            print("Example: ")
            print("python coprime_low_p_sampler.py p=1e-3 bias=3 L=7 p_known=0.1 P_known=0.02359 N=5000 N_cut=1000")
            print()
            print("#######################")

            sys.exit()

        try:
            N_Cut
        except NameError:
            N_cut = 0

        try:
            keep_sampled_errs
        except NameError:
            keep_sampled_errs = False

        test(p, bias, d, p_known, P_known, N, N_cut, keep_sampled_errs)