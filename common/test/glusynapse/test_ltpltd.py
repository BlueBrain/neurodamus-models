import neuron
import numpy as np
import numpy.testing as npt
from itertools import product

DEBUG = False
neuron.h.load_file("stdrun.hoc")

if DEBUG:
    import matplotlib.pyplot as plt


def get_ballstick():
    # Cerate a simple ball and stick model
    soma = neuron.h.Section()
    soma.insert("pas")
    soma.diam = 20.0
    soma.L = 20.0
    soma.Ra = 200.0
    soma.g_pas = 1e-5
    soma.cm = 1
    dend = neuron.h.Section()
    dend.insert("pas")
    dend.diam = 2.0
    dend.L = 100.0
    dend.Ra = 200.0
    dend.g_pas = 1e-5
    dend.cm = 1
    dend.connect(soma)
    return soma, dend


class TestLTPLTD(object):
    def setUp(self):
        # neuron.h.cvode.active(1)
        # neuron.h.cvode.atol(0.0001)
        # Create model
        seed = 1234
        self.soma, self.dend = get_ballstick()
        self.syn = neuron.h.GluSynapse_TM(0.5, sec=self.dend)
        self.syn.setRNG(seed, seed + 1, seed + 2)

    def finalizemodel(self):
        # Initialize neuron
        neuron.h.v_init = -70.0
        neuron.h.celsius = 34.0
        neuron.h.stdinit()

    def steady_state_rho(self, rho0, gamma_d, theta_x, locus):
        # Fixed parameters
        gamma_p = 100.0
        # Target rho and tau
        if theta_x == 0:
            target_rho = gamma_p / (gamma_p + gamma_d)
            tau = 10.0
        elif theta_x == np.inf and rho0 < neuron.h.rho_star_GB_GluSynapse_TM:
            target_rho = 0
            tau = 0.05
        elif theta_x == np.inf and rho0 == neuron.h.rho_star_GB_GluSynapse_TM:
            target_rho = neuron.h.rho_star_GB_GluSynapse_TM
            tau = 0.05
        elif theta_x == np.inf and rho0 > neuron.h.rho_star_GB_GluSynapse_TM:
            target_rho = 1
            tau = 0.05
        else:
            raise NotImplementedError()
        # All other parameters
        setattr(self.syn, "theta_d_%s_GB" % locus, theta_x)
        setattr(self.syn, "theta_p_%s_GB" % locus, theta_x)
        setattr(self.syn, "rho0_%s_GB" % locus, rho0)
        neuron.h.tau_ind_GB_GluSynapse_TM = tau
        setattr(neuron.h, "gamma_d_%s_GB_GluSynapse_TM" % locus, gamma_d)
        setattr(neuron.h, "gamma_p_%s_GB_GluSynapse_TM" % locus, gamma_p)
        # Recordings
        rho = neuron.h.Vector()
        t = neuron.h.Vector()
        rho.record(getattr(self.syn, "_ref_rho_%s_GB" % locus))
        t.record(neuron.h._ref_t)
        self.finalizemodel()
        # Run
        neuron.run(2000.0)
        if DEBUG:
            plt.figure()
            plt.plot(t, rho)
            plt.show()
        # Test
        rhonp = rho.as_numpy()
        tnp = t.as_numpy()
        npt.assert_almost_equal(rhonp[-1], target_rho, decimal=3)
        npt.assert_almost_equal(
            np.std(rhonp[np.searchsorted(tnp, 1500): -1]), 0, decimal=3
        )

    def test_generator_steady_state_rho(self):
        gamma_d_vec = [15, 100, 150]
        Use_vec = [0.1, 0.5, 0.9]
        theta_vec = [0, np.inf]
        locus_vec = ["pre", "post"]

        for Use, gamma_d, theta, locus in product(
            Use_vec, gamma_d_vec, theta_vec, locus_vec
        ):
            yield self.steady_state_rho, Use, gamma_d, theta, locus

    def test_use_convergence(self):
        # Parameters
        self.syn.theta_d_pre_GB = np.inf
        self.syn.theta_p_pre_GB = np.inf
        self.syn.Use_d_TM = 0.1
        self.syn.Use_p_TM = 0.9
        self.syn.Use0_TM = 0.1
        neuron.h.tau_exp_GB_GluSynapse_TM = 0.1
        # Recordings
        use = neuron.h.Vector()
        t = neuron.h.Vector()
        use.record(self.syn._ref_Use_TM)
        t.record(neuron.h._ref_t)
        self.finalizemodel()
        # Simulate rho change
        self.syn.rho_pre_GB = 1.0
        # Run
        neuron.run(2000.0)
        if DEBUG:
            plt.figure()
            plt.plot(t, use)
            plt.show()
        # Test
        usenp = use.as_numpy()
        npt.assert_almost_equal(usenp[-1], 0.9, decimal=5)

    def test_gmax_convergence(self):
        # Parameters
        self.syn.theta_d_post_GB = np.inf
        self.syn.theta_p_post_GB = np.inf
        self.syn.gmax_d_AMPA = 1.0
        self.syn.gmax_p_AMPA = 3.0
        self.syn.gmax0_AMPA = 0.1
        neuron.h.tau_exp_GB_GluSynapse_TM = 0.1
        # Recordings
        g = neuron.h.Vector()
        t = neuron.h.Vector()
        g.record(self.syn._ref_gmax_AMPA)
        t.record(neuron.h._ref_t)
        self.finalizemodel()
        # Simulate rho change
        self.syn.rho_post_GB = 1.0
        # Run
        neuron.run(2000.0)
        if DEBUG:
            plt.figure()
            plt.plot(t, g)
            plt.show()
        # Test
        gnp = g.as_numpy()
        npt.assert_almost_equal(gnp[-1], 3, decimal=5)
