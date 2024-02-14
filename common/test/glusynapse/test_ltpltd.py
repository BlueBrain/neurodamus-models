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
        # Create model
        seed = 1234
        self.soma, self.dend = get_ballstick()
        self.syn = neuron.h.GluSynapse(0.5, sec=self.dend)
        self.syn.setRNG(seed, seed + 1, seed + 2)

    def finalizemodel(self):
        # Initialize neuron
        neuron.h.v_init = -70.0
        neuron.h.celsius = 34.0
        neuron.h.stdinit()

    def steady_state_rho(self, rho0, gamma_d, theta_x):
        # Fixed parameters
        gamma_p = 100.0
        # Target rho and tau
        if theta_x == 0:
            target_rho = gamma_p / (gamma_p + gamma_d)
            tau = 10.0
        elif theta_x == np.inf and rho0 < neuron.h.rho_star_GB_GluSynapse:
            target_rho = 0
            tau = 0.05
        elif theta_x == np.inf and rho0 == neuron.h.rho_star_GB_GluSynapse:
            target_rho = neuron.h.rho_star_GB_GluSynapse
            tau = 0.05
        elif theta_x == np.inf and rho0 > neuron.h.rho_star_GB_GluSynapse:
            target_rho = 1
            tau = 0.05
        else:
            raise NotImplementedError()
        # All other parameters
        self.syn.theta_d_GB = theta_x
        self.syn.theta_p_GB = theta_x
        self.syn.rho0_GB = rho0
        neuron.h.tau_ind_GB_GluSynapse = tau
        neuron.h.gamma_d_GB_GluSynapse = gamma_d
        neuron.h.gamma_p_GB_GluSynapse = gamma_p
        # Recordings
        rho = neuron.h.Vector()
        t = neuron.h.Vector()
        rho.record(self.syn._ref_rho_GB)
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

        for Use, gamma_d, theta in product(
            Use_vec, gamma_d_vec, theta_vec
        ):
            yield self.steady_state_rho, Use, gamma_d, theta

    def use_convergence(self, rho_GB, Use):
        # Parameters
        self.syn.theta_d_GB = -1
        self.syn.theta_p_GB = -1
        self.syn.Use_d = 0.1
        self.syn.Use_p = 0.9
        self.syn.Use = Use
        neuron.h.tau_exp_GB_GluSynapse = 0.1
        # Recordings
        use = neuron.h.Vector()
        t = neuron.h.Vector()
        use.record(self.syn._ref_Use_GB)
        t.record(neuron.h._ref_t)
        self.finalizemodel()
        # Simulate rho change
        self.syn.rho_GB = rho_GB
        # Run
        neuron.run(2000.0)
        if DEBUG:
            plt.figure()
            plt.plot(t, use)
            plt.title(f"rho_GR {rho_GB}")
            plt.show()
        # Test
        usenp = use.as_numpy()
        use_target = self.syn.Use_d + rho_GB * (self.syn.Use_p - self.syn.Use_d)
        npt.assert_almost_equal(usenp[-1], use_target, decimal=2)

    def gmax_convergence(self, rho_GB, gmax0_AMPA):
        # Parameters
        self.syn.theta_d_GB = -1
        self.syn.theta_p_GB = -1
        self.syn.gmax_d_AMPA = 1.0
        self.syn.gmax_p_AMPA = 3.0
        self.syn.gmax0_AMPA = gmax0_AMPA
        neuron.h.tau_exp_GB_GluSynapse = 0.1
        # Recordings
        g = neuron.h.Vector()
        t = neuron.h.Vector()
        g.record(self.syn._ref_gmax_AMPA)
        t.record(neuron.h._ref_t)
        self.finalizemodel()
        # Simulate rho change
        self.syn.rho_GB = rho_GB
        # Run
        neuron.run(2000.0)
        if DEBUG:
            plt.figure()
            plt.plot(t, g)
            plt.title(f"rho_GR {rho_GB}")
            plt.show()
        # Test
        gnp = g.as_numpy()
        gmax_target = self.syn.gmax_d_AMPA + rho_GB * (self.syn.gmax_p_AMPA - self.syn.gmax_d_AMPA)
        npt.assert_almost_equal(gnp[-1], gmax_target, decimal=2)

    def test_convergence(self):
        rho_GB_vec = [0.0, 0.15, 0.50, 0.85, 1.00]
        Use_vec = [0.0, 0.05, 0.1, 0.15, 0.5, 0.85, 0.95, 1.00]
        gmax0_AMPA_vec = [0.0, 0.25, 1.00, 1.25, 2.0, 2.75, 3.00, 3.75, 4.00]

        for rho_GB, Use in product(
            rho_GB_vec, Use_vec
        ):
            yield self.use_convergence, rho_GB, Use

        for rho_GB, gmax0_AMPA in product(
            rho_GB_vec, gmax0_AMPA_vec
        ):
            yield self.gmax_convergence, rho_GB, gmax0_AMPA
