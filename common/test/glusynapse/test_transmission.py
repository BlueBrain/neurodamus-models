import neuron
import numpy as np
import numpy.testing as npt
from scipy import stats
from itertools import product

# Setup NEURON
neuron.h.load_file("stdrun.hoc")
DEBUG = False

# Constants
T = 310
R = 8.314
F = 96485
pm_pca = 1 / 1.38


def Pf(v, cao=1.6, m=120):
    return (4 * cao) / (
        4 * cao + pm_pca * m * (1.0 - np.exp(2 * v * 0.001 * F / (R * T)))
    )


def binf_VCR(v, T):
    """Vargas-Caballero and Robinson, 2003"""
    v05 = -13.0  # mV
    delta = 0.96
    z = 2  # Valence of magnesium ions
    return 1.0 / (1.0 + np.exp(-(v - v05) * z * delta * F / (1000.0 * R * T)))


class TestTransmission(object):
    def setUp(self):
        neuron.h.cvode.active(1)
        neuron.h.cvode.atol(1e-9)
        # Set params
        self.seed = 1234
        # Create model
        self.soma = neuron.h.Section()
        self.soma.insert("pas")
        self.soma.diam = 20.0
        self.soma.L = 20.0
        self.soma.Ra = 200.0
        self.soma.g_pas = 1e-5
        self.soma.cm = 1
        self.syn = neuron.h.GluSynapse(0.5, sec=self.soma)
        self.syn.setRNG(self.seed, self.seed)
        self.syn.verbose = DEBUG

    def finalizemodel(self):
        # Initialize neuron
        neuron.h.v_init = -70.0
        neuron.h.celsius = 34.0
        neuron.h.stdinit()

    def test_ampar_peakcond(self):
        """Test AMPAR peak conductance (location and magnitude).
        Peak conductance must be equal to gmax. Furthermore, time to peak must
        be prescribed by tau_r_AMPA and tau_d_AMPA."""
        # Activate synapse after 100 ms
        tspike = 100.0
        vecStim = neuron.h.VecStim()
        vec = neuron.h.Vector([tspike])
        vecStim.play(vec)
        netCon = neuron.h.NetCon(vecStim, self.syn)
        netCon.weight[0] = 1
        netCon.delay = 2  # ms
        self.syn.Use = 1
        # Disable AMPAR plasticity
        self.syn.theta_d_GB = np.inf
        self.syn.theta_p_GB = np.inf
        # Test 1 nS and 15 nS
        for gmax in (1.0, 15.0):
            self.syn.gmax0_AMPA = gmax
            # Set recordings
            i = neuron.h.Vector()
            i.record(self.syn._ref_i)
            t = neuron.h.Vector()
            t.record(neuron.h._ref_t)
            g = neuron.h.Vector()
            g.record(self.syn._ref_g_AMPA)
            # Run simulation
            self.finalizemodel()
            neuron.run(200.0)
            if DEBUG:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                ax.plot(t, g)
                plt.show()
            # Check peak conductance
            npt.assert_almost_equal(g.max(), 1e-3 * gmax, decimal=5)
            # Check peak location
            tpeak = t[np.argmax(g)]
            tpeak_expected = (
                tspike
                + (neuron.h.tau_r_AMPA_GluSynapse * self.syn.tau_d_AMPA)
                / (self.syn.tau_d_AMPA - neuron.h.tau_r_AMPA_GluSynapse)
                * np.log(self.syn.tau_d_AMPA / neuron.h.tau_r_AMPA_GluSynapse)
                + netCon.delay
            )
            npt.assert_almost_equal(tpeak, tpeak_expected, decimal=1)

    def test_nmdar_peakcond(self):
        """Test NMDAR peak conductance (location and magnitude), no magnesium.
        Peak conductance must be equal to gmax. Furthermore, time to peak must
        be prescribed by tau_r_NMDA and tau_d_NMDA."""
        # Activate synapse after 100 ms
        tspike = 100.0
        vecStim = neuron.h.VecStim()
        vec = neuron.h.Vector([tspike])
        vecStim.play(vec)
        netCon = neuron.h.NetCon(vecStim, self.syn)
        netCon.weight[0] = 1
        netCon.delay = 2  # ms
        self.syn.Use = 1
        # Magnesium free condition
        neuron.h.mg_GluSynapse = 0
        # Test 1 nS and 15 nS
        for gmax in (1.0, 15.0):
            self.syn.gmax_NMDA = gmax
            # Set recordings
            i = neuron.h.Vector()
            i.record(self.syn._ref_i)
            t = neuron.h.Vector()
            t.record(neuron.h._ref_t)
            g = neuron.h.Vector()
            g.record(self.syn._ref_g_NMDA)
            # Run simulation
            self.finalizemodel()
            neuron.run(200.0)
            if DEBUG:
                import matplotlib.pyplot as plt

                fig, ax = plt.subplots()
                ax.plot(t, g)
                plt.show()
            # Check peak conductance
            npt.assert_almost_equal(g.max(), 1e-3 * gmax, decimal=5)
            # Check peak location
            tpeak = t[np.argmax(g)]
            tpeak_expected = (
                tspike
                + (neuron.h.tau_r_NMDA_GluSynapse * neuron.h.tau_d_NMDA_GluSynapse)
                / (neuron.h.tau_d_NMDA_GluSynapse - neuron.h.tau_r_NMDA_GluSynapse)
                * np.log(neuron.h.tau_d_NMDA_GluSynapse / neuron.h.tau_r_NMDA_GluSynapse)
                + netCon.delay
            )
            npt.assert_almost_equal(tpeak, tpeak_expected, decimal=1)

    def release(self, u, Nrrp):
        # Activate synapse after 100 ms
        tspike = 100.0
        vecStim = neuron.h.VecStim()
        vec = neuron.h.Vector([tspike])
        vecStim.play(vec)
        netCon = neuron.h.NetCon(vecStim, self.syn)
        netCon.weight[0] = 1
        netCon.delay = 2  # ms
        # Disable AMPAR and Use plasticity
        self.syn.theta_d_GB = -1
        self.syn.theta_p_GB = -1
        # Set release probability
        self.syn.Use = u
        self.syn.Nrrp = Nrrp
        # Test many trials
        n_trials = 10000
        rel_events = np.zeros(Nrrp + 1)
        for trial in range(n_trials):
            # Set recordings
            i = neuron.h.Vector()
            i.record(self.syn._ref_i)
            t = neuron.h.Vector()
            t.record(neuron.h._ref_t)
            g = neuron.h.Vector()
            g.record(self.syn._ref_g_AMPA)
            # Run simulation
            self.finalizemodel()
            neuron.run(200.0)
            estimated_release = int(
                np.round(Nrrp * g.max() / (1e-3 * self.syn.gmax_AMPA))
            )
            rel_events[estimated_release] += 1
            if DEBUG:
                import matplotlib.pyplot as plt

                fig, ax = plt.subplots()
                ax.plot(t, g)
                plt.show()
        actual_prob = rel_events / float(n_trials)
        target_prob = stats.binom.pmf(np.array(range(Nrrp + 1)), Nrrp, u)
        npt.assert_array_almost_equal(target_prob, actual_prob, decimal=2)

    def test_release(self):
        """Test stochastic release.
        A single vesicle must be released with probability Use, assuming a long
        recovery time. Multiple vesicle release follows binomial statistics."""
        for Nrrp, u in product((1, 3), (0.2, 0.5, 0.8)):
            yield self.release, u, Nrrp

    def binomial(self, u, Nrrp):
        self.finalizemodel()
        n_trials = 1000000
        rel_events = np.zeros(Nrrp + 1)
        for trial in range(n_trials):
            res = int(self.syn.brand(Nrrp, u))
            rel_events[res] += 1
        actual_prob = rel_events / float(n_trials)
        target_prob = stats.binom.pmf(np.array(range(Nrrp + 1)), Nrrp, u)
        npt.assert_array_almost_equal(target_prob, actual_prob, decimal=3)

    def test_binomial(self):
        for Nrrp, u in product((1, 3), (0.2, 0.5, 0.8)):
            yield self.binomial, u, Nrrp

    def test_nmdar_calciumcurrent(self):
        """Test NMDAR-mediated calcium current against GHK estimate.
        As shown by Schneggenburger et al. 1993, the fractional calcium current
        through the NMDAR follows quite closely the GHK flux equation.
        """
        # Store current value of tau_d_NMDA and set new value to (almost) infinity
        tau_d_NMDA = neuron.h.tau_d_NMDA_GluSynapse
        neuron.h.tau_d_NMDA_GluSynapse = 1e12
        # Store current value of cao_CR, updated during the test
        cao_CR = neuron.h.cao_CR_GluSynapse
        # Activate synapse after 100 ms
        tspike = 100.0
        vecStim = neuron.h.VecStim()
        vec = neuron.h.Vector([tspike])
        vecStim.play(vec)
        netCon = neuron.h.NetCon(vecStim, self.syn)
        netCon.weight[0] = 1
        netCon.delay = 2  # ms
        self.syn.Use = 1
        # Disable AMPAR and Use plasticity
        self.syn.theta_d_GB = -1
        self.syn.theta_p_GB = -1
        # Extracellular calcium and voltage levels to test
        caovec = [1.2, 1.6, 2.0]
        vstepvec = np.linspace(-85.0, 35.0, 9)
        # Run tests
        for cao in caovec:
            if DEBUG:
                import matplotlib.pyplot as plt

                fig, ax = plt.subplots()
            for vstep in vstepvec:
                # Add voltage clamp
                vclamp = neuron.h.SEClamp(0.5, sec=self.soma)
                vclamp.amp1 = -70.0
                vclamp.amp2 = vstep
                vclamp.amp3 = -70.0
                vclamp.dur1 = 200
                vclamp.dur2 = 600
                vclamp.dur3 = 200
                # Set calcium level
                neuron.h.cao_CR_GluSynapse = cao
                # Set recordings
                t = neuron.h.Vector()
                v = neuron.h.Vector()
                g_NMDA = neuron.h.Vector()
                ica_NMDA = neuron.h.Vector()
                t.record(neuron.h._ref_t)
                v.record(self.soma(0.5)._ref_v)
                g_NMDA.record(self.syn._ref_g_NMDA)
                ica_NMDA.record(self.syn._ref_ica_NMDA)
                # Run sim
                self.finalizemodel()
                neuron.run(1000.0)
                # Compute ica/i ratio
                i_NMDA = np.array(g_NMDA) * (np.array(v) - neuron.h.E_NMDA_GluSynapse)
                tidx = np.searchsorted(t, 750.0)
                iratio = ica_NMDA[tidx] / i_NMDA[tidx]
                # GHK ica/i ratio
                iratio_ghk = Pf(vstep, cao=cao)
                # Check error
                if vstep < -20 or vstep > 20:
                    # Small error tolerance away from 0 mV
                    assert abs(iratio_ghk - iratio) < 0.06, "%f %f" % (
                        iratio_ghk,
                        iratio,
                    )
                else:
                    # Larger tolerance around 0 mV
                    assert abs(iratio_ghk - iratio) < 0.25, "%f %f" % (
                        iratio_ghk,
                        iratio,
                    )
                if DEBUG:
                    ax.plot(vstep, iratio, "bo")
                    ax.plot(vstep, iratio_ghk, "rs")
            if DEBUG:
                plt.show()
        # Restore tau_d_NMDA ad cao_CR
        neuron.h.tau_d_NMDA_GluSynapse = tau_d_NMDA
        neuron.h.cao_CR_GluSynapse = cao_CR

    def test_nmdar_magnesium_block(self):
        """Test fraction of unblocked NMDAR conductance.
        The model must reproduce the results in Fig. 1B of Vargas-Caballero and
        Robinson 2003.
        """
        # Store current value of tau_d_NMDA and set new value to (almost) infinity
        tau_d_NMDA = neuron.h.tau_d_NMDA_GluSynapse
        neuron.h.tau_d_NMDA_GluSynapse = 1e12
        # Activate synapse after 100 ms
        tspike = 100.0
        vecStim = neuron.h.VecStim()
        vec = neuron.h.Vector([tspike])
        vecStim.play(vec)
        netCon = neuron.h.NetCon(vecStim, self.syn)
        netCon.weight[0] = 1
        netCon.delay = 2  # ms
        self.syn.Use = 1
        # Disable AMPAR and Use plasticity
        self.syn.theta_d_GB = -1
        self.syn.theta_p_GB = -1
        # Voltage steps to test
        vvec = np.linspace(-70, 40, 15)
        scond = []
        for vstep in vvec:
            # Add voltage clamp
            vclamp = neuron.h.SEClamp(0.5, sec=self.soma)
            vclamp.amp1 = -70.0
            vclamp.amp2 = vstep
            vclamp.amp3 = -70.0
            vclamp.dur1 = 200
            vclamp.dur2 = 600
            vclamp.dur3 = 200
            # Set recordings
            t = neuron.h.Vector()
            v = neuron.h.Vector()
            g_NMDA = neuron.h.Vector()
            t.record(neuron.h._ref_t)
            v.record(self.soma(0.5)._ref_v)
            g_NMDA.record(self.syn._ref_g_NMDA)
            # Run sim
            self.finalizemodel()
            neuron.run(1000.0)
            # Get stationary conductance
            tidx = np.searchsorted(t, 750.0)
            scond.append(
                float(g_NMDA[tidx]) / (1e-3 * self.syn.gmax_NMDA)
            )  # Gmax NMDA needs to be manually computed
        scond = np.array(scond)
        # Compute target
        T_body = 309.15
        target_scond = binf_VCR(vvec, T_body)
        if DEBUG:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots()
            ax.plot(vvec, scond, label="NMODL")
            ax.plot(vvec, target_scond, label="Boltzmann distribution")
            ax.legend(loc=0)
            plt.show()
        # Test deviation from target
        npt.assert_array_almost_equal(scond, target_scond, decimal=3)
        # Restore tau_d_NMDA ad cao_CR
        neuron.h.tau_d_NMDA_GluSynapse = tau_d_NMDA
