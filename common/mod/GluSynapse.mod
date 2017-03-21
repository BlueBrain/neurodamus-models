COMMENT
/**
 * @file GluSynapse.mod
 * @brief Probabilistic synapse with short- and long-term plasticity
 * @author king, chindemi, rossert
 * @date 2017-03-16
 * @version 0.3.4
 * @remark Copyright © BBP/EPFL 2005-2017; All rights reserved.
           Do not distribute without further notice.
 */
ENDCOMMENT


TITLE Glutamatergic synapse


COMMENT
Glutamatergic synapse model featuring:
1) AMPA receptor with a dual-exponential conductance profile.
2) NMDA receptor  with a dual-exponential conductance profile and magnesium
   block as described in Jahr and Stevens 1990.
3) Tsodyks-Markram presynaptic short-term plasticity as Rahmon et al. 201x.
   Implementation based on the work of Eilif Muller, Michael Reimann and
   Srikanth Ramaswamy (Blue Brain Project, August 2011), who introduced the
   2-state Markov model of vesicle release. The new model is an extension of
   Fuhrmann et al. 2002, motivated by the following constraints:
        a) No consumption on failure
        b) No release until recovery
        c) Same ensemble averaged trace as canonical Tsodyks-Markram using same
           parameters determined from experiment.
   For a pre-synaptic spike or external spontaneous release trigger event, the
   synapse will only release if it is in the recovered state, and with
   probability u (which follows facilitation dynamics). If it releases, it will
   transition to the unrecovered state. Recovery is as a Poisson process with
   rate 1/Dep.
   John Rahmon and Giuseppe Chindemi introduced multi-vesicular release as an
   extension of the 2-state Markov model of vesicle release described above
   (Blue Brain Project, February 2017).
4) NMDAR-mediated calcium current. Fractional calcium current Pf_NMDA from
   Schneggenburger et al. 1993. Fractional NMDAR conductance treated as a
   calcium-only permeable channel with Erev = 40 mV independent of extracellular
   calcium concentration (see Jahr and Stevens 1993). Implemented by Christian
   Rossert and Giuseppe Chindemi (Blue Brain Project, 2016).
5) Spine
6) VDCC
7) Postsynaptic calcium dynamics.
8) Long-term synaptic plasticity. Calcium-based STDP model based on Graupner and
   Brunel 2012.
Model implementation, optimization and simulation curated by James King (Blue
Brain Project, 2017).
ENDCOMMENT


NEURON {
    THREADSAFE
    POINT_PROCESS GluSynapse

    : AMPA Receptor
    RANGE tau_r_AMPA, tau_d_AMPA, E_AMPA, gmax_AMPA
    RANGE factor_AMPA
    RANGE g_AMPA, i_AMPA        : Could be converted to LOCAL (performance)

    : NMDA Receptor
    RANGE tau_r_NMDA, tau_d_NMDA, E_NMDA
    RANGE factor_NMDA
    RANGE gmax_NMDA             : Could be converted to LOCAL (performance)
    RANGE g_NMDA, i_NMDA        : Could be converted to LOCAL (performance)

    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    RANGE Use, Dep, Fac, Nrrp
    RANGE u, Psurv              : Could be converted to LOCAL (performance)
    RANGE tsyn, unoccupied, occupied
    BBCOREPOINTER rng_rel

    : NMDAR-mediated calcium current
    RANGE Pf_NMDA, ica_NMDA

    : Spine
    RANGE volume_CR, area_CR

    : VDCC (R-type)
    RANGE gca_bar_VDCC, gca_bar_abs_VDCC, ica_VDCC, ljp_VDCC
    RANGE gca_VDCC              : Could be converted to LOCAL (performance)
    RANGE vhm_VDCC, km_VDCC, mtau_VDCC, minf_VDCC
    RANGE vhh_VDCC, kh_VDCC, htau_VDCC, hinf_VDCC
    RANGE tm_VDCC, th_VDCC

    : Postsynaptic Ca2+ dynamics
    RANGE gamma_ca_CR, tau_ca_CR, min_ca_CR, cao_CR

    : Long-term synaptic plasticity
    RANGE tau_GB, theta_d_GB, theta_p_GB, gamma_d_GB, gamma_p_GB
    RANGE rho_star_GB, rho0_GB
    RANGE enable_GB, depress_GB, potentiate_GB
    RANGE tau_Use_GB, Use_d_GB, Use_p_GB

    : Basic Synapse and legacy
    RANGE NMDA_ratio, w, mg
    RANGE g                     : Could be converted to LOCAL (performance)
    RANGE synapseID, verbose
    NONSPECIFIC_CURRENT i
}


UNITS {
    (nA)    = (nanoamp)
    (mV)    = (millivolt)
    (uS)    = (microsiemens)
    (nS)    = (nanosiemens)
    (pS)    = (picosiemens)
    (umho)  = (micromho)
    (um)    = (micrometers)
    (mM)    = (milli/liter)
    (uM)    = (micro/liter)
    FARADAY = (faraday) (coulomb)
    PI      = (pi)      (1)
    R       = (k-mole)  (joule/degC)
}


PARAMETER {
    celsius             (degC)

    : AMPA Receptor
    tau_r_AMPA  = 0.2   (ms)    : Tau rise, dual-exponential conductance profile
    tau_d_AMPA  = 1.7   (ms)    : Tau decay, IMPORTANT: tau_r < tau_d
    E_AMPA      = 0     (mV)    : Reversal potential
    gmax_AMPA   = 1.0   (nS)    : Maximal conductance

    : NMDA Receptor
    tau_r_NMDA   = 0.29 (ms)    : Tau rise, dual-exponential conductance profile
    tau_d_NMDA   = 43   (ms)    : Tau decay, IMPORTANT: tau_r < tau_d
    E_NMDA      = 0     (mV)    : Reversal potential

    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    Use         = 1.0   (1)     : Utilization of synaptic efficacy
    Dep         = 100   (ms)    : Relaxation time constant from depression
    Fac         = 10    (ms)    : Relaxation time constant from facilitation
    Nrrp        = 1     (1)     : Number of release sites for given contact

    : Spine
    volume_CR   = 0.087 (um3)   : From spine data by Ruth Benavides-Piccione
                                : (unpublished), value overwritten at runtime

    : VDCC (R-type)
    gca_bar_VDCC = 0.0744   (nS/um2)
    ljp_VDCC     = 0        (mV)
    vhm_VDCC     = -5.9     (mV)    : v 1/2 for act, Magee and Johnston 1995 (corrected for m*m)
    km_VDCC      = 9.5      (mV)    : act slope, Magee and Johnston 1995 (corrected for m*m)
    vhh_VDCC     = -39      (mV)    : v 1/2 for inact, Magee and Johnston 1995
    kh_VDCC      = -9.2     (mV)    : inact, Magee and Johnston 1995
    tm_VDCC      = 1        (ms)    : max time constant (guess)
    th_VDCC      = 27       (ms)    : max time constant 100*0.27

    : Postsynaptic Ca2+ dynamics
    gamma_ca_CR  = 0.04     (1)     : Percent of free calcium (not buffered), Sabatini et al 2002: kappa_e = 24+-11 (also 14 (2-31) or 22 (18-33))
    tau_ca_CR    = 12       (ms)    : Rate of removal of calcium, Sabatini et al 2002: 14ms (12-20ms)
    min_ca_CR    = 70e-6    (mM)    : Sabatini et al 2002: 70+-29 nM, per AP: 1.1 (0.6-8.2) uM = 1100 e-6 mM = 1100 nM
    cao_CR       = 2.0      (mM)    : Extracellular calcium concentration in slices

    : Long-term synaptic plasticity
    tau_GB       = 100      (s)
    theta_d_GB   = 0.006    (mM)
    theta_p_GB   = 0.001    (mM)
    gamma_d_GB   = 100      (1)
    gamma_p_GB   = 450      (1)
    rho_star_GB  = 0.5      (1)
    rho0_GB      = 0        (1)
    enable_GB    = 0        (1)
    tau_Use_GB   = 100      (s)
    Use_d_GB     = 0.2      (1)
    Use_p_GB     = 0.8      (1)

    : Basic Synapse and legacy
    NMDA_ratio  = 0.71  (1)     : In this model gmax_NMDA = gmax_AMPA*ratio_NMDA
    mg          = 1     (mM)    : Extracellular magnesium concentration
    synapseID   = 0
    verbose     = 0
}


COMMENT
The Verbatim block is needed to generate random nos. from a uniform distribution
between 0 and 1 for comparison with Pr to decide whether to activate the synapse
or not.
ENDCOMMENT
VERBATIM
// for MCellRan4
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// for random123
#include "nrnran123.h"

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM


ASSIGNED {
    : AMPA Receptor
    factor_AMPA (1)
    g_AMPA      (uS)
    i_AMPA      (nA)

    : NMDA Receptor
    factor_NMDA (1)
    gmax_NMDA   (nS)
    g_NMDA      (uS)
    i_NMDA      (nA)

    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    u           (1)     : Running release probability
    tsyn        (ms)    : Time of the last presynaptic spike
    Psurv       (1)     : Survival prob. of unrecovered state
    unoccupied  (1)     : Number of unoccupied release sites
    occupied    (1)     : Number of occupied release sites
    rng_rel             : Random Number Generator
    usingR123           : TEMPORARY until mcellran4 completely deprecated

    : NMDAR-mediated calcium current
    Pf_NMDA     (1)     : Fractional NMDAR calcium current
    ica_NMDA    (nA)

    : Spine
    area_CR     (um2)

    : VDCC (R-type)
    gca_VDCC            (uS)
    gca_bar_abs_VDCC    (nS)
    minf_VDCC           (1)
    hinf_VDCC           (1)
    mtau_VDCC           (ms)
    htau_VDCC           (ms)
    ica_VDCC            (nA)

    : Long-term synaptic plasticity
    depress_GB          (1)
    potentiate_GB       (1)

    : Basic Synapse and legacy
    w           (1)
    v           (mV)
    g           (uS)
    i           (nA)
}

STATE {
    : AMPA Receptor
    A_AMPA      (1)             : Decays with conductance tau_r_AMPA
    B_AMPA      (1)             : Decays with conductance tau_d_AMPA

    : NMDA Receptor
    A_NMDA      (1)             : Decays with conductance tau_r_NMDA
    B_NMDA      (1)             : Decays with conductance tau_d_NMDA

    : VDCC (R-type)
    m_VDCC      (1)
    h_VDCC      (1)

    : Postsynaptic Ca2+ dynamics
    cai_CR      (mM)    <1e-3>  : Intracellular calcium concentration

    : Long-term synaptic plasticity
    Rho_GB      (1)
    Use_GB      (1)
}

INITIAL{
    LOCAL tp_AMPA, tp_NMDA

    : AMPA Receptor
    A_AMPA      = 0
    B_AMPA      = 0
    : Time to peak of the conductance
    tp_AMPA     = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA)
    : Normalization factor - so that when t = tp_AMPA, g_AMPA = gmax_AMPA
    factor_AMPA = 1 / (-exp(-tp_AMPA/tau_r_AMPA)+exp(-tp_AMPA/tau_d_AMPA))

    : NMDA Receptor
    A_NMDA      = 0
    B_NMDA      = 0
    : Time to peak of the conductance
    tp_NMDA     = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA)
    : Normalization factor - so that when t = tp_NMDA, g_NMDA = gmax_NMDA
    factor_NMDA = 1 / (-exp(-tp_NMDA/tau_r_NMDA)+exp(-tp_NMDA/tau_d_NMDA))
    gmax_NMDA   = gmax_AMPA*NMDA_ratio

    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    tsyn        = 0
    Psurv       = 0
    u           = 0
    unoccupied  = 0
    occupied    = Nrrp

    : NMDAR-mediated calcium current
    Pf_NMDA     = (4*cao_CR) / (4*cao_CR + (1/1.38) * 120 (mM)) * 0.6

    : Spine
    UNITSOFF
    area_CR     = 4*PI*(3/4*volume_CR*1/PI)^(2/3)   : Assuming sphere for spine head
    UNITSON

    : VDCC (R-type)
    gca_bar_abs_VDCC    = gca_bar_VDCC * area_CR
    mtau_VDCC           = tm_VDCC
    htau_VDCC           = th_VDCC

    : Postsynaptic Ca2+ dynamics
    cai_CR      = min_ca_CR

    : Long-term synaptic plasticity
    Rho_GB          = rho0_GB
    Use_GB          = Use
    depress_GB      = 0
    potentiate_GB   = 0

    : Basic Synapse and legacy
    w           = 0

    : Initialize WATCH
    net_send(0, 1)
}


BREAKPOINT {
    LOCAL Eca_syn, mggate
    SOLVE state METHOD derivimplicit

    : AMPA Receptor
    g_AMPA = w*(1e-3)*gmax_AMPA*(B_AMPA-A_AMPA)
    i_AMPA = g_AMPA*(v-E_AMPA)

    : NMDA Receptor
    mggate = 1 / (1 + exp(0.062 (/mV) * -(v)) * (mg / 3.57 (mM)))
    g_NMDA = w*(1e-3)*gmax_NMDA*mggate*(B_NMDA-A_NMDA)
    i_NMDA = g_NMDA*(v-E_NMDA)

    : NMDAR-mediated calcium current
    ica_NMDA = Pf_NMDA*g_NMDA*(v-40.0)

    : VDCC (R-type)
    Eca_syn = nernst(cai_CR, cao_CR, 2)     : Ca reversal potential
    gca_VDCC = (1e-3) * gca_bar_abs_VDCC * m_VDCC * m_VDCC * h_VDCC
    ica_VDCC = gca_VDCC*(v-Eca_syn)
    minf_VDCC = 1 / (1 + exp(((vhm_VDCC - ljp_VDCC) - v) / km_VDCC))
    hinf_VDCC = 1 / (1 + exp(((vhh_VDCC - ljp_VDCC) - v) / kh_VDCC))

    : Update total g and inject current
    g = g_AMPA + g_NMDA
    i = i_AMPA + i_NMDA + ica_VDCC
}


DERIVATIVE state {
    : AMPA Receptor
    A_AMPA'     = -A_AMPA/tau_r_AMPA
    B_AMPA'     = -B_AMPA/tau_d_AMPA

    : NMDA Receptor
    A_NMDA'     = -A_NMDA/tau_r_NMDA
    B_NMDA'     = -B_NMDA/tau_d_NMDA

    : VDCC (R-type)
    m_VDCC'     = (minf_VDCC-m_VDCC)/mtau_VDCC
    h_VDCC'     = (hinf_VDCC-h_VDCC)/htau_VDCC

    : Postsynaptic Ca2+ dynamics
    cai_CR'     = -(1e-9)*(ica_NMDA + ica_VDCC)*gamma_ca_CR/((1e-15)*volume_CR*2*FARADAY) - (cai_CR - min_ca_CR)/tau_ca_CR

    : Long-term synaptic plasticity
    Rho_GB'     = ( - Rho_GB*(1-Rho_GB)*(rho_star_GB-Rho_GB)
                    + potentiate_GB*gamma_p_GB*(1-Rho_GB)
                    - depress_GB*gamma_d_GB*Rho_GB ) / ((1e3)*tau_GB)
    Use_GB'     = (Use_d_GB + Rho_GB*(Use_p_GB-Use_d_GB) - Use_GB) / ((1e3)*tau_Use_GB)
}


NET_RECEIVE (weight) {
    LOCAL result, ves, occu, Use_actual

    if(flag == 1) {
        : Flag 1, Initialize watch calls
        WATCH (cai_CR > theta_d_GB) 2
        WATCH (cai_CR < theta_d_GB) 3
        WATCH (cai_CR > theta_p_GB) 4
        WATCH (cai_CR < theta_p_GB) 5
        if(verbose > 0){
            UNITSOFF
            printf("Flag 1, Initialize watch calls\n")
            UNITSON
        }
    } else if(flag == 2) {
        depress_GB = 1
    } else if(flag == 3) {
        depress_GB = 0
    } else if(flag == 4) {
        potentiate_GB = 1
    } else if(flag == 5) {
        potentiate_GB = 0
    } else {
        : Do not perform any calculations if the synapse (netcon) is deactivated.
        : This avoids drawing from the random stream
        if(  !(weight > 0) ) {
            VERBATIM
            return;
            ENDVERBATIM
        }

        : Set synapse weight
        w = weight

        if(verbose > 0){
            UNITSOFF
            printf("t = %g, incoming spike at synapse %g\n", t, synapseID)
            UNITSON
        }

        : Select Use_GB, if long-term plasticity is enabled
        if(enable_GB == 1) {
            Use_actual = Use_GB
        } else {
            Use_actual = Use
        }

        : Calculate u at event
        if(Fac > 0) {
            : Update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
            u = u*exp(-(t - tsyn)/Fac)
            u = u + Use_actual*(1-u)
        } else {
            u = Use_actual
        }

        : Recovery
        FROM counter = 0 TO (unoccupied - 1) {
            : Iterate over all unoccupied sites and compute how many recover
            Psurv = exp(-(t-tsyn)/Dep)
            result = urand()
            if(result>Psurv) {
                occupied = occupied + 1   : recover a previously unoccupied site
                if(verbose > 0) {
                    UNITSOFF
                    printf("\tRecovered 1 vesicle, P = %g R = %g\n", Psurv, result)
                    UNITSON
                }
            }
        }

        ves = 0              : Initialize the number of released vesicles to 0
        occu = occupied - 1  : Store the number of occupied sites in a local var

        FROM counter = 0 TO occu {
            : iterate over all occupied sites and compute how many release
            result = urand()
            if(result<u) {
                : release a single site!
                occupied = occupied - 1  : Decrease the number of occupied sites
                ves = ves + 1            : Increase number of released vesicles
            }
        }

        : Update number of unoccupied sites
        unoccupied = Nrrp - occupied

        : Update tsyn
        : tsyn knows about all spikes, not only those that released
        : i.e. each spike can increase the u, regardless of recovered state.
        :      and each spike trigger an evaluation of recovery
        tsyn = t

        : Update state variables
        A_AMPA = A_AMPA + ves/Nrrp*factor_AMPA
        B_AMPA = B_AMPA + ves/Nrrp*factor_AMPA
        A_NMDA = A_NMDA + ves/Nrrp*factor_NMDA
        B_NMDA = B_NMDA + ves/Nrrp*factor_NMDA

        if ( verbose > 0 ) {
            UNITSOFF
            printf("\tReleased %g vesicles out of %g\n", ves, Nrrp)
            UNITSON
        }
    }
}


FUNCTION nernst(ci(mM), co(mM), z) (mV) {
    nernst = (1000) * R * (celsius + 273.15) / (z*FARADAY) * log(co/ci)
    if(verbose > 1) {
        UNITSOFF
        printf("nernst:%f R:%f celsius:%f \n", nernst, R, celsius)
        UNITSON
    }
}


PROCEDURE setRNG() {
    VERBATIM
    // For compatibility, allow for either MCellRan4 or Random123
    // Distinguish by the arg types
    // Object => MCellRan4, seeds (double) => Random123
    usingR123 = 0;
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_rel);
        uint32_t a2 = 0;

        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
        if (ifarg(2)) {
            a2 = (uint32_t)*getarg(2);
        }
        if (ifarg(3)) {
            *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, (uint32_t)*getarg(3));
        } else {
            *pv = nrnran123_newstream((uint32_t)*getarg(1), a2);
        }
        usingR123 = 1;
    } else if( ifarg(1) ) {   // not a double, so assume hoc object type
        void** pv = (void**)(&_p_rng_rel);
        *pv = nrn_random_arg(1);
    } else {  // no arg, so clear pointer
        void** pv = (void**)(&_p_rng_rel);
        *pv = (void*)0;
    }
    ENDVERBATIM
}


FUNCTION urand() {
    VERBATIM
    double value;
    if ( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng_rel);
    } else if (_p_rng_rel) {
        #if !defined(CORENEURON_BUILD)
        value = nrn_random_pick(_p_rng_rel);
        #endif
    } else {
        value = 0.0;
    }
    _lurand = value;
    ENDVERBATIM
}


FUNCTION toggleLTPlasticity() {
    enable_GB = 1-enable_GB
}


FUNCTION toggleVerbose() {
    verbose = 1-verbose
}


FUNCTION bbsavestate() {
        bbsavestate = 0
VERBATIM
#ifdef ENABLE_SAVE_STATE
        /* first arg is direction (0 save, 1 restore), second is array*/
        /* if first arg is -1, fill xdir with the size of the array */
        double *xdir, *xval, *hoc_pgetarg();
        long nrn_get_random_sequence(void* r);
        void nrn_set_random_sequence(void* r, int val);
        xdir = hoc_pgetarg(1);
        xval = hoc_pgetarg(2);
        if (_p_rng_rel) {
            // tell how many items need saving
            if (*xdir == -1) {  // count items
                if( usingR123 ) {
                    *xdir = 1.0;
                } else {
                    *xdir = 2.0;
                }
                return 0.0;
            } else if(*xdir ==0 ) {  // save
                if( usingR123 ) {
                    uint32_t seq;
                    char which;
                    nrnran123_getseq( (nrnran123_State*)_p_rng_rel, &seq, &which );
                    xval[0] = (double) seq;
                    xval[1] = (double) which;
                } else {
                    xval[0] = (double)nrn_get_random_sequence(_p_rng_rel);
                }
            } else {  // restore
                if( usingR123 ) {
                    nrnran123_setseq( (nrnran123_State*)_p_rng_rel, (uint32_t)xval[0], (char)xval[1] );
                } else {
                    nrn_set_random_sequence(_p_rng_rel, (long)(xval[0]));
                }
            }
        }
#endif
ENDVERBATIM
}

VERBATIM

static void bbcore_write(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
    // make sure offset array non-null
    if (iArray) {

        // get handle to random123 instance
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_rel);

        // get location for storing ids
        uint32_t* ia = ((uint32_t*)iArray) + *ioffset;

        // retrieve/store identifier seeds
        nrnran123_getids(*pv, ia, ia+1);
    }

    // increment integer offset (2 identifier), no double data
    *ioffset += 2;
    *doffset += 0;

}

static void bbcore_read(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {

    // make sure it's not previously set
    assert(!_p_rng_rel);

    uint32_t* ia = ((uint32_t*)iArray) + *ioffset;

    // make sure non-zero identifier seeds
    if (ia[0] != 0 || ia[1] != 0) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_rel);

        // get new stream
        *pv = nrnran123_newstream(ia[0], ia[1]);
    }

    // increment intger offset (2 identifiers), no double data
    *ioffset += 2;
    *doffset += 0;
}

ENDVERBATIM
