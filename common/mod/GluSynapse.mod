COMMENT
@file GluSynapse.mod
@brief Two state deterministic model of a Glutamatergic Synapse
@author chindemi, king, roessert
@date 07/07/2016
@version 0.3.4dev
@acknowledgement Michael Hines, for reviewing many important aspects of model
                 implementation, optimization, and all the technical support.
                 Francesco Cremonesi, for helping with DE simplification.
                 Christian Roessert, for implementing the VDCC and calcium
                 dynamics model, with spine volume normalization.
@remark Copyright © BBP/EPFL 2005-2016; All rights reserved.
        Do not distribute without further notice.
ENDCOMMENT


TITLE Glutamatergic Synapse


COMMENT
In brief:
    Two state deterministic model of a Glutamatergic Synapse.
             ----------------------------------
             |                                |
             |                                v
    Functional Synapse             Potential Synapse
             ^                                |
             |                                |
             ----------------------------------

Functional Synapse components:
    ## AMPAR ##
        Adapted from "ProbAMPANMDA_EMS.mod"
        Suffix: _AMPA

    ## NMDAR ##
        Adapted from "ProbAMPANMDA_EMS.mod"
        Fractional calcium current from Schneggenburger et al. 1993
        Magnesium dynamics from Jahr and Stevens 1990
        Suffix: _NMDA

    ## VDCC ##
        Christian Roessert, based on data from various sources, described inline
        Suffix: _VDCC

    ## Postsynaptic Ca2+ dynamics ##
        Christian Roessert, based on data from various sources, described inline
        Suffix: _CR

    ## Spine ##
        Christian Roessert, based on data from various sources, described inline
        Suffix: _CR

    ## Release ##
        Adapted from "ProbAMPANMDA_EMS.mod"
        See Fuhrmann et al. 2002

    ## Multi-vesicular release ##
        John Rahmon, unpublished
        Suffix: _MVR

    ## LTP/LTD ##
        See Graupner and Brunel 2012
        Suffix: _GB

    ## Elimination ##
        Test
        Suffix: _SE

Potential Synapse components:
    ## Synaptogenesis  ##
        Test
        Suffix: _SG

ENDCOMMENT


NEURON {
    THREADSAFE : TODO Check this keyword

    POINT_PROCESS GluSynapse

    : AMPAR range variables
    RANGE tau_r_AMPA, tau_d_AMPA, gmax_AMPA, factor_AMPA
    :RANGE i_AMPA, g_AMPA, q_AMPA

    : NMDAR range variables
    RANGE tau_r_NMDA, tau_d_NMDA, factor_NMDA
    :RANGE i_NMDA, g_NMDA, q_NMDA
    RANGE mg, mgnl_NMDA
    :RANGE mggate_NMDA
    RANGE Pf_NMDA, ica_NMDA

    : VDCC
    RANGE gca_bar_VDCC, gca_bar_abs_VDCC, ica_VDCC, ljp_VDCC
    :RANGE gca_VDCC
    RANGE vhm_VDCC, km_VDCC, mtau_VDCC, minf_VDCC
    RANGE vhh_VDCC, kh_VDCC, htau_VDCC, hinf_VDCC
    RANGE tm_VDCC, th_VDCC

    : Postsynaptic Ca2+ dynamics range variables
    RANGE cai_CR, cao_CR, gamma_ca_CR, tau_ca_CR, min_ca_CR

    : Spine range variable
    RANGE volume_CR, area_CR

    : Release range variables
    RANGE Use, u, Dep, Fac, u0, tsyn
    :RANGE Psurv
    POINTER rng_rel

    : Multi-vesicular release (MVR)
    RANGE unoccupied_MVR, occupied_MVR, N_MVR, q_MVR

    : LTP/LTD range variables
    RANGE tau_GB, rho_star_GB, rho0_GB
    RANGE gamma_p_GB, gamma_d_GB, theta_p_GB, theta_d_GB
    RANGE w0_GB, w1_GB
    RANGE Theta_d_GB, Theta_p_GB
    RANGE tau_Use_GB
    RANGE gamma_eff_GB

    : Elimination range variables
    :RANGE tau_SE, theta_e_SE

    : Synaptogenesis range variables
    :RANGE tau_SG, theta_g_SG

    : Shared range variables
    RANGE e, v, NMDA_ratio, vv, g
    :RANGE eca_syn

    : Other range variables
    RANGE synapseID, verboseLevel, LTPlasticity, synapseState
    :RANGE rewiring

    : Misc
    NONSPECIFIC_CURRENT i
}


UNITS {
    (nA)    = (nanoamp)
    (mV)    = (millivolt)
    (molar) = (1/liter)
    (um)    = (micrometers)
    (mM)    = (millimolar)
    (uS)    = (microsiemens)
    (nS)    = (nanosiemens)
    (S)     = (siemens)
    (pC)    = (picocoulomb)
    FARADAY = (faraday) (coulomb)
    PI      = (pi) (1)
    R       = (k-mole) (joule/degC)
}


PARAMETER {
    celsius                (degC): TODO

    : AMPAR parameters
    tau_r_AMPA   = 0.2     (ms)  : Dual-exponential conductance profile
    tau_d_AMPA   = 1.7     (ms)  : IMPORTANT: tau_r < tau_d
    gmax_AMPA    = 1.5     (nS)  : Overwritten by BlueBuilder assigned values
    :q_AMPA       = 1       (nS)  : Quantal size

    : NMDAR parameters
    tau_r_NMDA   = 0.29    (ms)  : Dual-exponential conductance profile
    tau_d_NMDA   = 43      (ms)  : IMPORTANT: tau_r < tau_d
    mg           = 1       (mM)  : Initial concentration of mg2+
    mgnl_NMDA    = 0.062   (/mV) : Nonlinearity factor
    Pf_NMDA      = 0.068   (1)   : Fractional NMDAR calcium current

    : VDCC
    ljp_VDCC     = 0       (mV)  :

    gca_bar_VDCC = 0.165   (nS/um2) :

    vhm_VDCC     = -5.9    (mV)  : v 1/2 for act, Magee and Johnston 1995 (corrected for m*m)
    km_VDCC      = 9.5     (mV)  : act slope, Magee and Johnston 1995 (corrected for m*m)
    vhh_VDCC     = -39     (mV)  : v 1/2 for inact, Magee and Johnston 1995
    kh_VDCC      = -9.2    (mV)  : inact, Magee and Johnston 1995

    tm_VDCC      = 1       (ms)  : max. time constant (guess)
    th_VDCC      = 27      (ms)  : max. time constant 100*0.27

    : Postsynaptic Ca2+ dynamics parameters (FIXED from Pub)
    gamma_ca_CR  = 0.04    (1)   : Percent of free calcium (not buffered), Sabatini et al. 2002: kappa_e = 24+-11 (also 14 (2-31) or 22 (18-33))
    tau_ca_CR    = 12      (ms)  : Rate of removal of calcium, Sabatini et al. 2002: 14ms (12-20ms)
    min_ca_CR    = 70e-6   (mM)  : Sabatini et al. 2002: 70+-29 nM, per AP: 1.1 (0.6-8.2) uM = 1100 e-6 mM = 1100 nM
    volume_CR    = 0.092   (um3) : From spine data by Ruth Benavides-Piccione (unpublished), default value overwritten at runtime time
    cao_CR       = 2.0     (mM)

    : Release parameters, just initial values! Use,
    : Dep and Fac are overwritten by BlueBuilder assigned values
    Use          = 1.0     (1)  : Utilization of synaptic efficacy
    Dep          = 100     (ms) : Relaxation time constant from depression
    Fac          = 10      (ms) : Relaxation time constant from facilitation
    u0           = 0       (1)  : Initial value of u, which is the running value
                                : of release probability

    : MVR parameters, just initial values!
    N_MVR        = 1       (1)  : Number of total release sites for given contact

    : LTP/LTD parameters
    w0_GB        = 0.2     (1)
    w1_GB        = 0.8     (1)
    rho_star_GB  = 0.5     (1)
    rho0_GB      = -1.0    (1)
    tau_GB       = 100.0   (s)
    theta_d_GB   = 0.008   (mM)
    theta_p_GB   = 0.012   (mM)
    gamma_d_GB   = 100.0   (1)
    gamma_p_GB   = 450.0   (1)
    tau_Use_GB   = 100.0   (s)

    : Elimination
    :tau_SE       = 1.0     (s)
    :theta_e_SE   = 0.2     (1)
    :factor_SE    = 0.5     (1)

    : Synaptogenesis
    :tau_SG       = 1.0     (s)
    :theta_g_SG   = 0.4     (1)
    :step_SG      = 0.01    (1)

    : Shared parameters
    e            = 0       (mV) : AMPA and NMDA reversal potential
    :eca_syn      = 136     (mV) : The value is updated in the BREAKPOINT
    NMDA_ratio   = 0.71    (1)  : The ratio of NMDA to AMPA

    : Misc
    synapseID    = 0       (1)
    verboseLevel = 0       (1)
    LTPlasticity = 0       (1)
    :rewiring     = 0       (1)
    synapseState = 2       (1)  : [0: Potential Synapse, 2: Functional Synapse]
}


COMMENT
The Verbatim block is needed to allow RNG.
ENDCOMMENT
VERBATIM
// for MCellRan4
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

// for random123
#include "nrnran123.h"

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM


ASSIGNED {
    : AMPAR assigned variables
    :g_AMPA      (uS)
    factor_AMPA (1)
    :i_AMPA      (nA)

    : NMDAR assigned variables
    :g_NMDA      (uS)
    :q_NMDA      (nS) : Computed in NET_RECEIVE, relative to gmax_AMPA. In the future we might want to add a delay to account for fast AMPA-mediated LTP/LTD
    factor_NMDA (1)
    :i_NMDA      (nA)
    ica_NMDA    (nA)
    :mggate_NMDA (1)

    : VDCC assigned variables
    :gca_VDCC    (uS)
    gca_bar_abs_VDCC (nS)
    ica_VDCC    (nA)
    minf_VDCC   (1)
    hinf_VDCC   (1)
    mtau_VDCC   (ms)
    htau_VDCC   (ms)

    : Spine assigned variables
    area_CR     (um2)

    : Release assigned variables
    tsyn        (ms) : the time of the last spike
    :Psurv       (1)
    u           (1)  : running release probability (attention: u is event based based, so only valid at incoming events)

    : MVR
    unoccupied_MVR (1) : no. of unoccupied sites following release event
    occupied_MVR   (1) : no. of occupied sites following one epoch of recovery

    : LTP/LTD assigned variables
    Theta_d_GB  (1)
    Theta_p_GB  (1)
    gamma_eff_GB (1)

    : Elimination assigned variables
    : None

    : Synaptogenesis assigned variables
    : None

    : Shared assigned variables
    dt          (ms)
    v           (mV)
    vv          (mV)
    i           (nA)
    g           (uS)

    : Misc
    rng_rel

    : temporary until mcellran4 completely deprecated
    usingR123
}

STATE {
    : AMPAR state variables to construct the dual-exponential profile
    A_AMPA       (nS) : Decays with conductance tau_r_AMPA
    B_AMPA       (nS) : Decays with conductance tau_d_AMPA

    : NMDAR state variables to construct the dual-exponential profile
    A_NMDA       (nS) : Decays with conductance tau_r_NMDA
    B_NMDA       (nS) : Decays with conductance tau_d_NMDA

    : VDCC
    m_VDCC       (1)
    h_VDCC       (1)

    : Postsynaptic Ca2+ dynamics state variables
    cai_CR       (mM) : Intracellular calcium concentration

    : LTP/LTD state variables
    rho_GB       (1)

    : Elimination state variables
    :integrity_SE (1)

    : Synaptogenesis state variables
    :contact_SG
}


INITIAL {
    LOCAL tp_AMPA, tp_NMDA

    if( synapseState == 0 ) {
        initialize_potential_synapse()
    } else if( synapseState == 2 ) {
        initialize_functional_synapse()
        if(rho0_GB == -1.0) {
            : Initial LTP
            rho_GB = (Use - w0_GB) / (w1_GB - w0_GB)
        } else {
            rho_GB = rho0_GB
        }
    }

    : Time to peak of the conductances
    tp_AMPA = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA)
    tp_NMDA = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA)

    : AMPA Normalization factor - so that when t = tp_AMPA, gsyn = gpeak
    factor_AMPA = -exp(-tp_AMPA/tau_r_AMPA)+exp(-tp_AMPA/tau_d_AMPA)
    factor_AMPA = 1/factor_AMPA

    : NMDA Normalization factor - so that when t = tp_NMDA, gsyn = gpeak
    factor_NMDA = -exp(-tp_NMDA/tau_r_NMDA)+exp(-tp_NMDA/tau_d_NMDA)
    factor_NMDA = 1/factor_NMDA

    gamma_eff_GB = 0

    : VDCC
    UNITSOFF
    area_CR = 4 * PI * (3/4 * volume_CR * 1/PI)^(2/3)  : convert volume to area, assuming sphere for spine head
    UNITSON
    gca_bar_abs_VDCC = gca_bar_VDCC * area_CR
    mtau_VDCC = tm_VDCC
    htau_VDCC = th_VDCC

    net_send(0, 1)
}


BEFORE BREAKPOINT {
    LOCAL step_gamma
    COMMENT
    This implementation of the Graupner model is slightly different from the
    original. The choice of precomputing Theta_d_GB/Theta_p_GB might affect the
    convergence of the model.
    ENDCOMMENT

    step_gamma = exp(dt*(( - 1.0 ) / 150.0))
    gamma_eff_GB = gamma_eff_GB*step_gamma + cai_CR * (1 - step_gamma)

    Theta_d_GB = Theta(cai_CR - theta_d_GB)*1000.0*gamma_eff_GB
    Theta_p_GB = Theta(cai_CR - theta_p_GB)*1000.0*gamma_eff_GB
}


BREAKPOINT {
    LOCAL g_AMPA, i_AMPA, mggate_NMDA, g_NMDA, i_NMDA, eca_syn, gca_VDCC

    SOLVE state METHOD cnexp

    : AMPAR
    g_AMPA = (1e-3)*(B_AMPA-A_AMPA)                               : compute time varying conductance as the difference of state variables B_AMPA and A_AMPA
    i_AMPA = g_AMPA*(v-e)                                         : compute the AMPA driving force based on the time varying conductance, membrane potential, and AMPA reversal

    :NMDAR
    mggate_NMDA = 1 / (1 + exp(mgnl_NMDA * -(v)) * (mg / 3.57 (mM))) : mggate kinetics - Jahr & Stevens 1990
    g_NMDA = (1e-3)*(B_NMDA-A_NMDA) * mggate_NMDA                 : compute time varying conductance as the difference of state variables B_NMDA and A_NMDA and mggate kinetics
    i_NMDA = g_NMDA*(v-e)                                         : compute the NMDA driving force based on the time varying conductance, membrane potential, and NMDA reversal
    ica_NMDA = Pf_NMDA * g_NMDA * (v-40.0)                        : compute calcium component of NMDAR-mediated current, note that calcium
                                                                  : reversal potential is independent of extracellular calcium concentration
                                                                  : (Jahr and Stevens 1993)

    g = g_AMPA + g_NMDA

    : VDCC
    eca_syn = nernst(cai_CR, cao_CR, 2)                           : Ca reversal potential
    gca_VDCC = (1e-3) * gca_bar_abs_VDCC * m_VDCC * m_VDCC * h_VDCC
    ica_VDCC = gca_VDCC*(v-eca_syn)
    minf_VDCC = 1 / (1 + exp(((vhm_VDCC - ljp_VDCC) - v) / km_VDCC))
    hinf_VDCC = 1 / (1 + exp(((vhh_VDCC - ljp_VDCC) - v) / kh_VDCC))

    : Total current
    i = i_AMPA + i_NMDA + ica_VDCC                                : TODO check if ica_VDCC corrupts the fitting

    : Spine voltage (unused)
    vv = v
}


AFTER SOLVE {
    LOCAL target_Use, step_Use
    if( LTPlasticity == 1 ) {
        target_Use = w0_GB + rho_GB*(w1_GB - w0_GB)
        step_Use = exp(dt*(( - 1.0 ) / ((1e3)*tau_Use_GB)))
        Use = Use*step_Use + target_Use * (1 - step_Use)
    }
}


DERIVATIVE state{
    LOCAL rho0, a, b

    : AMPAR
    A_AMPA' = -A_AMPA/tau_r_AMPA
    B_AMPA' = -B_AMPA/tau_d_AMPA

    : NMDAR
    A_NMDA' = -A_NMDA/tau_r_NMDA
    B_NMDA' = -B_NMDA/tau_d_NMDA

    : VDCC
    m_VDCC' = (minf_VDCC-m_VDCC)/mtau_VDCC
    h_VDCC' = (hinf_VDCC-h_VDCC)/htau_VDCC

    : Postsynaptic Ca2+ dynamics
    cai_CR' = -(1e-9) * (ica_NMDA + ica_VDCC) * gamma_ca_CR/((1e-15)*volume_CR*2*FARADAY) - (cai_CR - min_ca_CR)/tau_ca_CR

    : LTP/LTD
    rho0 = rho_GB
    a = ( -rho0*(1.0 - rho0)*(rho_star_GB - rho0)
          +gamma_p_GB*(1-rho0)*Theta_p_GB
          -gamma_d_GB*rho0*Theta_d_GB ) / ((1e3)*tau_GB)
    b = ( -3.0*rho0*rho0 + 2*(1 + rho_star_GB)*rho0 -rho_star_GB
          -gamma_p_GB*Theta_p_GB
          -gamma_d_GB*Theta_d_GB ) / ((1e3)*tau_GB)
    rho_GB' = a + b*(rho_GB - rho0)

    : Elimination
    :integrity_SE' = -integrity_SE/((1e3)*tau_SE*(1+rho0))

    : Synaptogenesis
    :contact_SG' = -contact_SG/((1e3)*tau_SG)
}


NET_RECEIVE (weight){
    LOCAL result, ves, occu, q_AMPA, q_NMDA, Psurv

    INITIAL{
        : TODO Check if this block is executed at the time of the first "true"
        :      spike or before the watch calls initialization. In the latter
        :      case, this is probably wrong.
    }

    if (flag == 1) {
        : Flag 1, Initialize watch calls
        :WATCH (integrity_SE < theta_e_SE) 4
        :WATCH (contact_SG > theta_g_SG) 5
        :printf("Flag 1, Initialize watch calls\n")

    :} else if(flag == 4 && rewiring == 1 && synapseState == 2){
    :    : Flag 4, eliminate synapse
    :    synapseState = 0
    :    initialize_potential_synapse()
    :    :printf("Flag 4, eliminate synapse\n")

    :} else if(flag == 5 && rewiring == 1){
    :    : Flag 5, new functional synapse
    :    synapseState = 2
    :    initialize_functional_synapse()
    :    :printf("Flag 5, new functional synapse\n")

    } else if (flag == 0 && synapseState == 2){
        : Flag 0 (default) and synapse in the functional state
        q_AMPA = gmax_AMPA / N_MVR
        q_NMDA = q_AMPA * NMDA_ratio : Here for backward compatibility

        : calc u at event
        if (Fac > 0) {
            : Update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
            u = u*exp(-(t - tsyn)/Fac)
            u = u + Use*(1-u)
        } else {
            u = Use
        }

        : recovery
        FROM counter = 0 TO (unoccupied_MVR - 1) {
            : Iterate over all unoccupied sites and compute how many recover
            Psurv = exp(-(t-tsyn)/Dep)
            result = urand()
            if (result>Psurv) {
                occupied_MVR = occupied_MVR + 1     : recover a previously unoccupied site
                if( verboseLevel > 0 ) {
                    UNITSOFF
                    printf( "Recovered! %f at time %g: Psurv = %g, urand=%g\n", synapseID, t, Psurv, result )
                    UNITSON
                }
            } else {
                if( verboseLevel > 0 ) {
                    UNITSOFF
                    printf( "Failed to recover! %f at time %g: Psurv = %g, urand=%g\n", synapseID, t, Psurv, result )
                    UNITSON
                }
            }
        }

        ves = 0                  : Initialize the number of released vesicles to 0
        occu = occupied_MVR - 1  : Store the number of occupied sites in a local variable

        FROM counter = 0 TO occu {
            : iterate over all occupied sites and compute how many release
            result = urand()
            if (result<u) {
                : release a single site!
                occupied_MVR = occupied_MVR - 1  : decrease the number of occupied sites by 1
                ves = ves + 1                    : increase number of relesed vesicles by 1

                if ( verboseLevel > 0 ) {
                    UNITSOFF
                    printf( "Release! %f at time %g: vals %g %g %g %g urand=%g\n", synapseID, t, A_AMPA, gmax_AMPA, factor_AMPA, weight, result )
                    UNITSON
                }
            } else {
                if( verboseLevel > 0 ) {
                    UNITSOFF
                    printf("Failure! %f at time %g: urand = %g\n", synapseID, t, result )
                    UNITSON
                }
            }
        }

        : Update number of unoccupied sites
        unoccupied_MVR = N_MVR - occupied_MVR

        : Update tsyn
        : tsyn knows about all spikes, not only those that released
        : i.e. each spike can increase the u, regardless of recovered state.
        :      and each spike trigger an evaluation of recovery
        tsyn = t

        if (ves > 0) { :no need to evaluate unless we have vesicle release
            A_AMPA = A_AMPA + ves*q_AMPA*factor_AMPA
            B_AMPA = B_AMPA + ves*q_AMPA*factor_AMPA
            A_NMDA = A_NMDA + ves*q_NMDA*factor_NMDA
            B_NMDA = B_NMDA + ves*q_NMDA*factor_NMDA

            : Update integrity variable
            :integrity_SE = integrity_SE + (1-integrity_SE)*factor_SE
        }

    } else if (flag == 0 && synapseState == 0){
        : Flag 0 (default) and synapse in the potential state
        :contact_SG = contact_SG + step_SG

    } else {
        : DO NOTHING
        : TODO Raise Exception
    }
}


PROCEDURE initialize_functional_synapse() {
    tsyn = 0
    u = u0

    A_AMPA = 0
    B_AMPA = 0

    A_NMDA = 0
    B_NMDA = 0

    : Postsynaptic Ca2+ concentration
    cai_CR = min_ca_CR

    : MVR
    unoccupied_MVR = 0
    occupied_MVR = N_MVR

    : LTP/LTD
    rho_GB = 0

    : Elimination
    :integrity_SE = 1

    :contact_SG = 0
}


PROCEDURE initialize_potential_synapse() {
    tsyn = 0
    u = u0

    A_AMPA = 0
    B_AMPA = 0

    A_NMDA = 0
    B_NMDA = 0

    : MVR
    unoccupied_MVR = N_MVR
    occupied_MVR = 0

    : Postsynaptic Ca2+ concentration
    cai_CR = min_ca_CR

    : LTP/LTD
    rho_GB = 0

    : Elimination
    :integrity_SE = 0

    :contact_SG = 0
}


FUNCTION Theta(x (mM)) (1) {
    if (x < 0.0) {
        Theta = 0.0
    } else {
        Theta = 1.0
    }
}


PROCEDURE setRNG() {
    VERBATIM
    // For compatibility, allow for either MCellRan4 or Random123.  Distinguish by the arg types
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


FUNCTION toggleVerbose() {
    verboseLevel = 1-verboseLevel
}


FUNCTION toggleLTPlasticity() {
    LTPlasticity = 1-LTPlasticity
}


:FUNCTION toggleRewiring() {
:    rewiring = 1-rewiring
:}


FUNCTION nernst(ci(mM), co(mM), z) (mV) {
    nernst = (1000) * R * (celsius + 273.15) / (z*FARADAY) * log(co/ci)

    if( verboseLevel > 0 ) {
        UNITSOFF
        :printf("nernst:%f R:%f celsius:%f \n", nernst, R, celsius)
        UNITSON
    }
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
