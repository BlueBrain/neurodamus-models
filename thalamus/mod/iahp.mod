TITLE Slow Ca-dependent potassium current
:
:   Ca++ dependent K+ current IC responsible for slow AHP
:   Differential equations
:
:   Model based on a first order kinetic scheme
:
:      <closed> + n cai <-> <open>	(alpha,beta)
:
:   Following this model, the activation fct will be half-activated at 
:   a concentration of Cai = (beta/alpha)^(1/n) = cac (parameter)
:
:   The mod file is here written for the case n=2 (2 binding sites)
:   ---------------------------------------------
:
:   This current models the "slow" IK[Ca] (IAHP): 
:      - potassium current
:      - activated by intracellular calcium
:      - NOT voltage dependent
:
:   A minimal value for the time constant has been added
:
:   Ref: Destexhe et al., J. Neurophysiology 72: 803-818, 1994.
:   See also: http://www.cnl.salk.edu/~alain , http://cns.fmed.ulaval.ca
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iahp
	USEION k READ ek WRITE ik
	USEION ca READ cai
        RANGE gkbar, m_inf, tau_m
	GLOBAL beta, cac
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}


PARAMETER {
	v		(mV)
	celsius	= 36	(degC)
	ek	= -95	(mV)
	cai 	= 2.4e-4 (mM)		: initial [Ca]i
	gkbar	= .01	(mho/cm2)
	beta	= 0.03	(1/ms)		: backward rate constant
	cac	= 0.025	(mM)		: middle point of activation fct
	taumin	= 0.1	(ms)		: minimal value of the time cst
}


STATE {
	m
}

ASSIGNED {
	ik	(mA/cm2)
	m_inf
	tau_m	(ms)
	tadj
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
    : added vshift to make slow AHP deep enough
	ik = gkbar * m*m * (v - ek + 15)
}

DERIVATIVE states { 
	evaluate_fct(v,cai)

	m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {
:
:  activation kinetics are assumed to be at 22 deg. C
:  Q10 is assumed to be 3
:
	tadj = 3 ^ ((celsius-22.0)/10)

	evaluate_fct(v,cai)
	m = m_inf
}

PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL car

	car = (cai/cac)^2

	m_inf = car / ( 1 + car )
	tau_m = 1 / beta / (1 + car) / tadj

        if(tau_m < taumin) { tau_m = taumin } 	: min value of time cst
}
UNITSON
