: Dynamics that track inside calcium concentration
: modified from Destexhe et al. 1994

NEURON	{
	SUFFIX CaDynamics_E3
	USEION ca READ ica WRITE cai
	RANGE decay, gamma, minCai, depth
}

UNITS	{
	(mV) = (millivolt)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulombs)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um)	= (micron)
}

PARAMETER	{

	depth = 0.1 (um) : depth of shell

	: Percent of free calcium (not buffered), Helmchen et al. 1996: (L5) kappa_e = 100
	gamma = 0.01 (1)

	: rate of removal of calcium:
	: L5 proximal dendrite: (Helmchen 1996) 70ms @ 7um diam
	decay = 70 (ms)

	: baseline calcium
	: CA1: (Sabatini 2002) 65 nM
	minCai = 6.5e-5 (mM)

}

ASSIGNED	{ica (mA/cm2)}

STATE	{
	cai (mM)
	}

BREAKPOINT	{ SOLVE states METHOD cnexp }

INITIAL {
	cai = minCai
}

DERIVATIVE states	{
	cai' = -(10000)*(ica*gamma/(2*FARADAY*depth)) - (cai - minCai)/decay
}
