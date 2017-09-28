:From Amarillo et al., 2014

NEURON	{
	SUFFIX TC_ih
	NONSPECIFIC_CURRENT ih
	RANGE gh_max, g_h, i_rec 
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gh_max = 2.2e-5(S/cm2) 
	e_h =  -43.0 (mV)
	celsius (degC)
}

ASSIGNED	{
	v	(mV)
	ih	(mA/cm2)
	g_h	(S/cm2)
	mInf
	mTau
	:tcorr		:Add temperature correction
	i_rec
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	g_h = gh_max*m
	ih = g_h*(v-e_h)
	i_rec = ih
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
	:tcorr = 4^((celsius-34)/10)  :EI: Recording temp. 34 C in Amarillo, 2014. q10 = 4 (Santoro et al., 2000).
                                      : Already corrected in equations 
}

PROCEDURE rates(){
	:UNITSOFF
	mInf = 1/(1+exp((v+82)/5.49))
	mTau = 1/((0.0008+0.0000035*exp(-0.05787*v)+exp(-1.87+0.0701*v)))
	:UNITSON
}
