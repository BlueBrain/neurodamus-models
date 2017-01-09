TITLE Leak potassium current
: From Amarillo et al., 2014

NEURON {
	SUFFIX TC_Kleak
	USEION  k READ ek WRITE ik
	RANGE g, i_rec
}

:CONSTANT {
	:Q10 = 3 (1) : To check, recordings at room temperature
:}

UNITS {
	(mA) = (milliamp)
	:(uA) = (microamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	:ek 			(mV) :EI: moved in assigned
	g = 1.0e-5	(S/cm2)	<0,1e9>
}

ASSIGNED {
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	:qt (1)
	i_rec
}


BREAKPOINT {
	ik = g*(v - ek)
	i_rec = ik
}







