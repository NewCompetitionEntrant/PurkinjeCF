TITLE P-type Calcium Channel in the Purkinje soma

COMMENT

This implements the P-type Ca2+ channel in the soma of the Purkinje cell.

This model is based on the following publications: 

1. Khaliq et al (2003). The contribution of resurgent sodium current to high-frequency firing in Purkinje neurons: an experimental and modelling study. J. Neuroscience, 23, 4899-4912.

2. Forrest Michael D., Intracellular calcium dynamics permit a Purkinje neuron model to perform toggle and gain computations upon its inputs, Frontiers in Computational Neuroscience, Vol 8, 2014     

	
ENDCOMMENT


NEURON {
    SUFFIX pc_ptype_ca
    USEION ca WRITE ica
    RANGE gmax, ica
}

UNITS {
    (mV)    = (millivolt)
    (mA)    = (milliamp)
    (nA)    = (nanoamp)
    (pA)    = (picoamp)
    (S)     = (siemens)
    (nS)    = (nanosiemens)
    (pS)    = (picosiemens)
    (um)    = (micron)
    (molar) = (1/liter)
    (mM)    = (millimolar)
}

CONSTANT {
	F       = 9.6485e4 (coulombs)
	R       = 8.3145 (joule/kelvin)
}

PARAMETER {
	gmax     = 0.00052 (S/cm2)
	init_cai = 0.0001 (mM)
	init_cao = 2 (mM)
	T        = 295 (kelvin)
	P_ca     = 5e-5 (cm/s)
}

ASSIGNED {
	v          (mV)
    ica        (mA/cm2)
	cai        (mM)
	m_inf
	tau_m      (ms)
	ghk_val    (coulombs/cm3)
	E          (volt)
}

STATE { m }

INITIAL {
	rates(v)
	m = m_inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	ica = (1e3) * gmax * m * ghk_val
}

DERIVATIVE states {
	rates(v)
	m' = (m_inf - m)/tau_m
}

FUNCTION ghk(v (mV), init_cai (mM), init_cao (mM), P_ca) (coulombs/cm3) {
	LOCAL z

	E = v * (1e-3)
	z = exp((-2*F*E)/(R*T))

	if (fabs(1 - z) > 1e-6) {
		ghk = (1e-3) * (4 * P_ca) * ((E*F^2)/(R*T))*((init_cai - (init_cao*z))/(1-z))
	} else {
		ghk = (1e-3) * (4 * P_ca) * ((E*F^2)/(R*T))*((init_cai - (init_cao*z))/(1e-6))
	}
}

FUNCTION taum( v (mV) ) (ms) {
	UNITSOFF
	if (v <= -50) {
		taum = 0.000264 + 0.128*exp(0.103*v)
	} else {
		taum = 0.000191 + 0.00376*exp(-(((v+41.9)/27.8)^2))
	}
	UNITSON
}

PROCEDURE rates( v (mV) ) {

	m_inf = 1 / (1+exp(-((v-(-19.0))/5.5)))

	tau_m = 1000 * taum(v)

	ghk_val = ghk(v, init_cai, init_cao, P_ca)
}


