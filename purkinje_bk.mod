TITLE BK-type K+ current channel in the Purkinje soma

COMMENT

This implements the Calcium activated K+ current channel (BK-type) in the soma of the Purkinje cell.

This model is based on the following publications: 

1. Khaliq et al (2003). The contribution of resurgent sodium current to high-frequency firing in Purkinje neurons: an experimental and modelling study. J. Neuroscience, 23, 4899-4912.

2. Forrest Michael D., Intracellular calcium dynamics permit a Purkinje neuron model to perform toggle and gain computations upon its inputs, Frontiers in Computational Neuroscience, Vol 8, 2014     

	
ENDCOMMENT


NEURON {
    SUFFIX pc_BK
    USEION k READ ek WRITE ik
    USEION ca READ cai
    RANGE gmax, ik
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
	tau_z   = 1 (ms)
}

PARAMETER {
	gmax     = 0.00728 (S/cm2)
}

ASSIGNED {
	v          (mV)
	ek         (mV)
    ik         (mA/cm2)
	cai        (mM)
	m_inf
	h_inf
	z_inf
	tau_m      (ms)
	tau_h      (ms)
}

STATE { 
	m FROM 0 TO 1
	h FROM 0 TO 1
	z FROM 0 TO 1
}

INITIAL {
	rates(v, cai)
	m = m_inf
	h = h_inf
	z = z_inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	ik = gmax * m^3 * z^2 * h * (v - ek)
}

DERIVATIVE states {
	rates(v, cai)
	m' = (m_inf - m)/tau_m
	h' = (h_inf - h)/tau_h
	z' = (z_inf - z)/tau_z
}

FUNCTION taum( v (mV) ) (ms) {
	UNITSOFF
	taum = 0.000505 + (1/(exp((v+86.4)/10.1) + exp((v-33.3)/(-10))))
	UNITSON
}

FUNCTION tauh( v (mV) ) (ms) {
	UNITSOFF
	tauh = 0.0019 + (1/(exp((v+48.5)/5.2) + exp((v-54.2)/(-12.9))))
	UNITSON
}

PROCEDURE rates( v_m (mV), cai (mM) ) {
    LOCAL v
	v = v_m + 5  : this accounts for junction potential

	m_inf = 1 / (1+exp(-((v-(-28.9))/6.2)))
	h_inf = 0.085 + ((1-0.085)/(1+exp(-((v-(-32))/(-5.8)))))
	z_inf = 1/(1 + (0.001/cai))

	tau_m = 1000 * taum(v)
	tau_h = 1000 * tauh(v)
}


