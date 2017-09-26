TITLE Slow K+ current channel in the Purkinje soma

COMMENT

This implements the TEA sensitive (slow) K+ channel in the soma of the Purkinje cell.

This model is based on the following publications: 

1. Khaliq et al (2003). The contribution of resurgent sodium current to high-frequency firing in Purkinje neurons: an experimental and modelling study. J. Neuroscience, 23, 4899-4912.

2. Forrest Michael D., Intracellular calcium dynamics permit a Purkinje neuron model to perform toggle and gain computations upon its inputs, Frontiers in Computational Neuroscience, Vol 8, 2014     

	
ENDCOMMENT


NEURON {
    SUFFIX pc_slow_K
    USEION k READ ek WRITE ik
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
}

PARAMETER {
	gmax     = 0.00416 (S/cm2)
}

ASSIGNED {
	v          (mV)
	ek         (mV)
    ik         (mA/cm2)
	m_inf
	tau_m      (ms)
}

STATE { m }

INITIAL {
	rates(v)
	m = m_inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	ik = gmax * m^4 * (v - ek)
}

DERIVATIVE states {
	rates(v)
	m' = (m_inf - m)/tau_m
}

FUNCTION taum( v (mV) ) (ms) {
	UNITSOFF
	taum = 0.000796 + (1/(exp((v+73.2)/11.7) + exp((v-306.7)/(-74.2))))
	UNITSON
}

PROCEDURE rates( v_m (mV) ) {
	LOCAL v
	v = v_m + 11  : this accounts for junction potential

	m_inf = 1 / (1+exp(-((v-(-16.5))/18.4)))

	tau_m = 1000 * taum(v)
}


