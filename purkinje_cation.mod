TITLE Hyperpolarization activated cation current channel in the Purkinje soma

COMMENT

This implements the hyperpolarization activated cation channel in the soma of the Purkinje cell.

This model is based on the following publications: 

1. Khaliq et al (2003). The contribution of resurgent sodium current to high-frequency firing in Purkinje neurons: an experimental and modelling study. J. Neuroscience, 23, 4899-4912.

2. Forrest Michael D., Intracellular calcium dynamics permit a Purkinje neuron model to perform toggle and gain computations upon its inputs, Frontiers in Computational Neuroscience, Vol 8, 2014     

	
ENDCOMMENT


NEURON {
    SUFFIX pc_cation
	NONSPECIFIC_CURRENT i
    RANGE gmax, eh
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
	gmax     = 0.000104 (S/cm2)
	eh       = -30 (mV)
}

ASSIGNED {
	v          (mV)
    i          (mA/cm2)
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

	i = gmax * m * (v - eh)
}

DERIVATIVE states {
	rates(v)
	m' = (m_inf - m)/tau_m
}

FUNCTION taum( v (mV) ) (ms) {
	UNITSOFF
	taum = 0.19 + 0.72*exp(-(((v+81.5)/11.9)^2))
	UNITSON
}

PROCEDURE rates( v (mV) ) {

	m_inf = 1 / (1+exp(-((v-(-90.1))/(-9.9))))

	tau_m = 1000 * taum(v)
}


