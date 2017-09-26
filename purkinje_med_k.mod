TITLE Medium K+ current channel in the Purkinje soma

COMMENT

This implements the moderately TEA sensitive (medium) K+ channel in the soma of the Purkinje cell.

This model is based on the following publications: 

1. Khaliq et al (2003). The contribution of resurgent sodium current to high-frequency firing in Purkinje neurons: an experimental and modelling study. J. Neuroscience, 23, 4899-4912.

2. Forrest Michael D., Intracellular calcium dynamics permit a Purkinje neuron model to perform toggle and gain computations upon its inputs, Frontiers in Computational Neuroscience, Vol 8, 2014     

	
ENDCOMMENT


NEURON {
    SUFFIX pc_medium_K
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
	gmax     = 0.00208 (S/cm2)
}

ASSIGNED {
	v          (mV)
	ek         (mV)
    ik        (mA/cm2)
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
	if (v < -20) {
		taum = 0.000688 + (1/(exp((v+64.2)/6.5) + exp((v-141.5)/(-34.8))))
	} else {
		taum = 0.00016 + 0.0008*exp(-0.0267*v)
	}
	UNITSON
}

PROCEDURE rates( v_m (mV) ) {
	LOCAL v
	v = v_m + 11  : this accounts for junction potential

	m_inf = 1 / (1+exp(-((v-(-24))/20.4)))

	tau_m = 1000 * taum(v)
}


