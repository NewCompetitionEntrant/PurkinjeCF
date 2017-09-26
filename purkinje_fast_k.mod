TITLE Fast K+ current channel in the Purkinje soma

COMMENT

This implements the highly TEA sensitive (fast) K+ channel in the soma of the Purkinje cell.

This model is based on the following publications: 

1. Khaliq et al (2003). The contribution of resurgent sodium current to high-frequency firing in Purkinje neurons: an experimental and modelling study. J. Neuroscience, 23, 4899-4912.

2. Forrest Michael D., Intracellular calcium dynamics permit a Purkinje neuron model to perform toggle and gain computations upon its inputs, Frontiers in Computational Neuroscience, Vol 8, 2014     

	
ENDCOMMENT


NEURON {
    SUFFIX pc_fast_K
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
    ik        (mA/cm2)
	m_inf
	h_inf
	tau_m      (ms)
	tau_h      (ms)
}

STATE { m h }

INITIAL {
	rates(v)
	m = m_inf
	h = h_inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	ik = gmax * m^3 * h * (v - ek)
}

DERIVATIVE states {
	rates(v)
	m' = (m_inf - m)/tau_m
	h' = (h_inf - h)/tau_h
}

FUNCTION taum( v (mV) ) (ms) {
	UNITSOFF
	if (v < -35) {
		taum = 0.000103 + (0.0149*exp(0.035*v))
	} else {
		taum = 0.000129 + (1/(exp((v+100.7)/12.9) + exp((v-56)/(-23.1))))
	}
	UNITSON
}

FUNCTION tauh( v (mV) ) (ms) {
	UNITSOFF
	if (v <= 0) {
		tauh = 1.22e-5 + 0.012*exp(-(((v+56.3)/49.6)^2))
	} else {
		tauh = 0.0012 + 0.0023*exp(-0.141*v)
	}
	UNITSON
}

PROCEDURE rates( v_m (mV) ) {
	LOCAL v
	v = v_m + 11 : this accounts for junction potential

	m_inf = 1 / (1+exp(-((v-(-24))/15.4)))
	h_inf = 0.31 + ((1-0.31)/(1+exp(-((v-(-5.8))/(-11.2)))))

	tau_m = 1000 * taum(v)
	tau_h = 1000 * tauh(v)
}


