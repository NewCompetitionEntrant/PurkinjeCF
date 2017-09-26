TITLE Calcium concentration in the Purkinje cell soma

COMMENT

This implements the Ca2+ concentration in the soma of the Purkinje cell.

This model is based on the following publications: 

1. Khaliq et al (2003). The contribution of resurgent sodium current to high-frequency firing in Purkinje neurons: an experimental and modelling study. J. Neuroscience, 23, 4899-4912.

2. Forrest Michael D., Intracellular calcium dynamics permit a Purkinje neuron model to perform toggle and gain computations upon its inputs, Frontiers in Computational Neuroscience, Vol 8, 2014     

	
ENDCOMMENT


NEURON {
    SUFFIX pc_cai
    USEION ca READ ica WRITE cai
    RANGE  beta, c_depth, c_area, cai
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
}

PARAMETER {
	beta      = 1 (/ms)           : Ca2+ diffusion constant
	c_depth   = 0.1 (um)         : depth in the cell
	c_area    = 1521 (um2)       : membrane surface area
	cai0      = 100e-6 (mM)
}

ASSIGNED {
	v          (mV)
    ica        (mA/cm2)
}

STATE { cai }

INITIAL {
	cai = cai0
}

BREAKPOINT {
	SOLVE integrate METHOD derivimplicit

	if (cai < 1e-4) {
		cai = 1e-4
	}
}

DERIVATIVE integrate {
	cai' = ((((-100 * ica)/(2*F*c_depth*c_area)) * 1e6) - (beta*cai))
}






