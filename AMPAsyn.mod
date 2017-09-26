TITLE AMPA synapse connecting the Climbing Fiber to the Purkinje Cell soma

COMMENT

This implements the AMPA synapse which connects the Climbing Fiber to the Purkinje Cell soma

This file also implements the mathematical model of the autophosphorylation of CaMKII in the soma of the Purkinje cell.

This model is based on the following publications: 

1. Zhabotinsky, A., Bistability in the Ca2+/Calmodulin dependent protein kinase-phosphatase system, BiosPhysical Journal, Vol 79, Novemeber 2000, pp 2211-2221

2. Wu, A.,  Lawrence, D. A., Monotonoicty and bistability of Calcium/Calmodulin-dependent protein kinase-phosphatase activation, 2010 American Control Conference, Baltimore, MD, USA, June 30-July 02, 2010

ENDCOMMENT


NEURON {
    POINT_PROCESS AMPAsyn
    RANGE g, ina, CaMKII, synapse_cai
    USEION na WRITE ina
    USEION ca READ cai
}

UNITS {
    (mV)    = (millivolt)
    (mA)    = (milliamp)
    (nA)    = (nanoamp)
    (pA)    = (picoamp)
    (S)     = (siemens)
    (mS)    = (millisiemens)
    (nS)    = (nanosiemens)
    (pS)    = (picosiemens)
    (um)    = (micron)
    (molar) = (1/liter)
    (mM)    = (millimolar)
}

CONSTANT {
}


PARAMETER {
	tau = 0.1       (ms)
	ginit_AMPAR = 0.00085 (S/cm2)
	e_AMPAR = 0 (mV)

	: Duration of the Ca2+ internal concentration change
	cai_duration = 14 (ms)

	: Steady state Ca2+
	cai_inf  = 100e-6 (mM)

	e_k = 20e-3 (mM)
	e_p_0 = 0.1e-3 (mM)
	I_0 = 0.1e-3 (mM)
	v_CaN = 1e-3 (/ms)
	v_PKA = 1e-3 (/ms)
	K_M = 0.4e-3 (mM)
	K_H1 = 4.0e-3 (mM)
	K_H2 = 0.7e-3 (mM)
	k1 = 0.5e-3 (/ms)
	k2 = 2.0e-3 (/ms)
	k3 = 1.0 (/ms)
	k4 = 0.001e-3 (/ms)
}

ASSIGNED {
	v             (millivolt)
    ina           (mA/cm2)
	cai           (mM)
    AMPAR_count
	CaMKII (mM)
    t_stimulus (ms)
	synapse_cai (mM)
	cai_inrush

	w_a
	w_d
	w_1
	w_2
	w_3
	w_4
	w_5
	w_6
	w_7
	w_8
	w_9
	alpha_0
	alpha_1
	alpha_2
	alpha_3
	alpha_4
	alpha_5
	alpha_6
	alpha_7
	alpha_8
	alpha_9
	delta_x
	I_tmp
}


STATE { 
	g (microsiemens) 
	P_0    (mM)
	P_1    (mM)
	P_2    (mM)
	P_3    (mM)
	P_4    (mM)
	P_5    (mM)
	P_6    (mM)
	P_7    (mM)
	P_8    (mM)
	P_9    (mM)
	P_f10   (mM)
	e_p   (mM)
	I     (mM)
}

INITIAL { 
	g = 0 
	P_0 = e_k
	P_1 = 0
	P_2 = 0
	P_3 = 0
	P_4 = 0
	P_5 = 0
	P_6 = 0
	P_7 = 0
	P_8 = 0
	P_9 = 0
	P_f10 = 0
	e_p = 0
	I = 0
	
    AMPAR_count = 60
	CaMKII = P_0
	t_stimulus = 0
	cai_inrush = 0
	synapse_cai=cai_inf
}


BREAKPOINT {

	: Only the Ca2+ concentration at this synapse will be used for
	: the CaMKII production
	if (cai_inrush == 1) {
		if ((t - t_stimulus) < cai_duration) {
			synapse_cai = cai
		} else {
			cai_inrush = 0
			synapse_cai = cai_inf
		}
	} else {
		synapse_cai = cai_inf
	}

	SOLVE states METHOD cnexp

	ina = g*AMPAR_count*(v - e_AMPAR)
	CaMKII = 1*P_1 + 2*P_2 + 3*P_3 + 4*P_4 + 5*P_5 + 6*P_6 + 7*P_7 + 8*P_8 + 9*P_9 + 10*P_f10

    : printf("P_0 = %g, P_1 = %g, P_2 = %g, P_3 = %g, P_4 = %g, CaMKII = %g, cai = %g\n", P_0, P_1, P_2, P_3, P_4, CaMKII, cai)
    : printf("P_5 = %g, P_6 = %g, P_7 = %g, P_8 = %g, P_9 = %g, P_f10 = %g, e_p = %g, I = %g\n\n", P_5, P_6, P_7, P_8, P_9, P_f10, e_p, I)
}

DERIVATIVE states {

	g' = -g/tau
	: printf("time: %g, g: %g\n", t, g)

	rates(synapse_cai)

	sum_Pi()
	P_1' = (-alpha_0*P_1) - (delta_x*w_d*(P_1-P_2)) + (alpha_0*P_0)

	sum_Pi()
	P_2' = (alpha_1*w_a*(P_1-P_2)) - (2*delta_x*w_d*(P_2-P_3))

	sum_Pi()
	P_3' = (alpha_2*w_a*(P_2-P_3)) - (3*delta_x*w_d*(P_3-P_4))

	sum_Pi()
	P_4' = (alpha_3*w_a*(P_3-P_4)) - (4*delta_x*w_d*(P_4-P_5))

	sum_Pi()
	P_5' = (alpha_4*w_a*(P_4-P_5)) - (5*delta_x*w_d*(P_5-P_6))

	sum_Pi()
	P_6' = (alpha_5*w_a*(P_5-P_6)) - (6*delta_x*w_d*(P_6-P_7))

	sum_Pi()
	P_7' = (alpha_6*w_a*(P_6-P_7)) - (7*delta_x*w_d*(P_7-P_8))

	sum_Pi()
	P_8' = (alpha_7*w_a*(P_7-P_8)) - (8*delta_x*w_d*(P_8-P_9))

	sum_Pi()
	P_9' = (alpha_8*w_a*(P_8-P_9)) - (9*delta_x*w_d*(P_9-P_f10))

	sum_Pi()
	P_f10' = (alpha_9*w_a*(P_9-P_f10)) - (10*delta_x*w_d*P_f10)

	e_p' = -k3*I*e_p + k4*(e_p_0 - e_p)
	I' = -k3*I*e_p + k4*(e_p_0 - e_p) + v_PKA*I_0 - I_tmp
}

PROCEDURE rates( cai_val (mM) ) {
	w_a = (((cai_val/K_H1)^4))/(1+((cai_val/K_H1)^4))
	w_d = e_p
	w_1 =  1.0
	w_2 = 1.8
	w_3 = 2.3
	w_4 = 2.7
	w_5 = 2.8
	w_6 = 2.7
	w_7 = 2.3
	w_8 = 1.8
	w_9 = 1.0
	alpha_0 = 10*k1*(w_a^2)
	alpha_1 = w_1 * k1
	alpha_2 = w_2 * k1
	alpha_3 = w_3 * k1
	alpha_4 = w_4 * k1
	alpha_5 = w_5 * k1
	alpha_6 = w_6 * k1
	alpha_7 = w_7 * k1
	alpha_8 = w_8 * k1
	alpha_9 = w_9 * k1
	I_tmp = v_CaN*I*((((cai_val/K_H2)^3))/(1+((cai_val/K_H2)^3)))
}

PROCEDURE sum_Pi() {
	LOCAL sum_Pi_val

	sum_Pi_val = P_1 + P_2 + P_3 + P_4 + P_5 + P_6 + P_7 + P_8 + P_9 + P_f10
	delta_x = k2/(K_M+sum_Pi_val)
}

NET_RECEIVE (weight (microsiemens)) {
	g = ginit_AMPAR
    AMPAR_count = AMPAR_count + (weight*10)

	: printf("time: %g, weight: %g, g: %g, AMPAR_count: %g\n", t, weight, g, AMPAR_count)

	t_stimulus = t   : save the neurotransmitter release start time
	cai_inrush = 1
}


