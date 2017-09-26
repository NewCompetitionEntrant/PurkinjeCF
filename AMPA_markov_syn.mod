TITLE AMPA synapse connecting the Climbing Fiber to the Purkinje Cell soma

COMMENT

This implements the AMPA synapse which connects the Climbing Fiber to the Purkinje Cell soma

This model is based on the following publications: 

P. Jonas, G. Major and B. Sakmann, Quantal components of unitary EPSCs at the mossy fibre synapse on CA3 Pyramidal cells
of the rat hippocampus, Journal of Physiology (1993), pp 615-663

ENDCOMMENT

NEURON {
    POINT_PROCESS AMPA_markov_syn
    USEION ca READ cai
    NONSPECIFIC_CURRENT i
    RANGE g, i, O, CaMKII, synapse_active, synapse_cai, N_AMPAR, CaMKII_thresh
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

PARAMETER {
    N_init_AMPAR = 325
	e_AMPAR = 0 (mV)

	: Maximum concentration of glutamate released by the climbing fiber
	: into the synapse
	max_glutamate = 3 (mM)

	: Initial single channel conductance of the AMPA receiver
	ginit_AMPAR = 8.5e-12 (S)

	: Maximum single channel conductance due to phophorylation
	gmax_AMPAR = 50e-12 (S)

	: Increase in AMPAR conductance due to phosphorylation
    : by CaMKII
	g_delta = 4e-12 (S)

	: Time constant for dephosphorylation of AMPAR
	tau_dephos=150 (ms)

	: Minimum CamKII concentration required for phosphorylation
	: of AMPARs
	CaMKII_thresh = 2.5e-13 (mM)

	: Duration of the Glutamate neurotransmitter release
    : into the synapse
    trans_duration = 0.1 (ms)

	: Duration of the Ca2+ internal concentration change
	cai_duration = 14 (ms)

	: Steady state Ca2+
	cai_inf  = 100e-6 (mM)

	: Parameters for the AMPAR-glutamate binding kinetics
	k_m1 = 4.26        (/ms)
	k_m2 = 3.26        (/ms)
	k_m3 = 0.0457      (/ms)
	alpha = 4.24       (/ms)
	beta = 0.90        (/ms)
	alpha_1 = 2.89     (/ms)
	beta_1 = 0.0392    (/ms)
	alpha_2 = 0.172    (/ms)
	beta_2 = 0.000727  (/ms)
	alpha_3 = 0.0177   (/ms)
	beta_3 = 0.004     (/ms)
	alpha_4 = 0.0168   (/ms)
	beta_4 = 0.1904    (/ms)

	: Parameters for the active CaMKII kinetics
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
	v          (millivolt)
    i          (nA)
	cai        (mM)
	glutamate  (mM)
    t_stimulus (ms)
	t_AMPAR_update (ms)
	synapse_cai (mM)
	synapse_active
	ca_inrush

	kp1_c      (/ms)
	kp2_c      (/ms)
	kp3_c      (/ms)

	CaMKII (mM)
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
	kinase_alpha_0
	kinase_alpha_1
	kinase_alpha_2
	kinase_alpha_3
	kinase_alpha_4
	kinase_alpha_5
	kinase_alpha_6
	kinase_alpha_7
	kinase_alpha_8
	kinase_alpha_9
	delta_x
	I_tmp

	g          (mS)
	g_stimulus (mS)
	N_AMPAR
	AMPAR_cntr
}


STATE { 
	C0 FROM 0 TO 1
	C1 FROM 0 TO 1
	C2 FROM 0 TO 1
	C3 FROM 0 TO 1
	C4 FROM 0 TO 1
	C5 FROM 0 TO 1
	O  FROM 0 TO 1

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
	g = ginit_AMPAR
    g_stimulus = 0
	t_stimulus = 0
	t_AMPAR_update = 0
	N_AMPAR=N_init_AMPAR
	AMPAR_cntr = 0
	glutamate = 0
	synapse_active = 0
	ca_inrush = 0
	synapse_cai=cai_inf
	rates()

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
	CaMKII = 0

	SOLVE init_ampa_kin
}

KINETIC ampa_kin {
	rates()

	~ C0 <-> C1 (kp1_c, k_m1)
	~ C1 <-> C2 (kp2_c, k_m2)
	~ C2 <-> O (alpha, beta)
	~ C1 <-> C3 (alpha_1, beta_1)
	~ C2 <-> C4 (alpha_2, beta_2)
	~ O <-> C5 (alpha_3, beta_3)
	~ C3 <-> C4 (kp3_c, k_m3)
	~ C4 <-> C5 (alpha_4, beta_4)

	CONSERVE C0+C1+C2+O+C3+C4+C5=1
}

LINEAR init_ampa_kin {
~                                                 C1*k_m1 - C0*kp1_c = 0
~       C0*kp1_c + C3*beta_1 + C2*k_m2 - C1*(k_m1 + kp2_c + alpha_1) = 0
~       C1*kp2_c + C4*beta_2 + O*beta - C2*(k_m2 + alpha_2 + alpha)  = 0
~                          C2*alpha + C5*beta_3 - O*(beta + alpha_3) = 0
~                         C1*alpha_1 + C4*k_m3 - C3*(beta_1 + kp3_c) = 0
~   C3*kp3_c + C2*alpha_2 + C5*beta_4 - C4*(k_m3 + beta_2 + alpha_4) = 0
: ~                      C4*alpha_4 + O*alpha_3 - C5*(beta_4 + beta_3) = 0

~ C0 + C1 + C2 + O + C3 + C4 + C5 = 1
}

BREAKPOINT {
	SOLVE ampa_kin METHOD sparse

	SOLVE kinase_states METHOD cnexp

	: Only the Ca2+ concentration at this synapse will be used for
	: the CaMKII phosphorylation
	if (ca_inrush == 1) {
		if ((t - t_stimulus) < cai_duration) {
			synapse_cai = cai
		} else {
			ca_inrush = 0
			synapse_cai = cai_inf
		}
	} else {
		synapse_cai = cai_inf
	}

	: Dephosphorylation of the AMPAR
	if (g > ginit_AMPAR) {
		g=g_stimulus*exp((t_stimulus-t)/tau_dephos)
	}

	if (g < ginit_AMPAR) {
		g=ginit_AMPAR
		g_stimulus=0
	}

	: AMPAR is an ionotropic receptor and it is permeable to Na+
    : (g is in S, v is in mV and i is in nA. Therefore, we have 
    : to multiply v by 10^6 to obtain the current in nA)
	i = 1e6*g*O*N_AMPAR*(v - e_AMPAR)

	: Active CaMKII is the sum of all N-fold phosphorylated holoenzymes
	CaMKII = 1*P_1 + 2*P_2 + 3*P_3 + 4*P_4 + 5*P_5 + 6*P_6 + 7*P_7 + 8*P_8 + 9*P_9 + 10*P_f10

	if ((t - t_AMPAR_update) > 1000) {
		update_AMPAR_count(CaMKII)
		t_AMPAR_update = t
	}
}

PROCEDURE update_AMPAR_count(CaMKII_val (mM)) {
	: These rates were determined by the MATLAB simulation
	if (CaMKII_val >= 0.02) {
		AMPAR_cntr = AMPAR_cntr + 0.5730
	} else if (CaMKII_val >= 0.01) {
		AMPAR_cntr = AMPAR_cntr + 0.2350
	} else if (CaMKII_val >= 0.005) {
		AMPAR_cntr = AMPAR_cntr + 0.0956
	} else if (CaMKII_val >= 0.0025) {
		AMPAR_cntr = AMPAR_cntr + 0.00242
	}

	if (AMPAR_cntr > 1.0) {
		N_AMPAR = N_AMPAR + 1
		AMPAR_cntr = 0
	}
}

PROCEDURE rates( ) {
	LOCAL delta_time

	delta_time = t - t_stimulus

	if (delta_time > trans_duration) {
		: If the time since the stimulus started is greater than the
        : synapse input pulse duration, this implies that no 
        : more glutamte is being released in the synapse
		glutamate = 0
		synapse_active = 0
	}

	: Set the AMPAR glutamate binding rates
	kp1_c = 4.59 * glutamate
	kp2_c = 28.4 * glutamate
	kp3_c = 1.27 * glutamate
}

DERIVATIVE kinase_states {

	kinase_rates(synapse_cai)

	sum_Pi()
	P_1' = (-kinase_alpha_0*P_1) - (delta_x*w_d*(P_1-P_2)) + (kinase_alpha_0*P_0)

	sum_Pi()
	P_2' = (kinase_alpha_1*w_a*(P_1-P_2)) - (2*delta_x*w_d*(P_2-P_3))
	
	sum_Pi()
	P_3' = (kinase_alpha_2*w_a*(P_2-P_3)) - (3*delta_x*w_d*(P_3-P_4))
	
	sum_Pi()
	P_4' = (kinase_alpha_3*w_a*(P_3-P_4)) - (4*delta_x*w_d*(P_4-P_5))
	
	sum_Pi()
	P_5' = (kinase_alpha_4*w_a*(P_4-P_5)) - (5*delta_x*w_d*(P_5-P_6))

	sum_Pi()
	P_6' = (kinase_alpha_5*w_a*(P_5-P_6)) - (6*delta_x*w_d*(P_6-P_7))
	
	sum_Pi()
	P_7' = (kinase_alpha_6*w_a*(P_6-P_7)) - (7*delta_x*w_d*(P_7-P_8))
	
	sum_Pi()
	P_8' = (kinase_alpha_7*w_a*(P_7-P_8)) - (8*delta_x*w_d*(P_8-P_9))
	
	sum_Pi()
	P_9' = (kinase_alpha_8*w_a*(P_8-P_9)) - (9*delta_x*w_d*(P_9-P_f10))
	
	sum_Pi()
	P_f10' = (kinase_alpha_9*w_a*(P_9-P_f10)) - (10*delta_x*w_d*P_f10)

	e_p' = -k3*I*e_p + k4*(e_p_0 - e_p)
	I' = -k3*I*e_p + k4*(e_p_0 - e_p) + v_PKA*I_0 - I_tmp
}

PROCEDURE kinase_rates( cai_val (mM) ) {
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
	kinase_alpha_0 = 10*k1*(w_a^2)
	kinase_alpha_1 = w_1 * k1
	kinase_alpha_2 = w_2 * k1
	kinase_alpha_3 = w_3 * k1
	kinase_alpha_4 = w_4 * k1
	kinase_alpha_5 = w_5 * k1
	kinase_alpha_6 = w_6 * k1
	kinase_alpha_7 = w_7 * k1
	kinase_alpha_8 = w_8 * k1
	kinase_alpha_9 = w_9 * k1
	I_tmp = v_CaN*I*((((cai_val/K_H2)^3))/(1+((cai_val/K_H2)^3)))
}

PROCEDURE sum_Pi() {
	LOCAL sum_Pi_val

	sum_Pi_val = P_1 + P_2 + P_3 + P_4 + P_5 + P_6 + P_7 + P_8 + P_9 + P_f10
	delta_x = k2/(K_M+sum_Pi_val)
}

NET_RECEIVE (weight (microsiemens)) {
	: printf("time: %g, weight: %g, g: %g\n", t, weight, g)

	: The amount of Glutamate neurotransmitter released depends 
	: on the strength of this synapse
	glutamate = (weight * max_glutamate)

	t_stimulus = t   : save the neurotransmitter release start time

	: Phosphorylation of AMPARs will occur if there is just enough
	: CaMKII
	if (CaMKII > CaMKII_thresh) {
		g=g+g_delta
		if (g > gmax_AMPAR) {
			g = gmax_AMPAR
		}
		g_stimulus = g
	}

	synapse_active = 1
	ca_inrush = 1
}


