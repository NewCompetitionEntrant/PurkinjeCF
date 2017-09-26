TITLE Resurgent Na+ Channel in the Purkinje soma

COMMENT

This implements the Resurgent Na+ channel in the soma of the Purkinje cell.

This model is based on the following publications: 

1. Khaliq et al (2003). The contribution of resurgent sodium current to high-frequency firing in Purkinje neurons: an experimental and modelling study. J. Neuroscience, 23, 4899-4912.

2. Forrest Michael D., Intracellular calcium dynamics permit a Purkinje neuron model to perform toggle and gain computations upon its inputs, Frontiers in Computational Neuroscience, Vol 8, 2014     

3. Raman, I.M., and Bean,B.P. (2001). Inactivation and recovery of sodium currents in cerebellar Purkinje neurons:evidence for twomechanisms. Biophys.j. 80, 729–737


	
ENDCOMMENT


NEURON {
    SUFFIX pc_res_Na
    USEION na READ ena WRITE ina
    RANGE gmax, ina
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

PARAMETER {
	gmax     = 0.0156 (S/cm2)
}

CONSTANT {
	Con        = 0.005
	Coff       = 0.5
	Oon        = 0.75
	: Oon        = 2.3
	Ooff       = 0.005

	: a = (Oon/Con)^(1/4)
	rt_a       = 3.4996
	: rt_a       = 4.6312

	: b = (Ooff/Coff)^(1/4)
	rt_b       = 0.3162
}

ASSIGNED {
	v          (mV)
	ena        (mV)
    ina        (mA/cm2)

	cf1         (/ms)
	cf2         (/ms)
	cf3         (/ms)
	cf4         (/ms)
	cf5         (/ms)
	cb1         (/ms)
	cb2         (/ms)
	cb3         (/ms)
	cb4         (/ms)
	cb5         (/ms)
	if1         (/ms)
	if2         (/ms)
	if3         (/ms)
	if4         (/ms)
	if5         (/ms)
	ib1         (/ms)
	ib2         (/ms)
	ib3         (/ms)
	ib4         (/ms)
	ib5         (/ms)
	icf1         (/ms)
	icf2         (/ms)
	icf3         (/ms)
	icf4         (/ms)
	icf5         (/ms)
	icb1         (/ms)
	icb2         (/ms)
	icb3         (/ms)
	icb4         (/ms)
	icb5         (/ms)
	of1          (/ms)
	ob1          (/ms)
	iof6         (/ms)
	iob6         (/ms)
}

STATE { 
	c1  FROM 0 TO 1 
	c2  FROM 0 TO 1
	c3  FROM 0 TO 1
	c4  FROM 0 TO 1
	c5  FROM 0 TO 1
	i1  FROM 0 TO 1
	i2  FROM 0 TO 1
	i3  FROM 0 TO 1
	i4  FROM 0 TO 1
	i5  FROM 0 TO 1
	i6  FROM 0 TO 1
	o   FROM 0 TO 1
	ob  FROM 0 TO 1
}

INITIAL { 
	rates(v)

	SOLVE init_kin
}

KINETIC kin {
	rates(v)
	~ c1 <-> c2 (cf1, cb1)
	~ c2 <-> c3 (cf2, cb2)
	~ c3 <-> c4 (cf3, cb3)
	~ c4 <-> c5 (cf4, cb4)
	~ c5 <-> o (cf5, cb5)

	~ o <-> ob (of1, ob1)

	~ i1 <-> i2 (if1, ib1)
	~ i2 <-> i3 (if2, ib2)
	~ i3 <-> i4 (if3, ib3)
	~ i4 <-> i5 (if4, ib4)
	~ i5 <-> i6 (if5, ib5)

	~ c1 <-> i1 (icf1, icb1)
	~ c2 <-> i2 (icf2, icb2)
	~ c3 <-> i3 (icf3, icb3)
	~ c4 <-> i4 (icf4, icb4)
	~ c5 <-> i5 (icf5, icb5)

	~ o <-> i6 (iof6, iob6)

	CONSERVE c1+c2+c3+c4+c5+i1+i2+i3+i4+i5+i6+ob+o=1
}

LINEAR init_kin {
 ~          i1*icb1 + c2*cb1  - c1*(    icf1+cf1) = 0
 ~ c1*cf1 + i2*icb2 + c3*cb2  - c2*(cb1+icf2+cf2) = 0
 ~ c2*cf2 + i3*icb3 + c4*cb3  - c3*(cb2+icf3+cf3) = 0
 ~ c3*cf3 + i4*icb4 + c5*cb4  - c4*(cb3+icf4+cf4) = 0
 ~ c4*cf4 + i5*icb5 + o*cb5   - c5*(cb4+icf5+cf5) = 0
 ~ c5*cf5 + ob*ob1  + i6*iob6 - o*(cb5+of1+iof6) = 0
 ~ o*of1  - ob*ob1 = 0

 ~          c1*icf1 + i2*ib1  - i1*(    icb1+if1) = 0
 ~ i1*if1 + c2*icf2 + i3*ib2  - i2*(ib1+icb2+if2) = 0
 ~ i2*if2 + c3*icf3 + i4*icb3 - i3*(ib2+icb3+if3) = 0
 ~ i3*if3 + c4*icf4 + i5*ib4  - i4*(ib3+icb4+if4) = 0
 ~ i4*if4 + c5*icf5 + i6*ib5  - i5*(ib4+icb5+if5) = 0
 
 ~ c1 + c2 + c3 + c4 + c5 + o + ob + i1 + i2 + i3 + i4 + i5 + i6 = 1
}

BREAKPOINT {
	SOLVE kin METHOD sparse

	ina = gmax * o * (v-ena)
}

PROCEDURE rates( v (mV) ) {
	LOCAL rt_alpha, rt_beta, rt_zeta, rt_gamma, rt_delta, rt_epsilon

	rt_alpha = 150 * exp(v/20)
	rt_beta = 3 * exp(v/(-20))
	rt_zeta = 0.03 * exp(v/(-25))
	rt_gamma   = 150 * exp(v/1e12)
	rt_delta   = 40 * exp(v/(-1e12))
	rt_epsilon = 1.75 * exp(v/1e12)
	: rt_epsilon = 1e-12 * exp(v/1e12)

	cf1 = 4 * rt_alpha
	cb1 = rt_beta

	cf2 = 3 * rt_alpha
	cb2 = 2 * rt_beta

	cf3 = 2 * rt_alpha
	cb3 = 3 * rt_beta

	cf4 = rt_alpha
	cb4 = 4 * rt_beta

	cf5 = rt_gamma
	cb5 = rt_delta

	of1 = rt_epsilon
	ob1 = rt_zeta

	if1 = 4 * rt_alpha * rt_a
	ib1 = rt_beta * rt_b
	
	if2 = 3 * rt_alpha * rt_a
	ib2 = 2 * rt_beta * rt_b

	if3 = 2 * rt_alpha * rt_a
	ib3 = 3 * rt_beta * rt_b

	if4 = rt_alpha * rt_a
	ib4 = 4 * rt_beta * rt_b

	if5 = rt_gamma
	ib5 = rt_delta

	icf1 = Con
	icb1 = Coff

	icf2 = Con * rt_a
	icb2 = Coff * rt_b

	icf3 = Con * (rt_a^2)
	icb3 = Coff * (rt_b^2)

	icf4 = Con * (rt_a^3)
	icb4 = Coff * (rt_b^3)

	icf5 = Con * (rt_a^4)
	icb5 = Coff * (rt_b^4)

	iof6 = Oon
	iob6 = Ooff
}









