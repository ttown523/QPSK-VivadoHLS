#include "qpsk.h"
#include <math.h>
#include "ap_shift_reg.h"
/*
			|
		01	|	11
	----------------
		00	|	10
			|
*/

/*
*	Matched filter element-by-element implementation
*/
//could implement the filter a different way since it is symmetric, but I
//for now I am going to leave it like this.
void simple_fir_filterI(double *y, double x){
	static ap_shift_reg<double,FILTER_TAPS> shiftRegister; //use predefined shift register class

	double acc = 0.0, data = 0.0;

	shiftRegister.shift(x);

	for (int i = FILTER_TAPS - 1; i >= 0; i--){
		acc += shiftRegister.read(i)*firCoeff[i];
	}

	*y = acc;
}

//
void simple_fir_filterQ(double *y, double x){
	static ap_shift_reg<double,FILTER_TAPS> shiftRegister; //use predefined shift register class

	double acc = 0.0, data = 0.0;

	shiftRegister.shift(x);

	for (int i = FILTER_TAPS - 1; i >= 0; i--){
		acc += shiftRegister.read(i)*firCoeff[i];
	}

	*y = acc;
}
/*
*	Need a second one for the static variable
*/
/*
void simple_fir_filterQ(double *y, double x){
	static ap_shift_reg<double,FILTER_TAPS> shiftRegister; //use predefined shift register class

	double acc = 0.0, data = 0.0;
	int j = FILTER_TAPS-1, i = 0;

	shiftRegister.shift(x);

	//96 is the filter coefficient that doesn't have a buddy
	//for (i = FILTER_TAPS - 1; i >= 0; i--){
	acc = shiftRegister.read(96);
	for (i = 0; i < 96; i++){
		acc += (shiftRegister.read(i)+shiftRegister.read(j))*firCoeff[i];
	}

	*y = acc;
}
*/


/************************************************
*												*
*		Timing and phase loop functions			*
*												*
*************************************************/

void timingPhaseCorrection(double MF_I, double MF_Q, double *ICorrected, double *QCorrected, bool *strobe){
	//registers used for this function
	static double ccwC_d = 0.0, ccwS_d = 0.0, ccwIOut_d = 0.0, ccwQOut_d = 0.0;
	static double interpolationIOut_d2 = 0.0, interpolationQOut_d2 = 0.0, timingControl_d = 0.0, mu_d = 0.0;
	static bool strobe_d = false;

	double ccwIOut = 0.0, ccwQOut = 0.0;
	double interpI1Out = 0.0, interpI2Out = 0.0, interpQ1Out = 0.0, interpQ2Out = 0.0;
	double tedOut = 0.0, pedOut = 0.0, timeLoopOut = 0.0, phaseLoopOut = 0.0;
	double interpolationControlOut = 0.0, ddsOutC = 0.0, ddsOutS = 0.0, mu = 0.0;
	bool underflow = false;

	*strobe = strobe_d;

	ccwRotation(MF_I, MF_Q, ccwC_d, ccwS_d, &ccwIOut, &ccwQOut);

	farrowInterpolationQuadraticI(ccwIOut, mu_d, &interpI1Out);
	farrowInterpolationQuadraticI2(ccwIOut_d, mu_d, &interpI2Out);
	farrowInterpolationQuadraticQ(ccwQOut, mu_d, &interpQ1Out);
	farrowInterpolationQuadraticQ2(ccwQOut_d, mu_d, &interpQ2Out);

	*ICorrected = interpI1Out;
	*QCorrected = interpQ1Out;

	TED(interpI1Out, interpI2Out, interpolationIOut_d2, interpQ1Out, interpQ2Out, interpolationQOut_d2, strobe_d, &tedOut);
	PED(interpI1Out, interpQ1Out, strobe_d, &pedOut);

	timingLoop(tedOut, &timeLoopOut);
	phaseLoop(pedOut, &phaseLoopOut);

	interpolationControl(timeLoopOut, &interpolationControlOut, &underflow);
	ddsFrequency(phaseLoopOut, &ddsOutC, &ddsOutS);

	mu = timingControl_d * 2;

	//Update all of the registers with their new values for next time the function is called
	if (underflow)
		mu_d = mu;

	ccwC_d = ddsOutC;
	ccwS_d = ddsOutS;

	ccwIOut_d = ccwIOut;
	ccwQOut_d = ccwQOut;

	if (strobe_d){
		interpolationIOut_d2 = interpI1Out;
		interpolationQOut_d2 = interpQ1Out;
	} //else leave the registers untouched

	timingControl_d = interpolationControlOut;
	strobe_d = underflow;
}

void ccwRotation(double I, double Q, double C, double S, double *Iprime, double *Qprime){
	*Iprime = I*C + Q*S;
	*Qprime = Q*C - I*S;
}

void farrowInterpolationQuadraticI(double in, double mu, double *out){
	static double inputShiftReg[2] = { 0, 0 }, scaledShiftReg[3] = { 0, 0, 0 };

	double nextScaled = 0.0, tmp1 = 0.0,  tmp2 = 0.0;

	nextScaled = in * -0.5;
	//consider breaking this down further in to each of the individual parts to be better in hardware
	tmp1 = inputShiftReg[0] + nextScaled - scaledShiftReg[0] + scaledShiftReg[1] + scaledShiftReg[2];
	tmp2 = scaledShiftReg[0] - nextScaled + scaledShiftReg[1] - scaledShiftReg[2];

	*out = inputShiftReg[1] + ((tmp2*mu + tmp1)*mu);

	//shift the inputs of the shift registers
	inputShiftReg[1] = inputShiftReg[0];
	inputShiftReg[0] = in;

	scaledShiftReg[2] = scaledShiftReg[1];
	scaledShiftReg[1] = scaledShiftReg[0];
	scaledShiftReg[0] = nextScaled;
}

void farrowInterpolationQuadraticI2(double in, double mu, double *out){
	static double inputShiftReg[2] = { 0, 0 }, scaledShiftReg[3] = { 0, 0, 0 };

	double nextScaled = 0.0, tmp1 = 0.0, tmp2 = 0.0;

	nextScaled = in * -0.5;
	//consider breaking this down further in to each of the individual parts to be better in hardware
	tmp1 = inputShiftReg[0] + nextScaled - scaledShiftReg[0] + scaledShiftReg[1] + scaledShiftReg[2];
	tmp2 = scaledShiftReg[0] - nextScaled + scaledShiftReg[1] - scaledShiftReg[2];

	*out = inputShiftReg[1] + ((tmp2*mu + tmp1)*mu);

	//shift the inputs of the shift registers
	inputShiftReg[1] = inputShiftReg[0];
	inputShiftReg[0] = in;

	scaledShiftReg[2] = scaledShiftReg[1];
	scaledShiftReg[1] = scaledShiftReg[0];
	scaledShiftReg[0] = nextScaled;
}

void farrowInterpolationQuadraticQ(double in, double mu, double *out){
	static double inputShiftReg[2] = { 0, 0 }, scaledShiftReg[3] = { 0, 0, 0 };

	double nextScaled = 0.0, tmp1 = 0.0,  tmp2 = 0.0;

	nextScaled = in * -0.5;
	//consider breaking this down further in to each of the individual parts to be better in hardware
	tmp1 = inputShiftReg[0] + nextScaled - scaledShiftReg[0] + scaledShiftReg[1] + scaledShiftReg[2];
	tmp2 = scaledShiftReg[0] - nextScaled + scaledShiftReg[1] - scaledShiftReg[2];

	*out = inputShiftReg[1] + ((tmp2*mu + tmp1)*mu);

	//shift the inputs of the shift registers
	inputShiftReg[1] = inputShiftReg[0];
	inputShiftReg[0] = in;

	scaledShiftReg[2] = scaledShiftReg[1];
	scaledShiftReg[1] = scaledShiftReg[0];
	scaledShiftReg[0] = nextScaled;

}

void farrowInterpolationQuadraticQ2(double in, double mu, double *out){
	static double inputShiftReg[2] = { 0, 0 }, scaledShiftReg[3] = { 0, 0, 0 };

	double nextScaled = 0.0, tmp1 = 0.0,  tmp2 = 0.0;

	nextScaled = in * -0.5;
	//consider breaking this down further in to each of the individual parts to be better in hardware
	tmp1 = inputShiftReg[0] + nextScaled - scaledShiftReg[0] + scaledShiftReg[1] + scaledShiftReg[2];
	tmp2 = scaledShiftReg[0] - nextScaled + scaledShiftReg[1] - scaledShiftReg[2];

	*out = inputShiftReg[1] + ((tmp2*mu + tmp1)*mu);

	//shift the inputs of the shift registers
	inputShiftReg[1] = inputShiftReg[0];
	inputShiftReg[0] = in;

	scaledShiftReg[2] = scaledShiftReg[1];
	scaledShiftReg[1] = scaledShiftReg[0];
	scaledShiftReg[0] = nextScaled;
}

void TED(double x, double xd1, double xd2, double y, double yd1, double yd2, bool strobe, double *e){
	Sign x_sign = 0, xd2_sign = 0, y_sign = 0, yd2_sign = 0;

	if (strobe){
		//possibly make a new function for this sign creation?
		x_sign = (x > 0.0) ? 1 : -1;
		xd2_sign = (xd2 > 0.0) ? 1 : -1;
		y_sign = (y > 0.0) ? 1 : -1;
		yd2_sign = (yd2 > 0.0) ? 1 : -1;

		*e = (xd1* (xd2_sign - x_sign)) + (yd1*(yd2_sign - y_sign));
	}
	else{
		*e = 0;
	}
}

void PED(double I, double Q, bool strobe, double *out){
	Sign signI = 0, signQ = 0;

	if (strobe){
		signI = (I > 0.0) ? 1 : -1;
		signQ = (Q > 0.0) ? 1 : -1;

		*out = (Q*signI) - (I*signQ);
	}
	else{
		*out = 0;
	}
}

void phaseLoop(double e, double *v){
	static double integratorFeedback = 0.0;

	double integratorResult = 0.0, proportionalResult = 0;

	proportionalResult = e*K1_PHASE;
	integratorResult = e*K2_PHASE;
	integratorFeedback += integratorResult;

	*v = proportionalResult + integratorFeedback;
}

void timingLoop(double e, double *v){
	static double integratorFeedback = 0.0;

	double integratorResult = 0.0, proportionalResult = 0;

	proportionalResult = e*K1_TIMING;
	integratorResult = e*K2_TIMING;
	integratorFeedback += integratorResult;

	*v = proportionalResult + integratorFeedback;
}

//decrementing mod-1 counter
void interpolationControl(double w, double *reg, bool *underflow){
	static double delay = 0.0;

	double remainder1 = 0;

	*underflow = (delay < 0.0);

	//finding the decimal part of the delay value
	remainder1 = (delay - (short)delay);

	if (remainder1 < 0.0)
		remainder1 += 1;

	*reg = remainder1;

	delay = remainder1 - (w + 0.5);

}

void ddsFrequency(double F, double *cosine, double *sine){
	static double delay = 0.0;

	double tmp = 0.0;

	//look into using cosine with Vivado HLS...it is probably supported but make sure to find this out...CORDIC algorithm for FPGA...
	*cosine = cos(delay);
	*sine = sin(delay);

	tmp = F + delay;

	if (tmp < TWO_PI){
		delay = tmp;
	}
	else {
		delay = tmp - TWO_PI;
	}
}
