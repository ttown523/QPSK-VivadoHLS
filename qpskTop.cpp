#include "qpskTop.h"

bool qpskElementDemodulatorTimingPhase(double sampleIn, Symbol *out){

	static DownsampleCounter downsampleCount = 0; //carrierIndex = 0,
	static TwoBitCounter carrierIndex = 0;

	double Irec = 0.0, Qrec = 0.0, Ifir = 0.0, Qfir = 0.0;
	double ICorrected = 0, QCorrected = 0;
	bool strobe = false;

	Irec = sampleIn*Icarrier[carrierIndex];
	Qrec = sampleIn*Qcarrier[carrierIndex];

	//update carrier index
	carrierIndex++; //This two bit counter will wrap and keep counting so no need for checks

	//compute the next filter sample...need two filters, because they each use static delay values within the function
	simple_fir_filterI(&Ifir, Irec);
	simple_fir_filterQ(&Qfir, Qrec);

	if (downsampleCount == 0){
		timingPhaseCorrection(Ifir, Qfir, &ICorrected, &QCorrected, &strobe);

		//decision block -- perhaps I should make this its own function...so the reader can know what its doing
		if (strobe){
			//std::cout << ICorrected <<  " " << QCorrected << std::endl;
			(*out)[1] = (ICorrected > 0.0 ? 1 : 0);
			(*out)[0] = (QCorrected > 0.0 ? 1 : 0);
		}
		downsampleCount++;

		return strobe;
	}
	else{
		downsampleCount++;
		return false;
	}
}
