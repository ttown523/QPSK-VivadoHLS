#include <iostream>
#include <stdio.h>
#include "qpskTop.h"
#include <fstream>
#include <cstdlib>

/************************************************************************************************
 *																								*
 * Testbench used for testing the correctness of the qpsk design in both hardware and software  *
 *																								*
 ************************************************************************************************/

int main () {

/*
	double data = 0, out1 = 0, out2 = 0;

	std::ifstream input("modulatedData.dat", std::ios::in);

	for(int i = 0; i < 1000; i++){
		input >> data;
		 simple_fir_filterI(&out1, data);
		 simple_fir_filterQ(&out2, data);

		 std::cout << out1 << ", " << out2 << std::endl;
	}

	return 0;
	*/

	//Open the file with the simulated modulation data
	std::ifstream input("modulatedData.dat", std::ios::in);
	//Create the result file
	std::ofstream outputFile("demodulatedMessage.dat", std::ios::trunc); // out.golden.dat

	if (!input.is_open()){ //make sure file is valid
		std::cout << "Error Opening File!" << std::endl;
	}

	double nextSample = 0.0;
	//short recievedSymbol[2] = { 0, 0 }, recv = 0;
	int fillFilter = 0;
	Symbol recievedSymbol;

	//Pass each sample from the input file into the demodulator
	std::cout << "Demodulating input signal..." << std::endl;
	while (input >> nextSample){
		if (qpskElementDemodulatorTimingPhase(nextSample, &recievedSymbol)) { //returns true if the next sample is valid
			//recv = (recievedSymbol[0] << 1) | recievedSymbol[1];
			if (fillFilter > 11){ //ignore the first few samples that we obtain, since we are waiting for the filter taps to fill up
				outputFile << recievedSymbol << std::endl;
			}
			else{
				fillFilter += 1;
			}
		}
	}

	//For some reason, we don't get the last symbol of the message in the previous loop.
	//So, we include this code so the last symbol is included and the file differencing will work
	nextSample = 0;
	while (true){
		if (qpskElementDemodulatorTimingPhase(nextSample, &recievedSymbol)) {
			//Convert the symbol to a number from 0-3
			//recv = (recievedSymbol[0] << 1) | recievedSymbol[1];
			outputFile << recievedSymbol << std::endl;

			break;
		}
	}

	input.close();
	outputFile.close();

	//Compare the generated message file with the golden file...FC is only valid in windows
	if (system("diff -w demodulatedMessage.dat golden.dat")){

		fprintf(stdout, "*******************************************\n");
		fprintf(stdout, "FAIL: Output DOES NOT match the golden output\n");
		fprintf(stdout, "*******************************************\n");

		return 1; //test failed
	}
	else {
		fprintf(stdout, "*******************************************\n");
		fprintf(stdout, "PASS: The output matches the golden output!\n");
		fprintf(stdout, "*******************************************\n");

		return 0; //test passed

	}
}
