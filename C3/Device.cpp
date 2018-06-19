#include "Device.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

Device::Device()
{
}

void Device::pathRestoration(string inputFile, string outputFile) {
	
	pack_input lastReadPack;
	pack_output lastOutputPack;

	int count = 0;

	ifstream file(inputFile);
	ofstream file_out(outputFile);

	file >> count;

	readPackOfFile(file, lastReadPack);
	lastReadPack.accel_z -= g;

	writePackToFile(file_out, lastOutputPack);

	for(int i = 1; i < count; i++) {
		pack_input readBuffer;

		readPackOfFile(file, readBuffer);
		readBuffer.accel_z -= g;

		lastOutputPack.v_x = readBuffer.v_x;

		///////////////
		algorithmPathRestoration(readBuffer, lastReadPack, lastOutputPack);
		///////////////

		writePackToFile(file_out, lastOutputPack);

		lastReadPack = readBuffer;
	}
	file.close();
	file_out.close();
}

void Device::algorithmPathRestoration(pack_input &inputNext, pack_input &inputFirst, pack_output &output) {
	//output.v_x

	double delta_time = inputNext.t - inputFirst.t;

	output.v_y += delta_time * inputFirst.accel_y;
	output.v_z += delta_time * inputFirst.accel_z;

	backRotate(inputNext, inputFirst, output);
}

void Device::writePackToFile(ofstream &out, pack_output &packOut) {
	out << packOut.t
		<< " " << packOut.x
		<< " " << packOut.y
		<< " " << packOut.z
		<< std::endl;
}

void Device::readPackOfFile(ifstream &inp, pack_input &packInp) {
	inp >> packInp.t;
	inp >> packInp.fi_x;
	inp >> packInp.fi_y;
	inp >> packInp.fi_z;
	inp >> packInp.accel_x;
	inp >> packInp.accel_y;
	inp >> packInp.accel_z;
	inp >> packInp.v_x;
}

void Device::backRotate(pack_input &inputNext, pack_input &inputFirst, pack_output &output) {
	double delta_time = inputNext.t - inputFirst.t;

	double alfa = inputFirst.fi_x,
		beta = inputFirst.fi_y,
		gamma = inputFirst.fi_z;

	//в неподвижной системе координат
	double v_x = output.v_x * (cos(alfa)*cos(gamma) - sin(alfa)*cos(beta)*sin(gamma))
		+ output.v_y * (sin(alfa) * cos(gamma) + cos(alfa) * cos(beta) * sin(gamma))
		+ output.v_z * (sin(beta) * sin(gamma));

	double v_y = output.v_x * (-cos(alfa)*sin(gamma) - sin(alfa) * cos(beta)*cos(gamma))
		+ output.v_y * (-sin(alfa) * sin(gamma) + cos(alfa)*cos(beta)*cos(gamma))
		+ output.v_z * (sin(beta) * cos(gamma));

	double v_z = output.v_x * (sin(alfa)*sin(beta))
		+ output.v_y * (-cos(alfa)*sin(beta))
		+ output.v_z * cos(beta);

	output.t = inputNext.t;
	output.x += delta_time * v_x;
	output.y += delta_time * v_y;
	output.z += delta_time * v_z;
}

Device::~Device()
{
}
