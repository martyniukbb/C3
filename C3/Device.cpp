#include "Device.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

Device::Device()
{
}

void Device::pathRestoration(string inputFile, string outputFile) {

	pack_input state_input;
	pack_output state;

	int count = 0;

	ifstream file(inputFile);
	ofstream file_out(outputFile);

	file >> count;

	readPackOfFile(file, state_input);
	writePackToFile(file_out, state);

	for(int i = 1; i < count;) {

		vector<pack_input> lastReadPack;
		vector<pack_output> lastOutputPack;

		for (int j = 0; (j < numberPack) && (i < count); j++, i++) {
			pack_input readBuffer;
			readPackOfFile(file, readBuffer);
			lastReadPack.push_back(readBuffer);
		}

		///////////////
		algorithmPathRestoration(lastReadPack, lastOutputPack, state, state_input);
		///////////////

		for (int j = 0; j < int(lastOutputPack.size()); j++) {
			writePackToFile(file_out, lastOutputPack[j]);
		}
	}
	file.close();
	file_out.close();
}

void Device::algorithmPathRestoration(vector<pack_input> &input, vector<pack_output> &output, pack_output &state, pack_input &state_input) {
	
	for (int i = 0; i < input.size(); i++) {
		input[i].accel_z -= g;
	}

	input = smoothing(input);

	for (int i = 0; i < input.size(); i++) {

		pack_output outputBuffer;

		double delta_time = input[i].t - state_input.t;

		state.v_x = input[i].v_x;

		state.v_y += delta_time * input[i].accel_y;
		state.v_z += delta_time * input[i].accel_z;

		backRotate(input[i], state_input, state);

		outputBuffer.t = state.t;
		outputBuffer.v_x = state.v_x;
		outputBuffer.v_y = state.v_y;
		outputBuffer.v_z = state.v_z;
		outputBuffer.x = state.x;
		outputBuffer.y = state.y;
		outputBuffer.z = state.z;

		output.push_back(outputBuffer);
		state_input = input[i];
	}
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

vector<pack_input> Device::smoothing(vector<pack_input> &input) {
	vector<pack_input> result;

	for (int i = 0; i < input.size()-1; i++) {

		pack_input push;

		push = input[i];

		push.accel_x = (input[i].accel_x + input[i + 1].accel_x) / 2.0;
		push.accel_y = (input[i].accel_y + input[i + 1].accel_y) / 2.0;
		push.accel_z = (input[i].accel_z + input[i + 1].accel_z) / 2.0;

		result.push_back(push);
	}
	result.push_back(input[input.size() - 1]);

	return result;
}

Device::~Device()
{
}
