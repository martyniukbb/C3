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
		vector<omega> lastOmegaPack;

		for (int j = 0; (j < numberPack) && (i < count); j++, i++) {
			pack_input readBuffer;
			readPackOfFile(file, readBuffer);
			lastReadPack.push_back(readBuffer);
		}

		for (int i = 0; i < lastReadPack.size()-1; i++) {
			omega pushBuffer;
			
			double deltaTime = lastReadPack[i + 1].t - lastReadPack[i].t;

			pushBuffer.omega_x = (lastReadPack[i + 1].fi_x - lastReadPack[i].fi_x)
				/ deltaTime;
			pushBuffer.omega_y = (lastReadPack[i + 1].fi_y - lastReadPack[i].fi_y)
				/ deltaTime;
			pushBuffer.omega_z = (lastReadPack[i + 1].fi_z - lastReadPack[i].fi_z)
				/ deltaTime;

			lastOmegaPack.push_back(pushBuffer);
		}

		lastOmegaPack = smoothing(lastOmegaPack);

		///////////////
		algorithmPathRestoration(lastReadPack, lastOutputPack, lastOmegaPack, state, state_input);
		///////////////

		for (int j = 0; j < int(lastOutputPack.size()); j++) {
			writePackToFile(file_out, lastOutputPack[j]);
		}
	}
	file.close();
	file_out.close();
}

void Device::algorithmPathRestoration(vector<pack_input> &input, vector<pack_output> &output, vector<omega> omegaPack, pack_output &state, pack_input &state_input) {
	
	for (int i = 0; i < input.size(); i++) {
		input[i].accel_z -= g;
	}

	input = smoothing(input);

	for (int i = 0; i < input.size(); i++) {

		pack_output outputBuffer;

		double delta_time = input[i].t - state_input.t;

		//state.v_x = input[i].v_x;

		state.v_x += delta_time * input[i].accel_x;
		state.v_y += delta_time * input[i].accel_y;
		state.v_z += delta_time * input[i].accel_z;

		//////////////
		if(i < input.size() - 1)
			subtractionCentrifugalForce(state, input[i], omegaPack[i], delta_time);
		//////////////

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
	/*double v_x = output.v_x * (cos(alfa)*cos(gamma) - sin(alfa)*cos(beta)*sin(gamma))
		+ output.v_y * (sin(alfa) * cos(gamma) + cos(alfa) * cos(beta) * sin(gamma))
		+ output.v_z * (sin(beta) * sin(gamma));

	double v_y = output.v_x * (-cos(alfa)*sin(gamma) - sin(alfa) * cos(beta)*cos(gamma))
		+ output.v_y * (-sin(alfa) * sin(gamma) + cos(alfa)*cos(beta)*cos(gamma))
		+ output.v_z * (sin(beta) * cos(gamma));

	double v_z = output.v_x * (sin(alfa)*sin(beta))
		+ output.v_y * (-cos(alfa)*sin(beta))
		+ output.v_z * cos(beta);*/

	double v_x = output.v_x * (cos(beta)*cos(gamma))
		+ output.v_y * (cos(alfa) * sin(gamma) + cos(gamma)*sin(alfa)*sin(beta))
		+ output.v_z * (sin(alfa)*sin(gamma) - cos(alfa)*cos(gamma)*sin(beta));

	double v_y = output.v_x * (-1 * cos(beta)*sin(gamma))
		+ output.v_y * (cos(alfa)*cos(gamma) - sin(alfa)*sin(beta)*sin(gamma))
		+ output.v_z * (cos(gamma)*sin(alfa) + cos(alfa)*sin(beta)*sin(gamma));

	double v_z = output.v_x * (sin(beta))
		+ output.v_y * (-1 * cos(beta)*sin(alfa))
		+ output.v_z * (cos(alfa)*cos(beta));

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

vector<omega> Device::smoothing(vector<omega> &input) {
	vector<omega> result;

	for (int i = 0; i < input.size() - 1; i++) {

		omega push;

		push = input[i];

		push.omega_x = (input[i].omega_x + input[i + 1].omega_x) / 2.0;
		push.omega_y = (input[i].omega_y + input[i + 1].omega_y) / 2.0;
		push.omega_z = (input[i].omega_z + input[i + 1].omega_z) / 2.0;

		result.push_back(push);
	}
	result.push_back(input[input.size() - 1]);

	return result;
}

void Device::subtractionCentrifugalForce(pack_output &state, pack_input &input, omega stateOmega, double deltaTime) {
	double vx = stateOmega.omega_y * state.v_z - stateOmega.omega_z * state.v_y,
		vy = -1 * (stateOmega.omega_x * state.v_z - stateOmega.omega_z * state.v_x),
		vz = stateOmega.omega_x * state.v_y - stateOmega.omega_y * state.v_x;
	state.v_x -= deltaTime * vx;
	state.v_y -= deltaTime * vy;
	state.v_z -= deltaTime * vz;
}

Device::~Device()
{
}
