#include "Device.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "kalman.h"

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
	file_out << count << endl;

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

void Device::determinationOfMeasurementErrors(string input, string inputTrajectory, string inputGPS) {

	vector<pack_output> trajectory;
	vector<pack_input> trajectory_first, trajectory_next;
	pack_output trajectory_start, trajectory_end;
	double count = 0;

	////////////////
	ifstream file(inputTrajectory);
	file >> count;

	for (int i = 0; i < count; i++) {
		pack_output push;

		file >> push.t;
		file >> push.x;
		file >> push.y;
		file >> push.z;

		trajectory.push_back(push);

		if (i == 0)
			trajectory_start = push;
		if ((i + 1) == count)
			trajectory_end = push;
	}

	file.close();
	////////////////

	/////////
	ifstream fileGPS(inputGPS);

	gps gps_start = readGPSPackOfFile(fileGPS);
	gps gps_end = readGPSPackOfFile(fileGPS);

	fileGPS.close();
	/////////

	vecXYZ pXYZ, pXYZ_correct; // 0 0 0 0 0 0

	resGPS resGPS_correct = gettingDifferenceGPS(gps_start, gps_end);

	pXYZ.x = trajectory_end.x - trajectory_start.x;
	pXYZ.y = trajectory_end.y - trajectory_start.y;
	pXYZ.z = trajectory_end.z - trajectory_start.z;

	pXYZ_correct.z = resGPS_correct.h2 - resGPS_correct.h1;

	double m_d2 = sqrt(pow(resGPS_correct.distance,2) - pow(pXYZ_correct.z,2));
	//double m_d2 = sqrt(pow(47.22, 2) - pow(pXYZ_correct.z, 2));

	vecXY pXY;

	pXY.x = pXYZ.x;
	pXY.y = pXYZ.y;

	double modPXY = sqrt(pow(pXY.x, 2) + pow(pXY.y, 2));

	pXY.x /= modPXY;
	pXY.y /= modPXY;

	pXY.x *= m_d2;
	pXY.y *= m_d2;

	pXYZ_correct.x = pXY.x;
	pXYZ_correct.y = pXY.y;

	/////////////////////////////////////
	double x_cor, y_cor, z_cor;

	x_cor = pXYZ_correct.x - pXYZ.x;
	y_cor = pXYZ_correct.y - pXYZ.y;
	z_cor = pXYZ_correct.z - pXYZ.z;

	x_cor /= count;
	y_cor /= count;
	z_cor /= count;

	for (int i = 0; i < count; i++) {
		trajectory[i].x += x_cor * (i + 1);
		trajectory[i].y += y_cor * (i + 1);
		trajectory[i].z += z_cor * (i + 1);
	}

	vector<double> trajectory_x, trajectory_y, trajectory_z;

	//trajectory_x.reserve(count);
	//trajectory_y.reserve(count);
	//trajectory_z.reserve(count);

	trajectory_x.resize(count);
	trajectory_y.resize(count);
	trajectory_z.resize(count);

	for (int i = 0; i < count; i++) {
		//trajectory_x.push_back[trajectory[i].x];
		//trajectory_y.push_back[trajectory[i].y];
		//trajectory_z.push_back[trajectory[i].z];
		trajectory_x[i] = trajectory[i].x;
		trajectory_y[i] = trajectory[i].y;
		trajectory_z[i] = trajectory[i].z;
	}

	trajectory_x = smoothingKalman(trajectory_x);
	trajectory_y = smoothingKalman(trajectory_y);
	trajectory_z = smoothingKalman(trajectory_z);

	for (int i = 0; i < count; i++) {
		trajectory[i].x = trajectory_x[i];
		trajectory[i].y = trajectory_y[i];
		trajectory[i].z = trajectory_z[i];
	}

	/////////////////////////////////////////////////////

	for (int i = 0; i < (count - 1); i++) {
		double deltaTime = trajectory[i + 1].t - trajectory[i].t;
		trajectory[i].v_x = (trajectory[i + 1].x - trajectory[i].x) / deltaTime;
		trajectory[i].v_y = -1 * (trajectory[i + 1].y - trajectory[i].y) / deltaTime;
		trajectory[i].v_z = -1 * (trajectory[i + 1].z - trajectory[i].z) / deltaTime;
	}
	trajectory[count - 1].v_x = trajectory[count - 2].v_x;
	trajectory[count - 1].v_y = trajectory[count - 2].v_y;
	trajectory[count - 1].v_z = trajectory[count - 2].v_z;

	ifstream fileAGV(input);

	fileAGV >> count;

	for (int i = 0; i < count; i++) {
		pack_input push;

		readPackOfFile(fileAGV, push);

		trajectory_first.push_back(push);
	}

	fileAGV.close();

	for (int i = 0; i < count; i++) {
		pack_input push;

		push.t = trajectory_first[i].t;
		push.fi_x = trajectory_first[i].fi_x;
		push.fi_y = trajectory_first[i].fi_y;
		push.fi_z = trajectory_first[i].fi_z;

		trajectory_next.push_back(push);
	}

	for (int i = 0; i < count; i++) {
		pack_output res = rotate(trajectory_next[i], trajectory[i]);

		trajectory[i].v_x = res.v_x;
		trajectory[i].v_y = res.v_y;
		trajectory[i].v_z = res.v_z;
	}

	for (int i = 0; i < (count-1); i++) {
		double deltaTime = trajectory[i + 1].t - trajectory[i].t;
		trajectory_next[i].accel_x = (trajectory[i + 1].v_x - trajectory[i].v_x) / deltaTime;
		trajectory_next[i].accel_y = (trajectory[i + 1].v_y - trajectory[i].v_y) / deltaTime;
		trajectory_next[i].accel_z = (trajectory[i + 1].v_z - trajectory[i].v_z) / deltaTime;
	}
	trajectory_next[count - 1].accel_x = trajectory_next[count - 2].accel_x;
	trajectory_next[count - 1].accel_y = trajectory_next[count - 2].accel_y;
	trajectory_next[count - 1].accel_z = trajectory_next[count - 2].accel_z;

	cout << "Quadratic Deviation: " << quadraticDeviation(x_cor, y_cor, z_cor, count);
	/////////////////////////////////////////////////////

	/////////////////////////////////////
}

void Device::algorithmPathRestoration(vector<pack_input> &input, vector<pack_output> &output, vector<omega> omegaPack, pack_output &state, pack_input &state_input) {
	
	for (int i = 0; i < input.size(); i++) {
		input[i].accel_z -= g;
	}

	input = smoothing(input);

	for (int i = 0; i < input.size(); i++) {

		pack_output outputBuffer;

		double delta_time = input[i].t - state_input.t;

		//state.v_x += delta_time * input[i].accel_x;
		state.v_x = input[i].v_x;
		state.v_y += delta_time * input[i].accel_y;
		state.v_z += delta_time * input[i].accel_z;

		//////////////
		if(i < input.size() - 1)
			subtractionCentrifugalForce(state, input[i], omegaPack[i], delta_time);
		//////////////

		pack_output resultRotate = backRotate(state_input, state);

		////////////
		state.t = input[i].t;
		state.x += delta_time * resultRotate.v_x;
		state.y -= delta_time * resultRotate.v_y;
		state.z -= delta_time * resultRotate.v_z;
		////////////

		outputBuffer = state;

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

gps Device::readGPSPackOfFile(ifstream &inp) {
	gps result;

	string buf;
	double bufDouble = 0;

	inp >> buf;
	inp >> buf;
	inp >> result.latitude;

	result.latitude /= 100;
	result.latitude *= pi / 180.0;

	inp >> buf;
	if (buf.size() == 1)
		if (buf[0] == 'S')
			result.latitude *= -1;

	inp >> result.longitude;

	result.longitude /= 100;
	result.longitude *= pi / 180.0;

	inp >> buf;
	if (buf.size() == 1)
		if (buf[0] == 'W')
			result.longitude *= -1;

	inp >> buf;
	inp >> buf;
	inp >> buf;

	inp >> result.altitude;
	inp >> buf;

	inp >> bufDouble;

	result.altitude -= bufDouble;

	inp >> buf;
	inp >> buf;
	inp >> buf;

	return result;
}

pack_output Device::rotate(pack_input &inputFirst, pack_output &output) {
	pack_output result;

	double alfa = inputFirst.fi_x,
		beta = inputFirst.fi_y,
		gamma = inputFirst.fi_z;

	double v_x = output.v_x * (cos(beta)*cos(gamma))
		+ output.v_y * (-1 * cos(beta)*sin(gamma))
		+ output.v_z * (sin(beta));

	double v_y = output.v_x * (cos(alfa) * sin(gamma) + cos(gamma)*sin(alfa)*sin(beta))
		+ output.v_y * (cos(alfa)*cos(gamma) - sin(alfa)*sin(beta)*sin(gamma))
		+ output.v_z * (-1 * cos(beta)*sin(alfa));

	double v_z = output.v_x * (sin(alfa)*sin(gamma) - cos(alfa)*cos(gamma)*sin(beta))
		+ output.v_y * (cos(gamma)*sin(alfa) + cos(alfa)*sin(beta)*sin(gamma))
		+ output.v_z * (cos(alfa)*cos(beta));

	result.v_x = v_x;
	result.v_y = v_y;
	result.v_z = v_z;

	return result;
}

pack_output Device::backRotate(pack_input &inputFirst, pack_output &output) {
	pack_output result;

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

	result.v_x = v_x;
	result.v_y = v_y;
	result.v_z = v_z;

	return result;
}

vector<pack_input> Device::smoothing(vector<pack_input> &input) {
	vector<pack_input> result;

	vector<double> accel_x, accel_y, accel_z;
	vector<double> fi_x, fi_y, fi_z;

	for (int i = 0; i < input.size(); i++) {
		accel_x.push_back(input[i].accel_x);
		accel_y.push_back(input[i].accel_y);
		accel_z.push_back(input[i].accel_z);

		fi_x.push_back(input[i].fi_x);
		fi_y.push_back(input[i].fi_y);
		fi_z.push_back(input[i].fi_z);
	}

	accel_x = smoothingKalman(accel_x);
	accel_y = smoothingKalman(accel_y);
	accel_z = smoothingKalman(accel_z);

	fi_x = smoothingKalman(fi_x);
	fi_y = smoothingKalman(fi_y);
	fi_z = smoothingKalman(fi_z);

	for (int i = 0; i < input.size(); i++) {
		pack_input push;

		push = input[i];

		push.accel_x = accel_x[i];
		push.accel_y = accel_y[i];
		push.accel_z = accel_z[i];

		push.fi_x = fi_x[i];
		push.fi_y = fi_y[i];
		push.fi_z = fi_z[i];

		result.push_back(push);
	}

	return result;
}

vector<omega> Device::smoothing(vector<omega> &input) {
	vector<omega> result;

	vector<double> omega_x, omega_y, omega_z;

	for (int i = 0; i < input.size(); i++) {
		omega_x.push_back(input[i].omega_x);
		omega_y.push_back(input[i].omega_y);
		omega_z.push_back(input[i].omega_z);
	}

	omega_x = smoothingKalman(omega_x);
	omega_y = smoothingKalman(omega_y);
	omega_z = smoothingKalman(omega_z);

	for (int i = 0; i < input.size(); i++) {
		omega push;

		push.omega_x = omega_x[i];
		push.omega_y = omega_y[i];
		push.omega_z = omega_z[i];

		result.push_back(push);
	}

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

vector<double> Device::smoothingKalman(vector<double> &measurements) {

	vector<double> result;

	////////////////////
	int n = 3; // Number of states
	int m = 1; // Number of measurements

	double dt = 0.01; // Time step

	Eigen::MatrixXd A(n, n); // System dynamics matrix
	Eigen::MatrixXd C(m, n); // Output matrix
	Eigen::MatrixXd Q(n, n); // Process noise covariance
	Eigen::MatrixXd R(m, m); // Measurement noise covariance
	Eigen::MatrixXd P(n, n); // Estimate error covariance

							 // Discrete LTI projectile motion, measuring position only
	A << 1, dt, 0, 0, 1, dt, 0, 0, 1;
	C << 1, 0, 0;

	// Reasonable covariance matrices
	Q << .05, .05, .0, .05, .05, .0, .0, .0, .0;
	R << 5;
	P << .1, .1, .1, .1, 10000, 10, .1, 10, 100;

	// Construct the filter
	KalmanFilter kf(dt, A, C, Q, R, P);

	// Best guess of initial states
	Eigen::VectorXd x0(n);
	x0 << measurements[0], 0, -9.81;
	kf.init(0, x0);

	// Feed measurements into filter, output estimated states
	double t = 0;
	Eigen::VectorXd y(m);
	result.push_back(kf.state().transpose()[0]);
	for (int i = 0; i < measurements.size(); i++) {
		t += dt;
		y << measurements[i];
		kf.update(y);
		result.push_back(kf.state().transpose()[0]);
	}

	return result;
}

double Device::quadraticDeviation(double x_cor, double y_cor, double z_cor, int count) {
	double result = 0.0;

	for (int i = 0; i < count; i++) {
		result += sqrt(pow(x_cor * (i + 1), 2) + pow(y_cor * (i + 1), 2) + pow(z_cor * (i + 1), 2));
	}

	result /= count;

	return result;
}

resGPS Device::gettingDifferenceGPS(gps gps_first, gps gps_next) {
	resGPS result;

	sphericalCS sphericalCS_first;
	sphericalCS sphericalCS_next;

	sphericalCS_first = geographicalToSpherical(gps_first);
	sphericalCS_next = geographicalToSpherical(gps_next);

	result.h1 = gps_first.altitude;
	result.h2 = gps_next.altitude;
	result.distance = sqrt(pow(sphericalCS_next.x - sphericalCS_first.x,2) 
		+ pow(sphericalCS_next.y - sphericalCS_first.y,2)
		+ pow(sphericalCS_next.z - sphericalCS_first.z, 2));

	return result;
}

sphericalCS Device::geographicalToSpherical(gps input) {
	sphericalCS result;

	double h = input.altitude;

	result.x = (r0 + h)*cos(input.latitude)*cos(input.longitude);
	result.y = (r0 + h)*cos(input.latitude)*sin(input.longitude);
	result.z = (r0 + h)*sin(input.latitude);

	return result;
}

Device::~Device()
{
}
