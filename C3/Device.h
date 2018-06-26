#pragma once

#include <vector>
#include <string>

#define g 9.81
#define numberPack 10
#define r0 6371000
#define pi 3.14

typedef struct vecXYZ {
	vecXYZ() {
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
	double x,
		y,
		z;
} vec;

typedef struct vecXY {
	vecXY() {
		x = 0.0;
		y = 0.0;
	}
	double x,
		y;
} vecXY;

typedef struct resGPS {
	double h1,
		h2,
		distance;
} resGPS;

typedef struct gps {
	double latitude, //широта
		longitude, //долгота
		altitude; //высота
} gps;

typedef struct sphericalCS {
	double x,
		y,
		z;
} sphericalCS;

typedef struct omega {
	omega() {
		omega_x = 0.0;
		omega_y = 0.0;
		omega_z = 0.0;
	}
	double omega_x,
		omega_y,
		omega_z;
} omega;

typedef struct pack_input {
	pack_input() {
		t = 0.0;
		fi_x = 0.0;
		fi_y = 0.0;
		fi_z = 0.0;
		accel_x = 0.0;
		accel_y = 0.0;
		accel_z = 0.0;
		v_x = 0.0;
	}
	double t,
		fi_x,
		fi_y,
		fi_z,
		accel_x,
		accel_y,
		accel_z,
		v_x;
} pack_input;

typedef struct pack_output {
	pack_output() {
		t = 0.0;
		x = 0.0;
		y = 0.0;
		z = 0.0;
		v_x = 0.0;
		v_y = 0.0;
		v_z = 0.0;
	}
	double t,
		x,
		y,
		z,
		v_x,
		v_y,
		v_z;
} pack_output;

using namespace std;

class Device
{
public:
	Device();
	static void pathRestoration(string inputFile, string outputFile);
	static void determinationOfMeasurementErrors(string input, string inputTrajectory, string inputGPS);
	static void algorithmPathRestoration(vector<pack_input> &input, vector<pack_output> &output, vector<omega> omegaPack, pack_output &state, pack_input &state_input);
	static void writePackToFile(ofstream &out, pack_output &packOut);
	static void readPackOfFile(ifstream &inp, pack_input &packInp);
	static gps readGPSPackOfFile(ifstream &inp);
	static pack_output rotate(pack_input &inputFirst, pack_output &output);
	static pack_output backRotate(pack_input &inputFirst, pack_output &output);
	static vector<pack_input> smoothing(vector<pack_input> &input);
	static vector<omega> smoothing(vector<omega> &input);
	static void subtractionCentrifugalForce(pack_output &state, pack_input &input, omega stateOmega, double deltaTime);
	static vector<double> smoothingKalman(vector<double> &input);
	static double quadraticDeviation(double x_cor, double y_cor, double z_cor, int count);

	static resGPS gettingDifferenceGPS(gps gps_first, gps gps_next);
	static sphericalCS geographicalToSpherical(gps input);
	~Device();
};

