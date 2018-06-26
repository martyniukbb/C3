#include <iostream>

#include "Device.h"

using namespace std;

int main(void) {

	Device::pathRestoration("in.dat", "out.dat");

	Device::determinationOfMeasurementErrors("in.dat", "out.dat", "inGPS.dat");

	system("PAUSE");

	return 0;
}