#include "ClothSimulator.h"

void ClothSimulator::explicitEuler() {
	double dts = Parameter::getInstance().getTimestep();

	// Store currents velocity and position
	for (int i = 0; i < nV; i++) {
		prev_pos[i] = pos[i];
		prev_vel[i] = vel[i];
	}

	// vertex constraint (horizontal)
	//inv_mass[6] = 0.0001;
	//inv_mass[119] = 0.0001;

	// vertex constraint (diagonal)
	//inv_mass[6] = 0.0001;
	//inv_mass[166] = 0.0001;

	 // vertex constraint (center)
	inv_mass[527] = 0.0001;

	//for (int i = 0; i < 100; i++) {
	//	inv_mass[i] = 0.0001;
	//}

	//inv_mass[100] = 0.0001;

	// velocity and predicted position update
	for (int i = 0; i < nV; i++) {
		vel[i] = prev_vel[i] + dts * inv_mass[i] * force[i];
		pred_pos[i] = prev_pos[i] + dts * vel[i];
	}
}

void ClothSimulator::integrationForPBD() {
	double dts = Parameter::getInstance().getTimestep();

	for (int i = 0; i < nV; i++) {
		vel[i] = (pred_pos[i] - pos[i]) / dts;
		pos[i] = pred_pos[i];
	}
}