#include "ClothSimulator.h"

void ClothSimulator::explicitEuler() {
	double dts = Parameter::getInstance().getTimestep();

	// vertex constraint (horizontal)
	//inv_mass[6] = 0.0001;
	//inv_mass[119] = 0.0001;

	// vertex constraint (diagonal)
	inv_mass[6] = 0.0001;
	inv_mass[166] = 0.0001;

	 // vertex constraint (center)
	//inv_mass[527] = 0.0001;

	//for (int i = 0; i < 100; i++) {
	//	inv_mass[i] = 0.0001;
	//}

	// velocity and predicted position update
	for (int i = 0; i < nV; i++) {
		pred_pos[i] = pos[i] + dts * vel[i] + dts * dts * inv_mass[i] * force[i];

		pos_itr[i] = pred_pos[i];
	}

	fill(lambda_itr.begin(), lambda_itr.end(), 0.0);

}

void ClothSimulator::integrationForPBD() {
	double dts = Parameter::getInstance().getTimestep();

	for (int i = 0; i < nV; i++) {
		vel[i] = (pos_itr[i] - pos[i]) / dts;
		pos[i] = pos_itr[i];
	}
}