#pragma once

class Parameter {
public:
	static Parameter& getInstance();
	void initParameters();
	double getTimestep();
	double getMass();
	double getGravity();
	double getStretchStiffness();
	double getBendingStiffness();
	int getIteration();

private:
	double timeStep;
	double mass;
	double gravity;
	double stretch;
	double bending;
	int iteration;

	Parameter();
	~Parameter();
	static Parameter* instance;
};