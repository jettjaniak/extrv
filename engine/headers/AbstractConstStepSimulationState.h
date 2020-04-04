#pragma once
#include "AbstractSimulationState.h"


struct AbstractConstStepSimulationState : AbstractSimulationState {

    double dt;

    AbstractConstStepSimulationState(double h_0, Parameters* p, unsigned int seed, double dt);
};