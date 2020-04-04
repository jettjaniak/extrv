#include "AbstractConstStepSimulationState.h"

AbstractConstStepSimulationState::AbstractConstStepSimulationState(
        double h_0, Parameters *p, unsigned int seed, double dt) :

        AbstractSimulationState(h_0, p, seed),
        dt(dt)

{
    try_dt = dt;
}
