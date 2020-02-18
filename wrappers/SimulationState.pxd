from libcpp.set cimport set
from SimulationSettings cimport SimulationSettings

cdef extern from '../engine/SimulationState.h':
    cdef cppclass SimulationState:
        # distance from sphere to surface
        double h
        # sphere's rotation
        double alpha_0

        SimulationSettings* settings

        SimulationState() except +
        SimulationState(double h_0, SimulationSettings* settings_, unsigned int seed)
    
        void simulate_one_step(double dt, double shear)
    
        void reseed(unsigned int seed)