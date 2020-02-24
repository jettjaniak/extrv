from libcpp.vector cimport vector
from Settings cimport Settings

cdef extern from '../engine/headers/SimulationState.h':
    cdef cppclass SimulationState:
        # distance from sphere to surface
        double h
        # sphere's rotation
        double rot

        Settings* settings

        SimulationState() except +
        SimulationState(double h_0, Settings* settings_, unsigned int seed)
    
        void simulate_one_step(double dt, double shear)
        void reseed(unsigned int seed)

    cdef cppclass Stats:
        vector[size_t] n_bd_lig_vec;

        Stats() except +
        Stats(SimulationState* s_) except +

        void update()