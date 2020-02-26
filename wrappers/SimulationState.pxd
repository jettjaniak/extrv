from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.set cimport set

from Settings cimport Settings
from types cimport xy_t

cdef extern from '../engine/headers/SimulationState.h':
    cdef cppclass SimulationState:
        # distance from sphere to surface
        double h
        # sphere's rotation
        double rot

        # TODO: documentation
        set[size_t] bd_lig_ind

        Settings* settings

        SimulationState() except +
        SimulationState(double h_0, Settings* settings_, unsigned int seed)
    
        void simulate_one_step(double dt, double shear)
        void reseed(unsigned int seed)

    cdef cppclass Stats:
        vector[size_t] n_bd_lig_vec;
        vector[pair[size_t, xy_t]] bd_ligs_ind_and_xy;

        Stats() except +
        Stats(SimulationState* s_) except +

        void update()