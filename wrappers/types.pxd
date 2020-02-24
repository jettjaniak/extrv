cdef extern from '../engine/headers/types.h':
    struct forces_t
    struct velocities_t
    struct xyz_t
    struct xy_t
    ctypedef generator_t
    ctypedef enum BondType:
        PSEL_BOND
        ESEL_BOND
        INTEGRIN_BOND