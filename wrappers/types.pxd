cdef extern from '../engine/types.h':
    struct forces_t
    struct velocities_t
    struct xyz_t
    struct xy_t
    ctypedef generator_t
    ctypedef enum LigandCategory:
        psgl
        integrin