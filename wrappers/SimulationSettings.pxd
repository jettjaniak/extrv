from libcpp.vector cimport vector
from libcpp.pair cimport pair


cdef extern from '../engine/types.h':
    ctypedef enum LigandCategory:
        psgl
        integrin


cdef extern from '../engine/SimulationSettings.h':
    cdef cppclass BondParameters:
        # optimal bond length in μm
        double lambda_
        # spring constant in kg/s^2
        double sigma

        # multiplicative constant in binding rate in μm^2/s, it doesn't include receptor density
        double k_f_0
        # density of bond receptors on the surface in 1/μm^2
        double rec_dens

        # reactive compliance for slip bond in μm
        double x1s
        # multiplicative constant in slip part of binding or rupture rate in 1/s
        double k01s

        # reactive compliance for catch bond in μm
        double x1c
        # multiplicative constant in catch part of rupture rate in 1/s
        double k01c

        BondParameters() except +
        BondParameters(double lambda__, double sigma_, double k_f_0_, double rec_dens_,
            double x1s_, double k01s_, double x1c_, double k01c_) except +

    cdef cppclass LigandParameters:
        vector[BondParameters*] bonds_p
        LigandCategory lig_category

        LigandParameters() except +
        LigandParameters(LigandCategory lig_category_) except +
        void add_bond_p(BondParameters* bond_p) except +

    # Workaround to use it in pair[ , ]
    ctypedef LigandParameters* LigandParameters_ptr

    cdef cppclass Parameters:
        # cell radius in μm
        double r_c
        # viscosity in kg/(μm s)
        double mu
        # temperature in K
        double temp
        # density difference in kg/μm^3
        double dens_diff
    
        # repulsive force coefficient in kg μm / s^2
        double f_rep_0
        # reciprocal length scale of repulsive force in Å
        double tau

        Parameters() except +
        Parameters(double r_c_, double mu_, double temp_, double dens_diff_, double f_rep_0_, double tau_) except +

    cdef cppclass SimulationSettings:
        Parameters* p;
        vector[pair[LigandParameters_ptr, size_t]] lig_types

        SimulationSettings() except +
        SimulationSettings(Parameters* p_) except +
        void add_lig_type(LigandParameters* lig_p, size_t n_of_lig) except +
