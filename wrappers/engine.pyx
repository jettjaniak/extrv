from libc cimport limits
cimport cython
import numpy as np
cimport numpy as np
from collections import namedtuple

from types cimport psgl, integrin  # LigandCategory enum values
from SimulationSettings cimport (
    BondParameters as BondParametersCpp,
    LigandParameters as LigandParametersCpp,
    Parameters as ParametersCpp,
    SimulationSettings as SimulationSettingsCpp,
)
from SimulationState cimport SimulationState as SimulationStateCpp

SimulationResult = namedtuple('SimulationResult', ['h', 'rot', 'n_bonds'])

# TODO: getters, setters
# TODO: documentation

cdef class BondParameters:
    cdef BondParametersCpp _bond_p_cpp

    def __init__(self, double lambda__, double sigma_, double k_f_0_, double rec_dens_,
                 double x1s_, double k01s_, double x1c_ = 0.0, double k01c_ = 0.0):
        self._bond_p_cpp = BondParametersCpp(lambda__, sigma_, k_f_0_, rec_dens_, x1s_, k01s_, x1c_ , k01c_)


cdef class LigandParameters:
    cdef LigandParametersCpp _lig_p_cpp
    cdef object bonds_p

    def __init__(self, str lig_category):
        if lig_category == 'psgl':
            self._lig_p_cpp = LigandParametersCpp(psgl)
        elif lig_category == 'integrin':
            self._lig_p_cpp = LigandParametersCpp(integrin)
        else:
            raise ValueError('Wrong ligand category!')
        self.bonds_p = []

    def add_bond_p(self, BondParameters bond_p):
        self._lig_p_cpp.add_bond_p(&bond_p._bond_p_cpp)
        self.bonds_p.append(bond_p)


cdef class Parameters:
    cdef ParametersCpp _p_cpp

    def __init__(self, double r_c_, double mu_, double temp_, double dens_diff_, double f_rep_0_, double tau_):
        self._p_cpp = ParametersCpp(r_c_, mu_, temp_, dens_diff_, f_rep_0_, tau_)


cdef class SimulationSettings:
    cdef SimulationSettingsCpp _settings_cpp
    cdef object p
    cdef object lig_types

    def __init__(self, Parameters p):
        self._settings_cpp = SimulationSettingsCpp(&(p._p_cpp))
        self.p = p
        self.lig_types = []

    def add_lig_type(self, LigandParameters lig_p, size_t n_of_lig):
        self._settings_cpp.add_lig_type(&lig_p._lig_p_cpp, n_of_lig)
        self.lig_types.append((lig_p, n_of_lig))


cdef class SimulationState:
    cdef SimulationStateCpp _ss_cpp

    def __init__(self, double h_0, SimulationSettings settings, unsigned int seed = 0):
        if seed == 0:
            seed = np.random.randint(limits.UINT_MAX, dtype='uint')
        self._ss_cpp = SimulationStateCpp(h_0, &settings._settings_cpp, seed)

    def simulate(self, size_t n_steps, double dt, double shear):
        cdef size_t i
        for i in range(n_steps):
            self._ss_cpp.simulate_one_step(dt, shear)

    @cython.cdivision(True)
    def simulate_with_history(self, size_t n_steps, double dt, double shear, int save_every):
        cdef size_t max_hist_size = 1 + (n_steps - 1) // save_every
        cdef int n_lig_types = self._ss_cpp.settings.lig_types.size()  # TODO: define it somewhere deeper
        cdef np.ndarray[np.int32_t, ndim=2] n_bonds = np.empty((max_hist_size, n_lig_types), dtype=np.int32)
        cdef np.ndarray[np.double_t, ndim=1] h = np.empty(max_hist_size, dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] rot = np.empty(max_hist_size, dtype=np.double)

        cdef hist_i = 0
        cdef size_t i
        for i in range(n_steps):
            self._ss_cpp.simulate_one_step(dt, shear)
            if i % save_every == 0:
                n_bonds[hist_i] = self._ss_cpp.stats.n_bd_lig_vec
                h[hist_i] = self._ss_cpp.h
                rot[hist_i] = self._ss_cpp.alpha_0
                hist_i += 1
        return SimulationResult(h, rot, n_bonds)

    @property
    def h(self):
        return self._ss_cpp.h

    @property
    def alpha_0(self):
        return self._ss_cpp.alpha_0