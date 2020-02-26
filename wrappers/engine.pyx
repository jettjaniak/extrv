from libc cimport limits
cimport cython
from cython.operator cimport dereference, postincrement
import numpy as np
cimport numpy as np
from collections import namedtuple

from types cimport PSEL_BOND, ESEL_BOND, INTEGRIN_BOND, xy_t
from Settings cimport (
    BondParameters as BondParametersCpp,
    LigandType as LigandTypeCpp,
    ModelParameters as ModelParametersCpp,
    Settings as SettingsCpp,
)
from SimulationState cimport SimulationState as SimulationStateCpp, Stats

SimulationResult = namedtuple('SimulationResult', ['h', 'rot'])#, 'n_bonds', 'bd_ligs_ind', 'bd_ligs_xy'])

# TODO: getters, setters
# TODO: documentation

cdef class BondParameters:
    cdef BondParametersCpp _bond_p_cpp

    def __init__(self, str bond_type_str, double lambda__, double sigma_, double k_f_0_, double rec_dens_,
                 double x1s_, double k01s_, double x1c_ = 0.0, double k01c_ = 0.0):
        # TODO: dict
        if bond_type_str == "esel":
            self._bond_p_cpp = BondParametersCpp(ESEL_BOND,
                                                 lambda__, sigma_, k_f_0_, rec_dens_, x1s_, k01s_, x1c_ , k01c_)
        elif bond_type_str == "psel":
            self._bond_p_cpp = BondParametersCpp(PSEL_BOND,
                                                 lambda__, sigma_, k_f_0_, rec_dens_, x1s_, k01s_, x1c_ , k01c_)
        elif bond_type_str == "integrin":
            self._bond_p_cpp = BondParametersCpp(INTEGRIN_BOND,
                                                 lambda__, sigma_, k_f_0_, rec_dens_, x1s_, k01s_, x1c_ , k01c_)
        else:
            raise ValueError(f"`{bond_type_str}` is not correct bond type name")


cdef class LigandType:
    cdef LigandTypeCpp _lig_type_cpp
    cdef object bonds_p

    def __init__(self):
        self._lig_type_cpp = LigandTypeCpp()
        self.bonds_p = []

    def add_bond_p(self, BondParameters bond_p):
        self._lig_type_cpp.add_bond_p(&bond_p._bond_p_cpp)
        self.bonds_p.append(bond_p)


cdef class ModelParameters:
    cdef ModelParametersCpp _p_cpp

    def __init__(self, double r_c_, double mu_, double temp_, double dens_diff_, double f_rep_0_, double tau_):
        self._p_cpp = ModelParametersCpp(r_c_, mu_, temp_, dens_diff_, f_rep_0_, tau_)


cdef class Settings:
    cdef SettingsCpp _settings_cpp
    cdef object p
    cdef object lig_types_and_nrs

    def __init__(self, ModelParameters p):
        self._settings_cpp = SettingsCpp(&(p._p_cpp))
        self.p = p
        self.lig_types_and_nrs = []

    def add_lig_type(self, LigandType lig_type, size_t n_of_lig):
        self._settings_cpp.add_lig_type(&lig_type._lig_type_cpp, n_of_lig)
        self.lig_types_and_nrs.append((lig_type, n_of_lig))

# TODO: extract to separate *.pyx
cdef class SimulationState:
    cdef SimulationStateCpp _ss_cpp

    def __init__(self, double h_0, Settings settings, unsigned int seed = 0):
        if seed == 0:
            seed = np.random.randint(limits.UINT_MAX, dtype='uint')
        self._ss_cpp = SimulationStateCpp(h_0, &settings._settings_cpp, seed)

    def simulate(self, size_t n_steps, double dt, double shear):
        cdef size_t i
        for i in range(n_steps):
            self._ss_cpp.simulate_one_step(dt, shear)

    # TODO: progress bar (tqdm too slow): https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    @cython.cdivision(True)
    def simulate_with_history(self, size_t n_steps, double dt, double shear, int save_every):
        cdef size_t max_hist_size = 1 + (n_steps - 1) // save_every

        cdef np.ndarray[np.double_t, ndim=1] h = np.empty(max_hist_size, dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] rot = np.empty(max_hist_size, dtype=np.double)

        cdef size_t hist_i = 0
        cdef size_t i
        cdef int n_all_bonds, j

        for i in range(n_steps):
            self._ss_cpp.simulate_one_step(dt, shear)
            if i % save_every == 0:
                h[hist_i] = self._ss_cpp.h
                rot[hist_i] = self._ss_cpp.rot
                hist_i += 1
                # TODO: expand, wrap and use History from SimulationState.h
        return SimulationResult(h, rot)

    @property
    def h(self):
        return self._ss_cpp.h

    @property
    def rot(self):
        return self._ss_cpp.rot