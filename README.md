# Installation

1. Install CMake. In Windows make sure installer adds it to PATH.
1. Run terminal / cmd / PowerShell and activate your Python environment (for Anaconda use Anaconda prompt).
1. Go to Git root (this directory).
1. Type `pip install ./engine/`

# Parameters of bond types:

All of them are real numbers (`double` / `float`).  
When passed to constructor, they are converted from {unit} and stored in (unit).

* **eq_bond_len**: equilibrium bond length {nm} (μm)
* **spring_const**: spring constant {dyn/cm} (kg/s^2)
* **binding_rate_0**: multiplicative constant in binding rate, it doesn't include receptor density {μm^2/s} (μm^2/s)
* **rec_dens**: density of bond receptors on the surface {1/μm^2} (1/μm^2)
* **react_compl_slip**: reactive compliance for binding rate and slip part of rupture rate {Å} (μm)
* **rup_rate_0_slip**: multiplicative constant in slip part of rupture rate {1/s} (1/s)
* **react_compl_catch**: reactive compliance for catch part of rupture rate {Å} (μm)
* **rup_rate_0_catch**: multiplicative constant in catch part of rupture rate {1/s} (1/s)

# Parameters of model

Again, all of them are real, { } and ( ) as above.

* **r_cell**: cell (sphere) radius {μm} (μm)
* **visc**: fluid viscosity {g/(cm s)} (kg/(μm s))
* **temp**: temperature {K} (K)
* **dens_diff**: density difference between cell and fluid (cell is more dense) {g/cm^3} (kg/μm^3)
* **f_rep_0**: repulsive force coefficient {N} ((kg μm)/s^2)
* **tau**: reciprocal length scale of repulsive force {Å} (Å)

# Parameters of simulation

All real, no unit conversion.

* **shear_rate**: shear rate in 1/s
* **rot**: sphere's rotation from initial position in radians
* **h**: height above endothelium (surface) in μm