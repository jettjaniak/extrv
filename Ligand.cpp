#include "Ligand.h"

#include "helpers.h"


Ligand::Ligand(double lig_x, double lig_y, LigandParameters* lig_p_) {
    auto r_alpha_pair = helpers::parametrize_ligand(lig_x, lig_y);
    r_cir = r_alpha_pair.first;
    alpha_inc = r_alpha_pair.second;

    lig_p = lig_p_;
}