#include "helpers.h"


namespace helpers {
    pair<double, double> parametrize_ligand(double lig_x, double lig_y) {
        double r_cir = sqrt(lig_x * lig_x + lig_y * lig_y);
        double alpha;
        // x = r sin(alpha)
        if (lig_y < 0)
            // We are in the bottom half, where alpha is in [-π/2, π/2],
            // which corresponds to asin return value range.
            alpha = asin(lig_x / r_cir);
        else {
            // We are in the top half, where alpha is in
            // [-π, -π/2] or [π/2, π].

            // sin(alpha - π) = - x / r
            alpha = PI - asin(lig_x / r_cir);
            // but we want alpha to stay in [-π, π]
            if (alpha > PI)
                alpha -= 2.0 * PI;
        }
        return {r_cir, alpha};
    }
}