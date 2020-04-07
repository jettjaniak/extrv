#pragma once

#include "AbstrSS.h"

struct AbstrProbSS : virtual AbstrSS {
    AbstrProbSS() = default;

    /// Constant
    double compute_dt_bonds() override;

    /// Compute probability and draw from binary distribution
    vector<size_t> compute_event_nrs(double dt) override;
};