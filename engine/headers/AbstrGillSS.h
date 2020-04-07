#pragma once

#include "AbstrSS.h"

struct AbstrGillSS : virtual AbstrSS {
    AbstrGillSS() = default;

    /// Gillespie algorithm - first event time
    double compute_dt_bonds() override;

    // Gillespie algorithm - which event first
    vector<size_t> compute_event_nrs(double) override;
};