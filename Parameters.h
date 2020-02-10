#pragma once

#include "types.h"

struct LigandParameters {
    double lambda_;
    double k_on_0;
    double sigma;
    double x1s;
};

struct Parameters {
    vector<LigandParameters*> lig_p;
};
