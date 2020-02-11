#pragma once

#include "types.h"

struct BondParameters {
    // optimal bond length in μm
    double lambda_ = 0.0;
    // spring constant in TODO: unit
    double sigma = 0.0;

    // reactive compliance for slip bond in μm
    double x1s = 0.0;
    // multiplicative constant in slip part of binding or rupture rate in 1/s
    double k01s = 0.0;

    // reactive compliance for catch bond in μm
    double x1c = 0.0;
    // multiplicative constant in catch part of binding or rupture rate in 1/s
    double k01c = 0.0;

};

struct LigandParameters {
    LigandType lig_type;
    vector<BondParameters*> bonds_p;
};

struct Parameters {
    vector<LigandParameters*> lig_p;
};
