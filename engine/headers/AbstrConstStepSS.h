#pragma once
#include "AbstrSS.h"
#include "AbstrGillSS.h"
#include "AbstrProbSS.h"


struct AbstrConstStepSS : virtual AbstrSS {

    double dt;

    explicit AbstrConstStepSS(double dt);
};