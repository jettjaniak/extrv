#include "AbstrProbSS.h"


double AbstrProbSS::compute_dt_bonds() {
    return -INFTY;
}

vector<size_t> AbstrProbSS::compute_event_nrs(double dt) {
    vector<size_t> event_nrs;
    size_t event_i = 0;
    double uniform, prob;
    for (const auto &rate : rates) {
        uniform = helpers::draw_from_uniform_dist(generator);
        prob = 1 - exp(- dt * rate);
        if (uniform < prob)
            event_nrs.push_back(event_i);
        event_i++;
    }
    return event_nrs;
}
