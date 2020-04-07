#include "AbstrGillSS.h"

#include <iostream>


//AbstrGillSS::AbstrGillSS(double h_0, Parameters *p, unsigned int seed) :
//    AbstrSS(h_0, p, seed) {}


double AbstrGillSS::compute_dt_bonds() {
    //////////////////////////////////////////////
    //  Gillespie algorithm - first event time  //
    //////////////////////////////////////////////

    double any_event_rate = 0.0;
    for (const auto &rate : rates)
        any_event_rate += rate;

    // TODO: just for debug, remove
    static constexpr double max_rate = 1e80;
    if (any_event_rate > max_rate)
        std::cout << "EVENT RATE > " << max_rate << std::endl;

    if (any_event_rate > 0.0) {
        std::exponential_distribution<double> dt_bonds_distribution(any_event_rate);
        return dt_bonds_distribution(generator);
    }
    return INFTY;
}


vector<size_t> AbstrGillSS::compute_event_nrs(double) {
    std::discrete_distribution<size_t> which_event_distribution(rates.begin(), rates.end());
    return {which_event_distribution(generator)};
}
