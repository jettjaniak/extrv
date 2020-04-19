#include "interpolated.h"

#include "helpers.h"

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>



namespace interpolated {
    double f_r(double h_over_r)
    {
        if (h_over_r >= F_R_DATA_X.front()) {
            double interp = helpers::linear_interpolation(F_R_DATA_X, LOG_F_R_DATA_Y, h_over_r);
            return exp(interp);
        }
        else return -(2.0 / 15.0) * log(h_over_r) - 0.2526;
    }

    double f_s(double h_over_r)
    {
        double interp = helpers::linear_interpolation(F_S_DATA_X, LOG_F_S_DATA_Y_MINUS_1, h_over_r);
        return exp(interp) + 1.0;
    }

    double f_t(double h_over_r)
    {
        if (h_over_r >= F_T_DATA_X.front())
        {
            double interp = helpers::linear_interpolation(F_T_DATA_X, LOG_MINUS_F_T_DATA_Y_MINUS_1, h_over_r);
            return - exp(interp) - 1.0;
        }
        else return (8.0 / 15.0) * log(h_over_r) - 0.9588;
    }

    double t_r(double h_over_r)
    {
        if (h_over_r >= T_R_DATA_X.front()) {
            double interp = helpers::linear_interpolation(T_R_DATA_X, LOG_MINUS_T_R_DATA_Y_MINUS_1, h_over_r);
            return - exp(interp) - 1.0;
        }
        else return 0.4 * log(h_over_r) - 0.3817;
    }

    double t_s(double h_over_r)
    {
        double interp = helpers::linear_interpolation(T_S_DATA_X, LOG_MINUS_T_S_DATA_Y_PLUS_1, h_over_r);
        return - exp(interp) + 1.0;
    }

    double t_t(double h_over_r)
    {
        if (h_over_r >= T_T_DATA_X.front()) {
            double interp = helpers::linear_interpolation(T_T_DATA_X, LOG_T_T_DATA_Y, h_over_r);
            return exp(interp);
        }
        else return - 0.1 * log(h_over_r) - 0.1895;
    }

    double lambda_fun(double h_over_r) {
        const static boost::math::interpolators::cardinal_cubic_b_spline<double>
                lambda_spline(LAMBDA_DATA.begin(), LAMBDA_DATA.end(), LAMBDA_T0, LAMBDA_STEP, LAMBDA_T0_PRIME);
        if (h_over_r < LAMBDA_T0)
            return 1 / h_over_r;
        else
            return lambda_spline(h_over_r);
    }
}