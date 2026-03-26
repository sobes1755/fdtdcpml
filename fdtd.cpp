#include "fdtdcpml.h"

#include <iostream>
#include <cmath>
#include <ranges>
#include <algorithm>
#include <execution>

// "Premature optimization is the root of all evil". Donald Knuth.

void FDTDCPML::fdtd_init_H()
{

    double txyz = (t_step_ * x_step_) / (y_step_ * z_step_);
    double tyzx = (t_step_ * y_step_) / (z_step_ * x_step_);
    double tzxy = (t_step_ * z_step_) / (x_step_ * y_step_);

    // Hx grid dims are p_, q_ - 1, r_ - 1...

    auto gridHx = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_),
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_ - 1));

    std::for_each(std::execution::par, gridHx.begin(), gridHx.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;
        auto [x, y, z] = Hx_pqr2xyz({p, q, r});

        double s_m = sigma_m(x, y, z);
        double mu = fdtdcpml::mu0 * mu_r(x, y, z);
        double tmp = (s_m * t_step_) / (2. * mu);

        HxHx_[p, q, r] = (1 - tmp) / (1. + tmp);
        HxEyEz_[p, q, r] = txyz / mu / (1. + tmp);

    });

    // Hy grid dims are p_ - 1, q_, r_ - 1...

    auto gridHy = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1),
        std::views::iota((size_t) 0, q_),
        std::views::iota((size_t) 0, r_ - 1));

    std::for_each(std::execution::par, gridHy.begin(), gridHy.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;
        auto [x, y, z] = Hy_pqr2xyz({p, q, r});

        double s_m = sigma_m(x, y, z);
        double mu = fdtdcpml::mu0 * mu_r(x, y, z);
        double tmp = (s_m * t_step_) / (2. * mu);

        HyHy_[p, q, r] = (1 - tmp) / (1. + tmp);
        HyEzEx_[p, q, r] = tyzx / mu / (1. + tmp);

    });

    // Hz grid dims are p_ - 1, q_ - 1, r_...

    auto gridHz = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1),
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_));

    std::for_each(std::execution::par, gridHz.begin(), gridHz.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;
        auto [x, y, z] = Hz_pqr2xyz({p, q, r});

        double s_m = sigma_m(x, y, z);
        double mu = fdtdcpml::mu0 * mu_r(x, y, z);
        double tmp = (s_m * t_step_) / (2. * mu);

        HzHz_[p, q, r] = (1 - tmp) / (1. + tmp);
        HzExEy_[p, q, r] = tzxy / mu / (1. + tmp);

    });

}

void FDTDCPML::fdtd_init_E()
{

    double txyz = (t_step_ * x_step_) / (y_step_ * z_step_);
    double tyzx = (t_step_ * y_step_) / (z_step_ * x_step_);
    double tzxy = (t_step_ * z_step_) / (x_step_ * y_step_);

    // Ex grid dims are p_ - 1, q_, r_, but Ex is zero on Y and Z ends (PEC)...

    auto gridEx = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_),
        std::views::iota((size_t) 1, q_ - 1),
        std::views::iota((size_t) 1, r_ - 1));

    std::for_each(std::execution::par, gridEx.begin(), gridEx.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;
        auto [x, y, z] = Ex_pqr2xyz({p, q, r});

        double s = sigma(x, y, z);
        double eps = fdtdcpml::eps0 * eps_r(x, y, z);
        double tmp = (s * t_step_) / (2 * eps);

        ExEx_[p, q, r] = (1 - tmp) / (1 + tmp);
        ExHyHz_[p, q, r] = txyz / eps / (1 + tmp);

    });

    // Ey grid dims are p_, q_ - 1, r_, but Ey is zero on Z and X ends (PEC)...

    auto gridEy = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1),
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 1, r_ - 1));

    std::for_each(std::execution::par, gridEy.begin(), gridEy.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;
        auto [x, y, z] = Ey_pqr2xyz({p, q, r});

        double s = sigma(x, y, z);
        double eps = fdtdcpml::eps0 * eps_r(x, y, z);
        double tmp = (s * t_step_) / (2. * eps);

        EyEy_[p, q, r] = (1 - tmp) / (1. + tmp);
        EyHzHx_[p, q, r] = tyzx / eps / (1. + tmp);

    });

    // Ez grid dims are p_, q_, r_ - 1, but Ez is zero on X and Y ends (PEC)...

    auto gridEz = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1),
        std::views::iota((size_t) 1, q_ - 1),
        std::views::iota((size_t) 0, r_ - 1));

    std::for_each(std::execution::par, gridEz.begin(), gridEz.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;
        auto [x, y, z] = Ez_pqr2xyz({p, q, r});

        double s = sigma(x, y, z);
        double eps = fdtdcpml::eps0 * eps_r(x, y, z);
        double tmp = (s * t_step_) / (2. * eps);

        EzEz_[p, q, r] = (1 - tmp) / (1. + tmp);
        EzHxHy_[p, q, r] = tzxy / eps / (1. + tmp);

    });

}

void FDTDCPML::fdtd_H()
{

    // Hx grid dims are p_, q_ - 1, r_ - 1...

    auto gridHx = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_),
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_ - 1));

    std::for_each(std::execution::par, gridHx.begin(), gridHx.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;

        Hx_[p, q, r] = HxHx_[p, q, r] * Hx_[p, q, r] +
            HxEyEz_[p, q, r] *
            ((Ey_[p, q, r + 1] - Ey_[p, q, r]) -
            (Ez_[p, q + 1, r ] - Ez_[p, q, r]));

    });

    // Hy grid dims are p_ - 1, q_, r_ - 1...

    auto gridHy = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1),
        std::views::iota((size_t) 0, q_),
        std::views::iota((size_t) 0, r_ - 1));

    std::for_each(std::execution::par, gridHy.begin(), gridHy.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;

        Hy_[p, q, r] = HyHy_[p, q, r] * Hy_[p, q, r] +
            HyEzEx_[p, q, r] *
            ((Ez_[p + 1, q, r] - Ez_[p, q, r]) -
            (Ex_[p, q, r + 1] - Ex_[p, q, r]));

    });

    // Hz grid dims are p_ - 1, q_ - 1, r_...

    auto gridHz = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1),
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_));

    std::for_each(std::execution::par, gridHz.begin(), gridHz.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;

        Hz_[p, q, r] = HzHz_[p, q, r] * Hz_[p, q, r] +
            HzExEy_[p, q, r] *
            ((Ex_[p, q + 1, r] - Ex_[p, q, r]) -
            (Ey_[p + 1, q, r] - Ey_[p, q, r]));

    });

}

void FDTDCPML::fdtd_E()
{

    // Ex grid dims are p_ - 1, q_, r_, but Ex is zero on Y and Z ends (PEC)...

    auto gridEx = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1),
        std::views::iota((size_t) 1, q_ - 1),
        std::views::iota((size_t) 1, r_ - 1));

    std::for_each(std::execution::par, gridEx.begin(), gridEx.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;

        Ex_[p, q, r] = ExEx_[p, q, r] * Ex_[p, q, r] +
            ExHyHz_[p, q, r] *
            ((Hz_[p, q, r] - Hz_[p, q - 1, r]) -
            (Hy_[p, q, r] - Hy_[p, q, r - 1]));

    });
  
    // Ey grid dims are p_, q_ - 1, r_, but Ey is zero on Z and X ends (PEC)...

    auto gridEy = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1),
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 1, r_ - 1));

    std::for_each(std::execution::par, gridEy.begin(), gridEy.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;

        Ey_[p, q, r] = EyEy_[p, q, r] * Ey_[p, q, r] +
            EyHzHx_[p, q, r] *
            ((Hx_[p, q, r] - Hx_[p, q, r - 1]) -
            (Hz_[p, q, r] - Hz_[p - 1, q, r]));

    });
  
    // Ez grid dims are p_, q_, r_ - 1, but Ez is zero on X and Y ends (PEC)...

    auto gridEz = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1),
        std::views::iota((size_t) 1, q_ - 1),
        std::views::iota((size_t) 0, r_ - 1));

    std::for_each(std::execution::par, gridEz.begin(), gridEz.end(), [&](auto tuple) {

        auto [p, q, r] = tuple;

        Ez_[p, q, r] = EzEz_[p, q, r] * Ez_[p, q, r] +
            EzHxHy_[p, q, r] *
            ((Hy_[p, q, r] - Hy_[p - 1, q, r]) -
            (Hx_[p, q, r] - Hx_[p, q - 1, r]));

    });

}                     
