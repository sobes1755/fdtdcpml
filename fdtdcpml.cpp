#include "fdtdcpml.h"

#include <iostream>
#include <cmath>
#include <ranges>
#include <algorithm>
#include <execution>
#include <print>

// "Premature optimization is the root of all evil". Donald Knuth.

void FDTDCPML::init()
{

    fdtd_init_H();
    fdtd_init_E();

    t_ = t_min_;    

}

void FDTDCPML::step() {

    fdtd_H(); 

//    H_source();

//    cpml_x_l_H(); cpml_x_h_H();
//    cpml_y_l_H(); cpml_y_h_H();
//    cpml_z_l_H(); cpml_z_h_H();

    fdtd_E();

//    E_source();

//    cpml_x_l_E(); cpml_x_h_E();
//    cpml_y_l_E(); cpml_y_h_E();
//    cpml_z_l_E(); cpml_z_h_E();

    t_ += t_step_;

    #ifdef DEBUG

    H_minmax();
    E_minmax();

    std::print("t = {},\n\n", t_);

    std::print("Hx_min = {} ({}, {}, {}), Hx_max = {} ({}, {}, {}),\n",
        Hx_min, Hx_min_pqr.p, Hx_min_pqr.q, Hx_min_pqr.r,
        Hx_max, Hx_max_pqr.p, Hx_max_pqr.q, Hx_max_pqr.r);
    std::print("Hy_min = {} ({}, {}, {}), Hy_max = {} ({}, {}, {}),\n",
        Hy_min, Hy_min_pqr.p, Hy_min_pqr.q, Hy_min_pqr.r,
        Hy_max, Hy_max_pqr.p, Hy_max_pqr.q, Hy_max_pqr.r);
    std::print("Hz_min = {} ({}, {}, {}), Hz_max = {} ({}, {}, {}),\n",
        Hz_min, Hz_min_pqr.p, Hz_min_pqr.q, Hz_min_pqr.r,
        Hz_max, Hz_max_pqr.p, Hz_max_pqr.q, Hz_max_pqr.r);

    std::print("Ex_min = {} ({}, {}, {}), Ex_max = {} ({}, {}, {}),\n",
        Ex_min, Ex_min_pqr.p, Ex_min_pqr.q, Ex_min_pqr.r,
        Ex_max, Ex_max_pqr.p, Ex_max_pqr.q, Ex_max_pqr.r);
    std::print("Ey_min = {} ({}, {}, {}), Ey_max = {} ({}, {}, {}),\n",
        Ey_min, Ey_min_pqr.p, Ey_min_pqr.q, Ey_min_pqr.r,
        Ey_max, Ey_max_pqr.p, Ey_max_pqr.q, Ey_max_pqr.r);
    std::print("Ez_min = {} ({}, {}, {}), Ez_max = {} ({}, {}, {}).\n\n",
        Ez_min, Ez_min_pqr.p, Ez_min_pqr.q, Ez_min_pqr.r,
        Ez_max, Ez_max_pqr.p, Ez_max_pqr.q, Ez_max_pqr.r);

    #endif

}

FDTDCPML::FDTDCPML(
    double x_min, double x_max, std::size_t p_int,
    double y_min, double y_max, std::size_t q_int,
    double z_min, double z_max, std::size_t r_int,
    double t_min, double t_max, std::size_t s_int,
    std::size_t p_abs_l, std::size_t p_abs_h,
    std::size_t q_abs_l, std::size_t q_abs_h,
    std::size_t r_abs_l, std::size_t r_abs_h) :

    p_int_{p_int},
    q_int_{q_int},
    r_int_{r_int},
    s_int_{s_int},

    p_abs_l_{p_abs_l}, p_abs_h_{p_abs_h},
    q_abs_l_{q_abs_l}, q_abs_h_{q_abs_h},
    r_abs_l_{r_abs_l}, r_abs_h_{r_abs_h},

    p_{1 + p_abs_l + p_int + p_abs_h + 1},  // The ones here...
    q_{1 + q_abs_l + q_int + q_abs_h + 1},  // are for the PEC (PMC)...
    r_{1 + r_abs_l + r_int + r_abs_h + 1},  // zeroes...

    x_min_{x_min}, x_max_{x_max},
    y_min_{y_min}, y_max_{y_max},
    z_min_{z_min}, z_max_{z_max},
    t_min_{t_min}, t_max_{t_max},

    x_step_{(x_max_ - x_min_) / (p_int_ - 1)},
    y_step_{(y_max_ - y_min_) / (q_int_ - 1)},
    z_step_{(z_max_ - z_min_) / (r_int_ - 1)},
    t_step_{std::min((t_max_ - t_min_) / (s_int_ - 1), 0.5 * t_step_default(x_step_, y_step_, z_step_))},

    rawHx_    (p_ * (q_ - 1) * (r_ - 1)),
    rawHxHx_  (p_ * (q_ - 1) * (r_ - 1)),
    rawHxEyEz_(p_ * (q_ - 1) * (r_ - 1)),
    rawHy_    ((p_ - 1) * q_ * (r_ - 1)),
    rawHyHy_  ((p_ - 1) * q_ * (r_ - 1)),
    rawHyEzEx_((p_ - 1) * q_ * (r_ - 1)),
    rawHz_    ((p_ - 1) * (q_ - 1) * r_),
    rawHzHz_  ((p_ - 1) * (q_ - 1) * r_),
    rawHzExEy_((p_ - 1) * (q_ - 1) * r_),

    rawEx_    ((p_ - 1) * q_ * r_),
    rawExEx_  ((p_ - 1) * q_ * r_),
    rawExHyHz_((p_ - 1) * q_ * r_),
    rawEy_    (p_ * (q_ - 1) * r_),
    rawEyEy_  (p_ * (q_ - 1) * r_),
    rawEyHzHx_(p_ * (q_ - 1) * r_),
    rawEz_    (p_ * q_ * (r_ - 1)),
    rawEzEz_  (p_ * q_ * (r_ - 1)),
    rawEzHxHy_(p_ * q_ * (r_ - 1)),

    Hx_    (rawHx_.data(),     p_, q_ - 1, r_ - 1),
    HxHx_  (rawHxHx_.data(),   p_, q_ - 1, r_ - 1),
    HxEyEz_(rawHxEyEz_.data(), p_, q_ - 1, r_ - 1),
    Hy_    (rawHy_.data(),     p_ - 1, q_, r_ - 1),
    HyHy_  (rawHyHy_.data(),   p_ - 1, q_, r_ - 1),
    HyEzEx_(rawHyEzEx_.data(), p_ - 1, q_, r_ - 1),
    Hz_    (rawHz_.data(),     p_ - 1, q_ - 1, r_),
    HzHz_  (rawHzHz_.data(),   p_ - 1, q_ - 1, r_),
    HzExEy_(rawHzExEy_.data(), p_ - 1, q_ - 1, r_),

    Ex_    (rawEx_.data(),     p_ - 1, q_, r_),
    ExEx_  (rawExEx_.data(),   p_ - 1, q_, r_),
    ExHyHz_(rawExHyHz_.data(), p_ - 1, q_, r_),
    Ey_    (rawEy_.data(),     p_, q_ - 1, r_),
    EyEy_  (rawEyEy_.data(),   p_, q_ - 1, r_),
    EyHzHx_(rawEyHzHx_.data(), p_, q_ - 1, r_),
    Ez_    (rawEz_.data(),     p_, q_, r_ - 1),
    EzEz_  (rawEzEz_.data(),   p_, q_, r_ - 1),
    EzHxHy_(rawEzHxHy_.data(), p_, q_, r_ - 1),

    psi_x_l_(q_ - 1, r_ - 1, p_abs_l_, x_step_, t_step_),
    psi_x_h_(q_ - 1, r_ - 1, p_abs_h_, x_step_, t_step_),
    psi_y_l_(p_ - 1, q_abs_l_, r_ - 1, y_step_, t_step_),
    psi_y_h_(p_ - 1, q_abs_h_, r_ - 1, y_step_, t_step_),
    psi_z_l_(p_ - 1, q_ - 1, r_abs_l_, z_step_, t_step_),
    psi_z_h_(p_ - 1, q_ - 1, r_abs_h_, z_step_, t_step_)

{

    #ifdef DEBUG

    std::print("x_min = {}, x_max = {}, y_min = {}, y_max = {}, z_min = {}, z_max = {};\n\n",
        x_min_, x_max_, y_min_, y_max_, z_min_, z_max_);

    std::print("p_ = {} = 1 + {} + {} + {} + 1,\n",
        p_, p_abs_l_, p_int_, p_abs_h_);
    std::print("q_ = {} = 1 + {} + {} + {} + 1,\n",
        q_, q_abs_l_, q_int_, q_abs_h_);
    std::print("r_ = {} = 1 + {} + {} + {} + 1;\n\n",
        r_, r_abs_l_, r_int_, r_abs_h_);

    std::print("x_step = {}, y_step = {}, z_step = {}, t_step = {},\n\n",
        x_step_, y_step_, z_step_, t_step_);

    H_minmax();
    E_minmax();

    std::print("t = {}, t_min = {}, t_max = {},\n\n", t_, t_min_, t_max_);

    std::print("Hx_min = {} ({}, {}, {}), Hx_max = {} ({}, {}, {}),\n",
        Hx_min, Hx_min_pqr.p, Hx_min_pqr.q, Hx_min_pqr.r,
        Hx_max, Hx_max_pqr.p, Hx_max_pqr.q, Hx_max_pqr.r);
    std::print("Hy_min = {} ({}, {}, {}), Hy_max = {} ({}, {}, {}),\n",
        Hy_min, Hy_min_pqr.p, Hy_min_pqr.q, Hy_min_pqr.r,
        Hy_max, Hy_max_pqr.p, Hy_max_pqr.q, Hy_max_pqr.r);
    std::print("Hz_min = {} ({}, {}, {}), Hz_max = {} ({}, {}, {}),\n",
        Hz_min, Hz_min_pqr.p, Hz_min_pqr.q, Hz_min_pqr.r,
        Hz_max, Hz_max_pqr.p, Hz_max_pqr.q, Hz_max_pqr.r);

    std::print("Ex_min = {} ({}, {}, {}), Ex_max = {} ({}, {}, {}),\n",
        Ex_min, Ex_min_pqr.p, Ex_min_pqr.q, Ex_min_pqr.r,
        Ex_max, Ex_max_pqr.p, Ex_max_pqr.q, Ex_max_pqr.r);
    std::print("Ey_min = {} ({}, {}, {}), Ey_max = {} ({}, {}, {}),\n",
        Ey_min, Ey_min_pqr.p, Ey_min_pqr.q, Ey_min_pqr.r,
        Ey_max, Ey_max_pqr.p, Ey_max_pqr.q, Ey_max_pqr.r);
    std::print("Ez_min = {} ({}, {}, {}), Ez_max = {} ({}, {}, {}).\n\n",
        Ez_min, Ez_min_pqr.p, Ez_min_pqr.q, Ez_min_pqr.r,
        Ez_max, Ez_max_pqr.p, Ez_max_pqr.q, Ez_max_pqr.r);

    // Test views...

    std::print("Ex_extent(0) = {}, p_ - 1 = {}, ", Ex_.extent(0), p_ - 1);
    std::print("Ex_extent(1) = {}, q_ - 1 = {}, ", Ex_.extent(1), q_ - 1);
    std::print("Ex_extent(2) = {}, r_ - 1 = {}.", Ex_.extent(2), r_ - 1);

    auto d0 = std::views::iota((size_t) 0, Ex_.extent(0));
    auto d1 = std::views::iota((size_t) 0, Ex_.extent(1));
    auto d2 = std::views::iota((size_t) 0, Ex_.extent(2));

    for (auto [i0, i1, i2] : std::views::cartesian_product(d0, d1, d2)) {
        Ex_[i0, i1, i2] = 0.;
        ExEx_[i0, i1, i2] = 0.;
        ExHyHz_[i0, i1, i2] = 0.;
        std::print("1 Ex_[{}, {}, {}] = {}, ", i0, i1, i2, Ex_[i0, i1, i2]);
        std::print("ExEx_[{}, {}, {}] = {}, ", i0, i1, i2, ExEx_[i0, i1, i2]);
        std::print("ExHyHz_[{}, {}, {}] = {}.\n", i0, i1, i2, ExHyHz_[i0, i1, i2]);
    }

    std::print("\n");

    auto gridEx = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1),
        std::views::iota((size_t) 1, q_ - 1),
        std::views::iota((size_t) 1, r_ - 1));

    std::for_each(gridEx.begin(), gridEx.end(), [&](auto tuple) {
        auto [i0, i1, i2] = tuple;
        Ex_[i0, i1, i2] = 0.;
        ExEx_[i0, i1, i2] = 0.;
        ExHyHz_[i0, i1, i2] = 0.;
        std::print("2 Ex_[{}, {}, {}] = {}, ", i0, i1, i2, Ex_[i0, i1, i2]);
        std::print("ExEx_[{}, {}, {}] = {}, ", i0, i1, i2, ExEx_[i0, i1, i2]);
        std::print("ExHyHz_[{}, {}, {}] = {}.\n", i0, i1, i2, ExHyHz_[i0, i1, i2]);
    });

    std::print("\n");

    #endif

}

void
FDTDCPML::H_source()
{
}

void
FDTDCPML::E_source()
{

    double x0 = std::lerp(x_min_, x_max_, 0.5);
    double y0 = std::lerp(y_min_, y_max_, 0.5);
    double z0 = std::lerp(z_min_, z_max_, 0.2);

    double f = 1.0e9;  // 1 GGz

    int l = 2.0;
    int n = 2.0;
    double w = std::min(x_max_ - x_min_, y_max_ - y_min_) / 8.0;

    E_source_LG(t_, x0, y0, z0, f, n, l, w);

}

void
FDTDCPML::E_source_LG(double t, double x0, double y0, double z0, double freq, int n, int l, double w)
{

//    std::cout << "x0 = " << x0 << ", y0 = " << y0 << ", z0 = " << z0 << std::endl;
//    std::cout << "freq = " << freq << ", n = " << n << ", l = " << l << ", w = " << w << std::endl;

    const double omega = 2.0 * std::numbers::pi * freq;
    const double period = 1.0 / freq;
    const double smooth = 1.0 - std::exp(- std::pow(t / (3.0 * period), 2));

    auto grid = std::views::cartesian_product(
        std::views::iota(p_abs_l_, p_abs_l_ + p_int_), 
        std::views::iota(q_abs_l_, q_abs_l_ + q_int_)
    );

    const auto [p0, q0, r0] = Ex_xyz2pqr({x0, y0, z0});

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {

        auto [p, q] = tuple;

        auto xyz = Ex_pqr2xyz({p, q, r0});

        double dx = xyz.x - x0;
        double dy = xyz.y - y0;
        double rho = std::hypot(dx, dy);
        double phi = std::atan2(dy, dx);
        
        double poly = std::assoc_laguerre(n, std::abs(l), 2.0 * std::pow(rho / w, 2));
        
        double amp = smooth
            * std::exp(- std::pow(rho / w, 2)) 
            * std::pow(std::sqrt(2.0) * rho / w, std::abs(l)) 
            * poly;

        double coeff = t_step_ / (fdtdcpml::eps0 * x_step_);  // Because we use Ex_ / x_step_.

        Ex_[p, q, r0] += amp * std::cos(omega * t - l * phi) * coeff;

    });

    // const auto [p0, q0, r0] = Ey_xyz2pqr({x0, y0, z0});  // Should give the same r0 value as Ex...

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {

        auto [p, q] = tuple;

        auto xyz = Ey_pqr2xyz({p, q, r0});

        double dx = xyz.x - x0;
        double dy = xyz.y - y0;
        double rho = std::hypot(dx, dy);
        double phi = std::atan2(dy, dx);
        
        double poly = std::assoc_laguerre(n, std::abs(l), 2.0 * std::pow(rho / w, 2));

        double amp = smooth
            * std::exp(- std::pow(rho / w, 2)) 
            * std::pow(std::sqrt(2.0) * rho / w, std::abs(l)) 
            * poly;

        double coeff = t_step_ / (fdtdcpml::eps0 * y_step_);  // Because we use Ey_ / y_step_.

        Ey_[p, q, r0] += amp * std::sin(omega * t - l * phi) * coeff;

    });

}
