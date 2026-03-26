#include "fdtdcpml.h"

#include <iostream>
#include <cmath>
#include <ranges>
#include <algorithm>
#include <execution>

// "Premature optimization is the root of all evil". Donald Knuth.

void FDTDCPML::cpml_x_l_H() {

    auto& cpml = psi_x_l_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Hy_x = cpml.psi_H1_;
    auto& psi_Hz_x = cpml.psi_H2_;

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_l_), 
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_ - 1));

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = 1 + p; 

        double dEz_dx = Ez_[gp + 1, q, r] - Ez_[gp, q, r];
        psi_Hy_x[q, r, p] = b[p] * psi_Hy_x[q, r, p] + c[p] * dEz_dx;
        Hy_[gp, q, r] += HyEzEx_[gp, q, r] * (dEz_dx * (k_inv[p] - 1.0) + psi_Hy_x[q, r, p]);

        double dEy_dx = Ey_[gp + 1, q, r] - Ey_[gp, q, r];
        psi_Hz_x[q, r, p] = b[p] * psi_Hz_x[q, r, p] + c[p] * dEy_dx;
        Hz_[gp, q, r] -= HzExEy_[gp, q, r] * (dEy_dx * (k_inv[p] - 1.0) + psi_Hz_x[q, r, p]);
    });

}

void FDTDCPML::cpml_x_h_H() {

    auto& cpml = psi_x_h_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Hy_x = cpml.psi_H1_;
    auto& psi_Hz_x = cpml.psi_H2_;

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_h_), 
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_ - 1));

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = p_ - 2 - p; 

        double dEz_dx = Ez_[gp + 1, q, r] - Ez_[gp, q, r];
        psi_Hy_x[q, r, p] = b[p] * psi_Hy_x[q, r, p] + c[p] * dEz_dx;
        Hy_[gp, q, r] += HyEzEx_[gp, q, r] * (dEz_dx * (k_inv[p] - 1.0) + psi_Hy_x[q, r, p]);

        double dEy_dx = Ey_[gp + 1, q, r] - Ey_[gp, q, r];
        psi_Hz_x[q, r, p] = b[p] * psi_Hz_x[q, r, p] + c[p] * dEy_dx;
        Hz_[gp, q, r] -= HzExEy_[gp, q, r] * (dEy_dx * (k_inv[p] - 1.0) + psi_Hz_x[q, r, p]);
    });

}

void FDTDCPML::cpml_y_l_H() {

    auto& cpml = psi_y_l_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Hx_y = cpml.psi_H1_;
    auto& psi_Hz_y = cpml.psi_H2_;

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t)0, p_ - 1), 
        std::views::iota((size_t)0, q_abs_l_),
        std::views::iota((size_t)0, r_ - 1));

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = 1 + q; 

        double dEz_dy = Ez_[p, gq + 1, r] - Ez_[p, gq, r];
        psi_Hx_y[p, r, q] = b[q] * psi_Hx_y[p, r, q] + c[q] * dEz_dy;
        Hx_[p, gq, r] -= HxEyEz_[p, gq, r] * (dEz_dy * (k_inv[q] - 1.0) + psi_Hx_y[p, r, q]);

        double dEx_dy = Ex_[p, gq + 1, r] - Ex_[p, gq, r];
        psi_Hz_y[p, r, q] = b[q] * psi_Hz_y[p, r, q] + c[q] * dEx_dy;
        Hz_[p, gq, r] += HzExEy_[p, gq, r] * (dEx_dy * (k_inv[q] - 1.0) + psi_Hz_y[p, r, q]);
    });

}

void FDTDCPML::cpml_y_h_H() {

    auto& cpml = psi_y_h_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Hx_y = cpml.psi_H1_;
    auto& psi_Hz_y = cpml.psi_H2_;

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1), 
        std::views::iota((size_t) 0, q_abs_h_),
        std::views::iota((size_t) 0, r_ - 1));

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = q_ - 2 - q; 

        double dEz_dy = Ez_[p, gq + 1, r] - Ez_[p, gq, r];
        psi_Hx_y[p, r, q] = b[q] * psi_Hx_y[p, r, q] + c[q] * dEz_dy;
        Hx_[p, gq, r] -= HxEyEz_[p, gq, r] * (dEz_dy * (k_inv[q] - 1.0) + psi_Hx_y[p, r, q]);

        double dEx_dy = Ex_[p, gq + 1, r] - Ex_[p, gq, r];
        psi_Hz_y[p, r, q] = b[q] * psi_Hz_y[p, r, q] + c[q] * dEx_dy;
        Hz_[p, gq, r] += HzExEy_[p, gq, r] * (dEx_dy * (k_inv[q] - 1.0) + psi_Hz_y[p, r, q]);
    });

}

void FDTDCPML::cpml_z_l_H() {

    auto& cpml = psi_z_l_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Hx_z = cpml.psi_H1_;
    auto& psi_Hy_z = cpml.psi_H2_;

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t)0, p_ - 1), 
        std::views::iota((size_t)0, q_ - 1),
        std::views::iota((size_t)0, r_abs_l_));

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = 1 + r; 

        double dEy_dz = Ey_[p, q, gr + 1] - Ey_[p, q, gr];
        psi_Hx_z[p, q, r] = b[r] * psi_Hx_z[p, q, r] + c[r] * dEy_dz;
        Hx_[p, q, gr] += HxEyEz_[p, q, gr] * (dEy_dz * (k_inv[r] - 1.0) + psi_Hx_z[p, q, r]);

        double dEx_dz = Ex_[p, q, gr + 1] - Ex_[p, q, gr];
        psi_Hy_z[p, q, r] = b[r] * psi_Hy_z[p, q, r] + c[r] * dEx_dz;
        Hy_[p, q, gr] -= HyEzEx_[p, q, gr] * (dEx_dz * (k_inv[r] - 1.0) + psi_Hy_z[p, q, r]);
    });

}

void FDTDCPML::cpml_z_h_H() {

    auto& cpml = psi_z_h_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Hx_z = cpml.psi_H1_;
    auto& psi_Hy_z = cpml.psi_H2_;

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1), 
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_abs_h_));

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = r_ - 2 - r; 

        double dEy_dz = Ey_[p, q, gr + 1] - Ey_[p, q, gr];
        psi_Hx_z[p, q, r] = b[r] * psi_Hx_z[p, q, r] + c[r] * dEy_dz;
        Hx_[p, q, gr] += HxEyEz_[p, q, gr] * (dEy_dz * (k_inv[r] - 1.0) + psi_Hx_z[p, q, r]);

        double dEx_dz = Ex_[p, q, gr + 1] - Ex_[p, q, gr];
        psi_Hy_z[p, q, r] = b[r] * psi_Hy_z[p, q, r] + c[r] * dEx_dz;
        Hy_[p, q, gr] -= HyEzEx_[p, q, gr] * (dEx_dz * (k_inv[r] - 1.0) + psi_Hy_z[p, q, r]);
    });

}

void FDTDCPML::cpml_x_l_E() {

    auto& cpml = psi_x_l_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Ey_x = cpml.psi_E1_; // Влияет на Ey через dHz/dx
    auto& psi_Ez_x = cpml.psi_E2_; // Влияет на Ez через dHy/dx

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_l_), 
        std::views::iota((size_t) 1, q_ - 1),
        std::views::iota((size_t) 1, r_ - 1)
    );

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = 1 + p; 

        double dHz_dx = Hz_[gp, q, r] - Hz_[gp - 1, q, r];
        psi_Ey_x[q, r, p] = b[p] * psi_Ey_x[q, r, p] + c[p] * dHz_dx;
        Ey_[gp, q, r] -= EyHzHx_[gp, q, r] * (dHz_dx * (k_inv[p] - 1.0) + psi_Ey_x[q, r, p]);

        double dHy_dx = Hy_[gp, q, r] - Hy_[gp - 1, q, r];
        psi_Ez_x[q, r, p] = b[p] * psi_Ez_x[q, r, p] + c[p] * dHy_dx;
        Ez_[gp, q, r] += EzHxHy_[gp, q, r] * (dHy_dx * (k_inv[p] - 1.0) + psi_Ez_x[q, r, p]);
    });

}

void FDTDCPML::cpml_x_h_E() {

    auto& cpml = psi_x_h_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Ey_x = cpml.psi_E1_;  // Ey = ... - dHz/dx
    auto& psi_Ez_x = cpml.psi_E2_;  // Ez = ... + dHy/dx

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_h_), 
        std::views::iota((size_t) 1, q_ - 1),
        std::views::iota((size_t) 1, r_ - 1)
    );

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = p_ - 2 - p; 

        double dHz_dx = Hz_[gp, q, r] - Hz_[gp - 1, q, r];
        psi_Ey_x[q, r, p] = b[p] * psi_Ey_x[q, r, p] + c[p] * dHz_dx;
        Ey_[gp, q, r] -= EyHzHx_[gp, q, r] * (dHz_dx * (k_inv[p] - 1.0) + psi_Ey_x[q, r, p]);

        double dHy_dx = Hy_[gp, q, r] - Hy_[gp - 1, q, r];
        psi_Ez_x[q, r, p] = b[p] * psi_Ez_x[q, r, p] + c[p] * dHy_dx;
        Ez_[gp, q, r] += EzHxHy_[gp, q, r] * (dHy_dx * (k_inv[p] - 1.0) + psi_Ez_x[q, r, p]);
    });

}

void FDTDCPML::cpml_y_l_E() {

    auto& cpml = psi_y_l_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Ex_y = cpml.psi_E1_;  // psi_E1_ отвечает за dHz/dy (влияет на Ex)
    auto& psi_Ez_y = cpml.psi_E2_;  // psi_E2_ отвечает за dHx/dy (влияет на Ez)

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1), 
        std::views::iota((size_t) 0, q_abs_l_),
        std::views::iota((size_t) 1, r_ - 1)
    );

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = 1 + q;  // global mesh q index (пропускаем PEC на индексе 0)

        double dHz_dy = Hz_[p, gq, r] - Hz_[p, gq - 1, r];
        psi_Ex_y[p, r, q] = b[q] * psi_Ex_y[p, r, q] + c[q] * dHz_dy;
        Ex_[p, gq, r] += ExHyHz_[p, gq, r] * (dHz_dy * (k_inv[q] - 1.0) + psi_Ex_y[p, r, q]);

        double dHx_dy = Hx_[p, gq, r] - Hx_[p, gq - 1, r];
        psi_Ez_y[p, r, q] = b[q] * psi_Ez_y[p, r, q] + c[q] * dHx_dy;
        Ez_[p, gq, r] -= EzHxHy_[p, gq, r] * (dHx_dy * (k_inv[q] - 1.0) + psi_Ez_y[p, r, q]);
    });

}

void FDTDCPML::cpml_y_h_E() {

    auto& cpml = psi_y_h_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Ex_y = cpml.psi_E1_;  // psi_E1_ для Y-грани — это dHz/dy (влияет на Ex)
    auto& psi_Ez_y = cpml.psi_E2_;  // psi_E2_ для Y-грани — это dHx/dy (влияет на Ez)

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1), 
        std::views::iota((size_t) 0, q_abs_h_),
        std::views::iota((size_t) 1, r_ - 1)
    );

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = q_ - 2 - q; 

        double dHz_dy = Hz_[p, gq, r] - Hz_[p, gq - 1, r];
        psi_Ex_y[p, r, q] = b[q] * psi_Ex_y[p, r, q] + c[q] * dHz_dy;
        Ex_[p, gq, r] += ExHyHz_[p, gq, r] * (dHz_dy * (k_inv[q] - 1.0) + psi_Ex_y[p, r, q]);

        double dHx_dy = Hx_[p, gq, r] - Hx_[p, gq - 1, r];
        psi_Ez_y[p, r, q] = b[q] * psi_Ez_y[p, r, q] + c[q] * dHx_dy;
        Ez_[p, gq, r] -= EzHxHy_[p, gq, r] * (dHx_dy * (k_inv[q] - 1.0) + psi_Ez_y[p, r, q]);
    });

}

void FDTDCPML::cpml_z_l_E() {

    auto& cpml = psi_z_l_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Ex_z = cpml.psi_E1_;  // psi_E1_ отвечает за dHy/dz (влияет на Ex)
    auto& psi_Ey_z = cpml.psi_E2_;  // psi_E2_ отвечает за dHx/dz (влияет на Ey)

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1), 
        std::views::iota((size_t) 1, q_ - 1),
        std::views::iota((size_t) 0, r_abs_l_)
    );

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = 1 + r;  // global mesh r index (пропускаем PEC на индексе 0)

        double dHy_dz = Hy_[p, q, gr] - Hy_[p, q, gr - 1];
        psi_Ex_z[p, q, r] = b[r] * psi_Ex_z[p, q, r] + c[r] * dHy_dz;
        Ex_[p, q, gr] -= ExHyHz_[p, q, gr] * (dHy_dz * (k_inv[r] - 1.0) + psi_Ex_z[p, q, r]);

        double dHx_dz = Hx_[p, q, gr] - Hx_[p, q, gr - 1];
        psi_Ey_z[p, q, r] = b[r] * psi_Ey_z[p, q, r] + c[r] * dHx_dz;
        Ey_[p, q, gr] += EyHzHx_[p, q, gr] * (dHx_dz * (k_inv[r] - 1.0) + psi_Ey_z[p, q, r]);
    });

}

void FDTDCPML::cpml_z_h_E() {

    auto& cpml = psi_z_h_;

    auto& b = cpml.b_;
    auto& c = cpml.c_;
    auto& k_inv = cpml.k_inv_;

    auto& psi_Ex_z = cpml.psi_E1_;  // psi_E1_ отвечает за dHy/dz (влияет на Ex)
    auto& psi_Ey_z = cpml.psi_E2_;  // psi_E2_ отвечает за dHx/dz (влияет на Ey)

    auto grid = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1), 
        std::views::iota((size_t) 1, q_ - 1),
        std::views::iota((size_t) 0, r_abs_h_)
    );

    std::for_each(std::execution::par, grid.begin(), grid.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = r_ - 2 - r; 

        double dHy_dz = Hy_[p, q, gr] - Hy_[p, q, gr - 1];
        psi_Ex_z[p, q, r] = b[r] * psi_Ex_z[p, q, r] + c[r] * dHy_dz;
        Ex_[p, q, gr] -= ExHyHz_[p, q, gr] * (dHy_dz * (k_inv[r] - 1.0) + psi_Ex_z[p, q, r]);

        double dHx_dz = Hx_[p, q, gr] - Hx_[p, q, gr - 1];
        psi_Ey_z[p, q, r] = b[r] * psi_Ey_z[p, q, r] + c[r] * dHx_dz;
        Ey_[p, q, gr] += EyHzHx_[p, q, gr] * (dHx_dz * (k_inv[r] - 1.0) + psi_Ey_z[p, q, r]);
    });

}
