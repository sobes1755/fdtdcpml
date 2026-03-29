#include "fdtdcpml.h"

#include <iostream>
#include <cmath>
#include <ranges>
#include <algorithm>
#include <execution>

// "Premature optimization is the root of all evil". Donald Knuth.

void FDTDCPML::cpml_x_l_H() {

    auto& cpml = psi_x_l_;

    auto& b = cpml.b_H_;
    auto& c = cpml.c_H_;
    auto& k_inv = cpml.k_inv_H_;

    auto& psi_Hy = cpml.psi_H1_;
    auto& psi_Hz = cpml.psi_H2_;

    auto grid_Hy = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_l_), 
        std::views::iota((size_t) 0, q_),
        std::views::iota((size_t) 0, r_ - 1));
    auto grid_Hz = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_l_), 
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_));

    std::for_each(std::execution::par, grid_Hy.begin(), grid_Hy.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = 1 + p; 

        double dEz_dx = Ez_[gp + 1, q, r] - Ez_[gp, q, r];
        psi_Hy[q, r, p] = b[p] * psi_Hy[q, r, p] + c[p] * dEz_dx;
        Hy_[gp, q, r] += HyEzEx_[gp, q, r] * (dEz_dx * (k_inv[p] - 1.0) + psi_Hy[q, r, p]);
    });
    std::for_each(std::execution::par, grid_Hz.begin(), grid_Hz.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = 1 + p; 

        double dEy_dx = Ey_[gp + 1, q, r] - Ey_[gp, q, r];
        psi_Hz[q, r, p] = b[p] * psi_Hz[q, r, p] + c[p] * dEy_dx;
        Hz_[gp, q, r] -= HzExEy_[gp, q, r] * (dEy_dx * (k_inv[p] - 1.0) + psi_Hz[q, r, p]);
    });

}

void FDTDCPML::cpml_x_h_H() {

    auto& cpml = psi_x_h_;

    auto& b = cpml.b_H_;
    auto& c = cpml.c_H_;
    auto& k_inv = cpml.k_inv_H_;

    auto& psi_Hy = cpml.psi_H1_;
    auto& psi_Hz = cpml.psi_H2_;

    auto grid_Hy = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_h_), 
        std::views::iota((size_t) 0, q_),
        std::views::iota((size_t) 0, r_ - 1));
    auto grid_Hz = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_h_), 
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_));

    std::for_each(std::execution::par, grid_Hy.begin(), grid_Hy.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = p_ - 2 - p; 

        double dEz_dx = Ez_[gp + 1, q, r] - Ez_[gp, q, r];
        psi_Hy[q, r, p] = b[p] * psi_Hy[q, r, p] + c[p] * dEz_dx;
        Hy_[gp, q, r] += HyEzEx_[gp, q, r] * (dEz_dx * (k_inv[p] - 1.0) + psi_Hy[q, r, p]);
    });
    std::for_each(std::execution::par, grid_Hz.begin(), grid_Hz.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = p_ - 2 - p; 

        double dEy_dx = Ey_[gp + 1, q, r] - Ey_[gp, q, r];
        psi_Hz[q, r, p] = b[p] * psi_Hz[q, r, p] + c[p] * dEy_dx;
        Hz_[gp, q, r] -= HzExEy_[gp, q, r] * (dEy_dx * (k_inv[p] - 1.0) + psi_Hz[q, r, p]);
    });

}

void FDTDCPML::cpml_y_l_H() {

    auto& cpml = psi_y_l_;

    auto& b = cpml.b_H_;
    auto& c = cpml.c_H_;
    auto& k_inv = cpml.k_inv_H_;

    auto& psi_Hx = cpml.psi_H1_;
    auto& psi_Hz = cpml.psi_H2_;

    auto grid_Hx = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_), 
        std::views::iota((size_t) 0, q_abs_l_),
        std::views::iota((size_t) 0, r_ - 1));
    auto grid_Hz = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1), 
        std::views::iota((size_t) 0, q_abs_l_),
        std::views::iota((size_t) 0, r_));

    std::for_each(std::execution::par, grid_Hx.begin(), grid_Hx.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = 1 + q; 

        double dEz_dy = Ez_[p, gq + 1, r] - Ez_[p, gq, r];
        psi_Hx[p, r, q] = b[q] * psi_Hx[p, r, q] + c[q] * dEz_dy;
        Hx_[p, gq, r] -= HxEyEz_[p, gq, r] * (dEz_dy * (k_inv[q] - 1.0) + psi_Hx[p, r, q]);
    });
    std::for_each(std::execution::par, grid_Hz.begin(), grid_Hz.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = 1 + q; 

        double dEx_dy = Ex_[p, gq + 1, r] - Ex_[p, gq, r];
        psi_Hz[p, r, q] = b[q] * psi_Hz[p, r, q] + c[q] * dEx_dy;
        Hz_[p, gq, r] += HzExEy_[p, gq, r] * (dEx_dy * (k_inv[q] - 1.0) + psi_Hz[p, r, q]);
    });

}

void FDTDCPML::cpml_y_h_H() {

    auto& cpml = psi_y_h_;

    auto& b = cpml.b_H_;
    auto& c = cpml.c_H_;
    auto& k_inv = cpml.k_inv_H_;

    auto& psi_Hx = cpml.psi_H1_;
    auto& psi_Hz = cpml.psi_H2_;

    auto grid_Hx = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_), 
        std::views::iota((size_t) 0, q_abs_h_),
        std::views::iota((size_t) 0, r_ - 1));
    auto grid_Hz = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1), 
        std::views::iota((size_t) 0, q_abs_h_),
        std::views::iota((size_t) 0, r_));

    std::for_each(std::execution::par, grid_Hx.begin(), grid_Hx.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = q_ - 2 - q; 

        double dEz_dy = Ez_[p, gq + 1, r] - Ez_[p, gq, r];
        psi_Hx[p, r, q] = b[q] * psi_Hx[p, r, q] + c[q] * dEz_dy;
        Hx_[p, gq, r] -= HxEyEz_[p, gq, r] * (dEz_dy * (k_inv[q] - 1.0) + psi_Hx[p, r, q]);
    });
    std::for_each(std::execution::par, grid_Hz.begin(), grid_Hz.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = q_ - 2 - q; 

        double dEx_dy = Ex_[p, gq + 1, r] - Ex_[p, gq, r];
        psi_Hz[p, r, q] = b[q] * psi_Hz[p, r, q] + c[q] * dEx_dy;
        Hz_[p, gq, r] += HzExEy_[p, gq, r] * (dEx_dy * (k_inv[q] - 1.0) + psi_Hz[p, r, q]);
    });

}

void FDTDCPML::cpml_z_l_H() {

    auto& cpml = psi_z_l_;

    auto& b = cpml.b_H_;
    auto& c = cpml.c_H_;
    auto& k_inv = cpml.k_inv_H_;

    auto& psi_Hx = cpml.psi_H1_;
    auto& psi_Hy = cpml.psi_H2_;

    auto grid_Hx = std::views::cartesian_product(
        std::views::iota((size_t)0, p_), 
        std::views::iota((size_t)0, q_ - 1),
        std::views::iota((size_t)0, r_abs_l_));
    auto grid_Hy = std::views::cartesian_product(
        std::views::iota((size_t)0, p_ - 1), 
        std::views::iota((size_t)0, q_),
        std::views::iota((size_t)0, r_abs_l_));

    std::for_each(std::execution::par, grid_Hx.begin(), grid_Hx.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = 1 + r; 

        double dEy_dz = Ey_[p, q, gr + 1] - Ey_[p, q, gr];
        psi_Hx[p, q, r] = b[r] * psi_Hx[p, q, r] + c[r] * dEy_dz;
        Hx_[p, q, gr] += HxEyEz_[p, q, gr] * (dEy_dz * (k_inv[r] - 1.0) + psi_Hx[p, q, r]);
    });
    std::for_each(std::execution::par, grid_Hy.begin(), grid_Hy.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = 1 + r; 

        double dEx_dz = Ex_[p, q, gr + 1] - Ex_[p, q, gr];
        psi_Hy[p, q, r] = b[r] * psi_Hy[p, q, r] + c[r] * dEx_dz;
        Hy_[p, q, gr] -= HyEzEx_[p, q, gr] * (dEx_dz * (k_inv[r] - 1.0) + psi_Hy[p, q, r]);
    });

}

void FDTDCPML::cpml_z_h_H() {

    auto& cpml = psi_z_h_;

    auto& b = cpml.b_H_;
    auto& c = cpml.c_H_;
    auto& k_inv = cpml.k_inv_H_;

    auto& psi_Hx = cpml.psi_H1_;
    auto& psi_Hy = cpml.psi_H2_;

    auto grid_Hx = std::views::cartesian_product(
        std::views::iota((size_t)0, p_), 
        std::views::iota((size_t)0, q_ - 1),
        std::views::iota((size_t)0, r_abs_h_));
    auto grid_Hy = std::views::cartesian_product(
        std::views::iota((size_t)0, p_ - 1), 
        std::views::iota((size_t)0, q_),
        std::views::iota((size_t)0, r_abs_h_));

    std::for_each(std::execution::par, grid_Hx.begin(), grid_Hx.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = r_ - 2 - r; 

        double dEy_dz = Ey_[p, q, gr + 1] - Ey_[p, q, gr];
        psi_Hx[p, q, r] = b[r] * psi_Hx[p, q, r] + c[r] * dEy_dz;
        Hx_[p, q, gr] += HxEyEz_[p, q, gr] * (dEy_dz * (k_inv[r] - 1.0) + psi_Hx[p, q, r]);
    });
    std::for_each(std::execution::par, grid_Hy.begin(), grid_Hy.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = r_ - 2 - r; 

        double dEx_dz = Ex_[p, q, gr + 1] - Ex_[p, q, gr];
        psi_Hy[p, q, r] = b[r] * psi_Hy[p, q, r] + c[r] * dEx_dz;
        Hy_[p, q, gr] -= HyEzEx_[p, q, gr] * (dEx_dz * (k_inv[r] - 1.0) + psi_Hy[p, q, r]);
    });

}

void FDTDCPML::cpml_x_l_E() {

    auto& cpml = psi_x_l_;

    auto& b = cpml.b_E_;
    auto& c = cpml.c_E_;
    auto& k_inv = cpml.k_inv_E_;

    auto& psi_Ey = cpml.psi_E1_;  // Влияет на Ey через dHz/dx
    auto& psi_Ez = cpml.psi_E2_;  // Влияет на Ez через dHy/dx

    auto grid_Ey = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_l_), 
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 1, r_ - 1)  // PEC
    );
    auto grid_Ez = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_l_), 
        std::views::iota((size_t) 1, q_ - 1),  // PEC
        std::views::iota((size_t) 0, r_ - 1) 
    );

    std::for_each(std::execution::par, grid_Ey.begin(), grid_Ey.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = 1 + p; 

        double dHz_dx = Hz_[gp, q, r] - Hz_[gp - 1, q, r];
        psi_Ey[q, r, p] = b[p] * psi_Ey[q, r, p] + c[p] * dHz_dx;
        Ey_[gp, q, r] -= EyHzHx_[gp, q, r] * (dHz_dx * (k_inv[p] - 1.0) + psi_Ey[q, r, p]);
    });
    std::for_each(std::execution::par, grid_Ez.begin(), grid_Ez.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = 1 + p; 

        double dHy_dx = Hy_[gp, q, r] - Hy_[gp - 1, q, r];
        psi_Ez[q, r, p] = b[p] * psi_Ez[q, r, p] + c[p] * dHy_dx;
        Ez_[gp, q, r] += EzHxHy_[gp, q, r] * (dHy_dx * (k_inv[p] - 1.0) + psi_Ez[q, r, p]);
    });

}

void FDTDCPML::cpml_x_h_E() {

    auto& cpml = psi_x_h_;

    auto& b = cpml.b_E_;
    auto& c = cpml.c_E_;
    auto& k_inv = cpml.k_inv_E_;

    auto& psi_Ey = cpml.psi_E1_;  // Влияет на Ey через dHz/dx
    auto& psi_Ez = cpml.psi_E2_;  // Влияет на Ez через dHy/dx

    auto grid_Ey = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_h_), 
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 1, r_ - 1)  // PEC
    );
    auto grid_Ez = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_abs_h_), 
        std::views::iota((size_t) 1, q_ - 1),  // PEC
        std::views::iota((size_t) 0, r_ - 1) 
    );

    std::for_each(std::execution::par, grid_Ey.begin(), grid_Ey.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = p_ - 2 - p; 

        double dHz_dx = Hz_[gp, q, r] - Hz_[gp - 1, q, r];
        psi_Ey[q, r, p] = b[p] * psi_Ey[q, r, p] + c[p] * dHz_dx;
        Ey_[gp, q, r] -= EyHzHx_[gp, q, r] * (dHz_dx * (k_inv[p] - 1.0) + psi_Ey[q, r, p]);
    });
    std::for_each(std::execution::par, grid_Ez.begin(), grid_Ez.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gp = p_ - 2 - p; 

        double dHy_dx = Hy_[gp, q, r] - Hy_[gp - 1, q, r];
        psi_Ez[q, r, p] = b[p] * psi_Ez[q, r, p] + c[p] * dHy_dx;
        Ez_[gp, q, r] += EzHxHy_[gp, q, r] * (dHy_dx * (k_inv[p] - 1.0) + psi_Ez[q, r, p]);
    });

}

void FDTDCPML::cpml_y_l_E() {

    auto& cpml = psi_y_l_;

    auto& b = cpml.b_E_;
    auto& c = cpml.c_E_;
    auto& k_inv = cpml.k_inv_E_;

    auto& psi_Ex = cpml.psi_E1_;  // Влияет на Ex через dHz/dy
    auto& psi_Ez = cpml.psi_E2_;  // Влияет на Ez через dHx/dy

    auto grid_Ex = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1), 
        std::views::iota((size_t) 0, q_abs_l_),
        std::views::iota((size_t) 1, r_ - 1)  // PEC
    );
    auto grid_Ez = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1),  // PEC
        std::views::iota((size_t) 0, q_abs_l_),
        std::views::iota((size_t) 0, r_ - 1) 
    );

    std::for_each(std::execution::par, grid_Ex.begin(), grid_Ex.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = 1 + q;

        double dHz_dy = Hz_[p, gq, r] - Hz_[p, gq - 1, r];
        psi_Ex[p, r, q] = b[q] * psi_Ex[p, r, q] + c[q] * dHz_dy;
        Ex_[p, gq, r] += ExHyHz_[p, gq, r] * (dHz_dy * (k_inv[q] - 1.0) + psi_Ex[p, r, q]);
    });
    std::for_each(std::execution::par, grid_Ez.begin(), grid_Ez.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = 1 + q;

        double dHx_dy = Hx_[p, gq, r] - Hx_[p, gq - 1, r];
        psi_Ez[p, r, q] = b[q] * psi_Ez[p, r, q] + c[q] * dHx_dy;
        Ez_[p, gq, r] -= EzHxHy_[p, gq, r] * (dHx_dy * (k_inv[q] - 1.0) + psi_Ez[p, r, q]);
    });

}

void FDTDCPML::cpml_y_h_E() {

    auto& cpml = psi_y_h_;

    auto& b = cpml.b_E_;
    auto& c = cpml.c_E_;
    auto& k_inv = cpml.k_inv_E_;

    auto& psi_Ex = cpml.psi_E1_;  // Влияет на Ex через dHz/dy
    auto& psi_Ez = cpml.psi_E2_;  // Влияет на Ez через dHx/dy

    auto grid_Ex = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1), 
        std::views::iota((size_t) 0, q_abs_h_),
        std::views::iota((size_t) 1, r_ - 1)  // PEC
    );
    auto grid_Ez = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1),  // PEC
        std::views::iota((size_t) 0, q_abs_h_),
        std::views::iota((size_t) 0, r_ - 1) 
    );

    std::for_each(std::execution::par, grid_Ex.begin(), grid_Ex.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = q_ - 2 - q; 

        double dHz_dy = Hz_[p, gq, r] - Hz_[p, gq - 1, r];
        psi_Ex[p, r, q] = b[q] * psi_Ex[p, r, q] + c[q] * dHz_dy;
        Ex_[p, gq, r] += ExHyHz_[p, gq, r] * (dHz_dy * (k_inv[q] - 1.0) + psi_Ex[p, r, q]);
    });
    std::for_each(std::execution::par, grid_Ez.begin(), grid_Ez.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gq = q_ - 2 - q; 

        double dHx_dy = Hx_[p, gq, r] - Hx_[p, gq - 1, r];
        psi_Ez[p, r, q] = b[q] * psi_Ez[p, r, q] + c[q] * dHx_dy;
        Ez_[p, gq, r] -= EzHxHy_[p, gq, r] * (dHx_dy * (k_inv[q] - 1.0) + psi_Ez[p, r, q]);
    });

}

void FDTDCPML::cpml_z_l_E() {

    auto& cpml = psi_z_l_;

    auto& b = cpml.b_E_;
    auto& c = cpml.c_E_;
    auto& k_inv = cpml.k_inv_E_;

    auto& psi_Ex = cpml.psi_E1_;  // Влияет на Ex через dHy/dz
    auto& psi_Ey = cpml.psi_E2_;  // Влияет на Ey через dHx/dz

    auto grid_Ex = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1), 
        std::views::iota((size_t) 1, q_ - 1),  // PEC
        std::views::iota((size_t) 0, r_abs_l_)
    );
    auto grid_Ey = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1),  // PEC
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_abs_l_) 
    );

    std::for_each(std::execution::par, grid_Ex.begin(), grid_Ex.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = 1 + r;

        double dHy_dz = Hy_[p, q, gr] - Hy_[p, q, gr - 1];
        psi_Ex[p, q, r] = b[r] * psi_Ex[p, q, r] + c[r] * dHy_dz;
        Ex_[p, q, gr] -= ExHyHz_[p, q, gr] * (dHy_dz * (k_inv[r] - 1.0) + psi_Ex[p, q, r]);
    });
    std::for_each(std::execution::par, grid_Ey.begin(), grid_Ey.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = 1 + r;

        double dHx_dz = Hx_[p, q, gr] - Hx_[p, q, gr - 1];
        psi_Ey[p, q, r] = b[r] * psi_Ey[p, q, r] + c[r] * dHx_dz;
        Ey_[p, q, gr] += EyHzHx_[p, q, gr] * (dHx_dz * (k_inv[r] - 1.0) + psi_Ey[p, q, r]);
    });

}

void FDTDCPML::cpml_z_h_E() {

    auto& cpml = psi_z_h_;

    auto& b = cpml.b_E_;
    auto& c = cpml.c_E_;
    auto& k_inv = cpml.k_inv_E_;

    auto& psi_Ex = cpml.psi_E1_;  // Влияет на Ex через dHy/dz
    auto& psi_Ey = cpml.psi_E2_;  // Влияет на Ey через dHx/dz

    auto grid_Ex = std::views::cartesian_product(
        std::views::iota((size_t) 0, p_ - 1), 
        std::views::iota((size_t) 1, q_ - 1),  // PEC
        std::views::iota((size_t) 0, r_abs_h_)
    );
    auto grid_Ey = std::views::cartesian_product(
        std::views::iota((size_t) 1, p_ - 1),  // PEC
        std::views::iota((size_t) 0, q_ - 1),
        std::views::iota((size_t) 0, r_abs_h_) 
    );

    std::for_each(std::execution::par, grid_Ex.begin(), grid_Ex.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = r_ - 2 - r; 

        double dHy_dz = Hy_[p, q, gr] - Hy_[p, q, gr - 1];
        psi_Ex[p, q, r] = b[r] * psi_Ex[p, q, r] + c[r] * dHy_dz;
        Ex_[p, q, gr] -= ExHyHz_[p, q, gr] * (dHy_dz * (k_inv[r] - 1.0) + psi_Ex[p, q, r]);
    });
    std::for_each(std::execution::par, grid_Ey.begin(), grid_Ey.end(), [&](auto tuple) {
        auto [p, q, r] = tuple;
        size_t gr = r_ - 2 - r; 

        double dHx_dz = Hx_[p, q, gr] - Hx_[p, q, gr - 1];
        psi_Ey[p, q, r] = b[r] * psi_Ey[p, q, r] + c[r] * dHx_dz;
        Ey_[p, q, gr] += EyHzHx_[p, q, gr] * (dHx_dz * (k_inv[r] - 1.0) + psi_Ey[p, q, r]);
    });

}
