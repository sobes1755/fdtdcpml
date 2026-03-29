#pragma once

struct CPML {

  size_t size_E1H2_;
  size_t size_E2H1_;

  std::vector<double> b_E_, c_E_, k_inv_E_;
  std::vector<double> b_H_, c_H_, k_inv_H_;

  std::vector<double> raw_psi_E1_, raw_psi_E2_;
  std::vector<double> raw_psi_H1_, raw_psi_H2_;

  std::mdspan<double, std::dextents<size_t, 3>> psi_E1_, psi_E2_;
  std::mdspan<double, std::dextents<size_t, 3>> psi_H1_, psi_H2_;

  CPML(size_t size1, size_t size2, size_t size3, double step3, double stepT) :

      size_E1H2_((size1 - 1) * size2 * size3),
      size_E2H1_(size1 * (size2 - 1) * size3),

      b_E_(size3), c_E_(size3), k_inv_E_(size3),
      b_H_(size3), c_H_(size3), k_inv_H_(size3),

      raw_psi_E1_(size_E1H2_), raw_psi_E2_(size_E2H1_),
      raw_psi_H1_(size_E2H1_), raw_psi_H2_(size_E1H2_),

      psi_E1_(raw_psi_E1_.data(), size1 - 1, size2, size3),
      psi_E2_(raw_psi_E2_.data(), size1, size2 - 1, size3),
      psi_H1_(raw_psi_H1_.data(), size1, size2 - 1, size3),
      psi_H2_(raw_psi_H2_.data(), size1 - 1, size2, size3)

  {

    double R = 1e-6;  // dB = 10 * log_10(R^2) = -120
    double depth = size3 * step3;

    double deg = 3.;
    double kappa_max = 3.0;
    double alpha_max = 0.1;  // complex frequency shift (CFS)
    double sigma_max = 0.8 * ((deg + 1.0) * std::log(1.0 / R)) / (2.0 * fdtdcpml::eta0 * depth);

    // pqr == 0 -- вблизи PEC
    // pqr == size3 - 1 -- вблизи расчётной области

    for (size_t pqr = 0; pqr < size3; ++pqr) {

      double rho = static_cast<double>(size3 - pqr) / static_cast<double>(size3);
      double kappa = 1.0 + (kappa_max - 1.0) * std::pow(rho, deg);
      double alpha = alpha_max * (1.0 - rho);
      double sigma = sigma_max * std::pow(rho, deg);

      b_E_[pqr] = std::exp(-(sigma / kappa + alpha) * stepT / fdtdcpml::eps0);
      c_E_[pqr] = sigma * (b_E_[pqr] - 1.0) / (sigma * kappa + alpha * std::pow(kappa, 2));
      k_inv_E_[pqr] = 1.0 / kappa;

    }

    for (size_t pqr = 0; pqr < size3; ++pqr) {

      double rho = static_cast<double>(size3 - pqr - 0.5) / static_cast<double>(size3);
      double kappa = 1.0 + (kappa_max - 1.0) * std::pow(rho, deg);
      double alpha = alpha_max * (1.0 - rho);
      double sigma = sigma_max * std::pow(rho, deg);

      b_H_[pqr] = std::exp(-(sigma / kappa + alpha) * stepT / fdtdcpml::eps0);
      c_H_[pqr] = sigma * (b_H_[pqr] - 1.0) / (sigma * kappa + alpha * std::pow(kappa, 2));
      k_inv_H_[pqr] = 1.0 / kappa;

    }

  }

};
