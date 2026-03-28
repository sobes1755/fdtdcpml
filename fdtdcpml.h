#include <vector>
#include <experimental/mdspan>

#include "fdtd.h"
#include "cpml.h"
#include "mdspanio.h"

class FDTDCPML {
private:

    // FDTD interior region grid sizes...

    std::size_t p_int_;
    std::size_t q_int_;
    std::size_t r_int_;
    std::size_t s_int_;

    // FDTD absorbtion (CPML) region grid depths...

    std::size_t p_abs_l_, p_abs_h_;
    std::size_t q_abs_l_, q_abs_h_;
    std::size_t r_abs_l_, r_abs_h_;    

    // FDTD grid (both interior region and absorbtion region) sizes...

    std::size_t p_;
    std::size_t q_;
    std::size_t r_;
    std::size_t s_; 

    // FDTD interior region size (physical)...

    double x_min_, x_max_;
    double y_min_, y_max_;
    double z_min_, z_max_;
    double t_min_, t_max_;

    // FDTD spatial and temporal step sizes (physical)...

    double x_step_;
    double y_step_;
    double z_step_;
    double t_step_;

    // Memory for H, E and H, E update coefficients...

    using GridRaw = std::vector<double>;

    GridRaw rawHx_, rawHxHx_, rawHxEyEz_;
    GridRaw rawHy_, rawHyHy_, rawHyEzEx_;
    GridRaw rawHz_, rawHzHz_, rawHzExEy_;
    GridRaw rawEx_, rawExEx_, rawExHyHz_;
    GridRaw rawEy_, rawEyEy_, rawEyHzHx_;
    GridRaw rawEz_, rawEzEz_, rawEzHxHy_;

    // Views for H and E componenets...

    using GridView = std::mdspan<double, std::dextents<size_t, 3>>;

    GridView Hx_, HxHx_, HxEyEz_;
    GridView Hy_, HyHy_, HyEzEx_;
    GridView Hz_, HzHz_, HzExEy_;
    GridView Ex_, ExEx_, ExHyHz_;
    GridView Ey_, EyEy_, EyHzHx_;
    GridView Ez_, EzEz_, EzHxHy_;

    // Memory (vector<double>) and views for CPML psi values (absorbtion region)...

    CPML psi_x_l_;
    CPML psi_x_h_;
    CPML psi_y_l_;
    CPML psi_y_h_;
    CPML psi_z_l_;
    CPML psi_z_h_;

    //

    double t_ = 0.0;

public:

    // Getters...

    double x_step() const noexcept { return x_step_; }
    double y_step() const noexcept { return y_step_; }
    double z_step() const noexcept { return z_step_; }
    double t_step() const noexcept { return t_step_; }

    double t() const noexcept { return t_; }

    auto Hx() const noexcept { return Hx_; }
    auto Hy() const noexcept { return Hy_; }
    auto Hz() const noexcept { return Hz_; }

    auto Ex() const noexcept { return Ex_; }
    auto Ey() const noexcept { return Ey_; }
    auto Ez() const noexcept { return Ez_; }

    //

    enum class ECOMP { Hx, Hy, Hz, Ex, Ey, Ez };
    enum class EPQRS { P, Q, R, S };
    enum class EXYZT { X, Y, Z, T };

    //

    using ScaledSpan3D = std::mdspan<double, std::dextents<size_t, 3>, std::layout_right, fdtdcpml::ScalingAccessor<double>>;
    using SlicedSpan2D = std::mdspan<double, std::dextents<size_t, 2>, std::layout_stride, fdtdcpml::ScalingAccessor<double>>;

    SlicedSpan2D slice(ECOMP cmp, EPQRS pqr, size_t v) {

        ScaledSpan3D scaled;
        SlicedSpan2D sliced;

        if (cmp == ECOMP::Hx) {
            scaled = std::mdspan(Hx_.data_handle(), Hx_.mapping(), fdtdcpml::ScalingAccessor<double>{x_step_});
        } else if (cmp == ECOMP::Hy) {
            scaled = std::mdspan(Hy_.data_handle(), Hy_.mapping(), fdtdcpml::ScalingAccessor<double>{y_step_});
        } else if (cmp == ECOMP::Hz) {
            scaled = std::mdspan(Hz_.data_handle(), Hz_.mapping(), fdtdcpml::ScalingAccessor<double>{z_step_});
        } else if (cmp == ECOMP::Ex) {
            scaled = std::mdspan(Ex_.data_handle(), Ex_.mapping(), fdtdcpml::ScalingAccessor<double>{x_step_});
        } else if (cmp == ECOMP::Ey) {
            scaled = std::mdspan(Ey_.data_handle(), Ey_.mapping(), fdtdcpml::ScalingAccessor<double>{y_step_});
        } else if (cmp == ECOMP::Ez) {
            scaled = std::mdspan(Ez_.data_handle(), Ez_.mapping(), fdtdcpml::ScalingAccessor<double>{z_step_});
        }

        if (pqr == EPQRS::P) {
            sliced = SlicedSpan2D(std::submdspan(scaled, v, std::full_extent, std::full_extent));
        } else if (pqr == EPQRS::Q) {
            sliced = SlicedSpan2D(std::submdspan(scaled, std::full_extent, v, std::full_extent));
        } else if (pqr == EPQRS::R) {
            sliced = SlicedSpan2D(std::submdspan(scaled, std::full_extent, std::full_extent, v));
        }

        return sliced;

    }

    // Calculate default time step...

    static constexpr double t_step_default(const double x_step, const double y_step, const double z_step) {

        double reciprocal = std::sqrt(
            1. / (x_step * x_step) +
            1. / (y_step * y_step) +
            1. / (z_step * z_step));

        return 1. / (fdtdcpml::c0 * reciprocal);

    }

    // Simple functions which describe medium properties...

    using FunTXYZ = double (*)(double t, double x, double y, double z);
    using FunT = double (*)(double t);
    using FunXYZ = double (*)(double x, double y, double z);

    FunXYZ sigma = [](double, double, double) -> double { return 0.0; };
    FunXYZ sigma_m = [](double, double, double) -> double { return 0.0; };
    FunXYZ mu_r = [](double, double, double) -> double { return 1.0; };
    FunXYZ eps_r = [](double, double, double) -> double { return 1.0; };

    // 

    FDTDCPML(
        double x_min = 0, double x_max = 1, std::size_t p_int = 257,
        double y_min = 0, double y_max = 1, std::size_t q_int = 257,
        double z_min = 0, double z_max = 1, std::size_t r_int = 257,
        double t_min = 0, double t_max = 1, std::size_t s_int = 257,
        std::size_t p_abs_l = 10, std::size_t p_abs_h = 10,
        std::size_t q_abs_l = 10, std::size_t q_abs_h = 10,
        std::size_t r_abs_l = 10, std::size_t r_abs_h = 10
    );

    //

    void init();
    void step();

    //

    struct XYZ { double x, y, z; };
    struct PQR { size_t p, q, r; };

    XYZ Hx_pqr2xyz(const PQR s) const noexcept {
        double x = x_min_ + (std::clamp(s.p, p_abs_l_, p_abs_l_ + p_int_) - p_abs_l_) * x_step_;
        double y = y_min_ + (std::clamp(s.q, q_abs_l_, q_abs_l_ + q_int_) - q_abs_l_ + 0.5) * y_step_;
        double z = z_min_ + (std::clamp(s.r, r_abs_l_, r_abs_l_ + r_int_) - r_abs_l_ + 0.5) * z_step_;
        return { x, y, z };
    }
    XYZ Hy_pqr2xyz(const PQR s) const noexcept {
        double x = x_min_ + (std::clamp(s.p, p_abs_l_, p_abs_l_ + p_int_) - p_abs_l_ + 0.5) * x_step_;
        double y = y_min_ + (std::clamp(s.q, q_abs_l_, q_abs_l_ + q_int_) - q_abs_l_) * y_step_;
        double z = z_min_ + (std::clamp(s.r, r_abs_l_, r_abs_l_ + r_int_) - r_abs_l_ + 0.5) * z_step_;
        return { x, y, z};
    }
    XYZ Hz_pqr2xyz(const PQR s) const noexcept {
        double x = x_min_ + (std::clamp(s.p, p_abs_l_, p_abs_l_ + p_int_) - p_abs_l_ + 0.5) * x_step_;
        double y = y_min_ + (std::clamp(s.q, q_abs_l_, q_abs_l_ + q_int_) - q_abs_l_ + 0.5) * y_step_;
        double z = z_min_ + (std::clamp(s.r, r_abs_l_, r_abs_l_ + r_int_) - r_abs_l_) * z_step_;
        return { x, y, z};
    }

    XYZ Ex_pqr2xyz(const PQR s) const noexcept {
        double x = x_min_ + (std::clamp(s.p, p_abs_l_, p_abs_l_ + p_int_) - p_abs_l_ + 0.5) * x_step_;
        double y = y_min_ + (std::clamp(s.q, q_abs_l_, q_abs_l_ + q_int_) - q_abs_l_) * y_step_;
        double z = z_min_ + (std::clamp(s.r, r_abs_l_, r_abs_l_ + r_int_) - r_abs_l_) * z_step_;
        return { x, y, z };
    }
    XYZ Ey_pqr2xyz(const PQR s) const noexcept {
        double x = x_min_ + (std::clamp(s.p, p_abs_l_, p_abs_l_ + p_int_) - p_abs_l_) * x_step_;
        double y = y_min_ + (std::clamp(s.q, q_abs_l_, q_abs_l_ + q_int_) - q_abs_l_ + 0.5) * y_step_;
        double z = z_min_ + (std::clamp(s.r, r_abs_l_, r_abs_l_ + r_int_) - r_abs_l_) * z_step_;
        return { x, y, z };
    }
    XYZ Ez_pqr2xyz(const PQR s) const noexcept {
        double x = x_min_ + (std::clamp(s.p, p_abs_l_, p_abs_l_ + p_int_) - p_abs_l_) * x_step_;
        double y = y_min_ + (std::clamp(s.q, q_abs_l_, q_abs_l_ + q_int_) - q_abs_l_) * y_step_;
        double z = z_min_ + (std::clamp(s.r, r_abs_l_, r_abs_l_ + r_int_) - r_abs_l_ + 0.5) * z_step_;
        return { x, y, z };
    }

    PQR Ex_xyz2pqr(const XYZ s) const noexcept {
        size_t p = static_cast<size_t>(std::round((s.x - x_min_) / x_step_ - 0.5)) + p_abs_l_;
        size_t q = static_cast<size_t>(std::round((s.y - y_min_) / y_step_)) + q_abs_l_;
        size_t r = static_cast<size_t>(std::round((s.z - z_min_) / z_step_)) + r_abs_l_;
        return { p, q, r };
    }
    PQR Ey_xyz2pqr(const XYZ s) const noexcept {
        size_t p = static_cast<size_t>(std::round((s.x - x_min_) / x_step_)) + p_abs_l_;
        size_t q = static_cast<size_t>(std::round((s.y - y_min_) / y_step_ - 0.5)) + q_abs_l_;
        size_t r = static_cast<size_t>(std::round((s.z - z_min_) / z_step_)) + r_abs_l_;
        return { p, q, r };
    }
    PQR Ez_xyz2pqr(const XYZ s) const noexcept {
        size_t p = static_cast<size_t>(std::round((s.x - x_min_) / x_step_)) + p_abs_l_;
        size_t q = static_cast<size_t>(std::round((s.y - y_min_) / y_step_)) + q_abs_l_;
        size_t r = static_cast<size_t>(std::round((s.z - z_min_) / z_step_ - 0.5)) + r_abs_l_;
        return { p, q, r };
    }

private:

    // Debug

    double Hx_min, Hx_max;
    double Hy_min, Hy_max;
    double Hz_min, Hz_max;
    double Ex_min, Ex_max;
    double Ey_min, Ey_max;
    double Ez_min, Ez_max;

    PQR Hx_min_pqr, Hx_max_pqr;
    PQR Hy_min_pqr, Hy_max_pqr;    
    PQR Hz_min_pqr, Hz_max_pqr;
    PQR Ex_min_pqr, Ex_max_pqr;
    PQR Ey_min_pqr, Ey_max_pqr;    
    PQR Ez_min_pqr, Ez_max_pqr;

    //

    void H_minmax() {
        std::tie(Hx_min, Hx_max, Hx_min_pqr, Hx_max_pqr) = view_minmax(rawHx_, Hx_); 
        std::tie(Hy_min, Hy_max, Hy_min_pqr, Hy_max_pqr) = view_minmax(rawHy_, Hy_);
        std::tie(Hz_min, Hz_max, Hz_min_pqr, Hz_max_pqr) = view_minmax(rawHz_, Hz_);
    }

    void E_minmax() {
        std::tie(Ex_min, Ex_max, Ex_min_pqr, Ex_max_pqr) = view_minmax(rawEx_, Ex_); 
        std::tie(Ey_min, Ey_max, Ey_min_pqr, Ey_max_pqr) = view_minmax(rawEy_, Ey_);
        std::tie(Ez_min, Ez_max, Ez_min_pqr, Ez_max_pqr) = view_minmax(rawEz_, Ez_);
    }

    std::tuple<double, double, PQR, PQR> view_minmax(GridRaw v1d, GridView v3d) {

        auto [ptr_min, ptr_max] = std::minmax_element(v1d.begin(), v1d.end());

        double v_min = *ptr_min;
        double v_max = *ptr_max;

        size_t offset_min = std::distance(v1d.begin(), ptr_min);
        size_t offset_max = std::distance(v1d.begin(), ptr_max);

        size_t p_min = offset_min / (v3d.extent(2) * v3d.extent(1));
        size_t q_min = (offset_min / v3d.extent(2)) % v3d.extent(1);
        size_t r_min = offset_min % v3d.extent(2);
        size_t p_max =  offset_max / (v3d.extent(2) * v3d.extent(1));
        size_t q_max = (offset_max / v3d.extent(2)) % v3d.extent(1);
        size_t r_max = offset_max % v3d.extent(2);

        return {v_min, v_max, PQR{p_min, q_min, r_min}, PQR{p_max, q_max, r_max}};

    }

    //

    void fdtd_init_H();
    void fdtd_init_E();

    void fdtd_H();
    void fdtd_E();

    void cpml_x_l_H(); void cpml_x_h_H();
    void cpml_y_l_H(); void cpml_y_h_H();
    void cpml_z_l_H(); void cpml_z_h_H();

    void cpml_x_l_E(); void cpml_x_h_E();
    void cpml_y_l_E(); void cpml_y_h_E();
    void cpml_z_l_E(); void cpml_z_h_E();

    void H_source();
    void E_source();

    void E_source_LG(double t, double x0, double y0, double z0, double freq, int n, int l, double w);

};
