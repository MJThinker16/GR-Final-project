// Geodesic.cpp
// Geodesic.cpp
#include "Geodesic.hpp"
#include <cmath>
#include <iostream>

static constexpr double PI = 3.141592653589793238462643383279502884197169399375105820974;

Geodesic::Geodesic(double length, double mass, double throatRadius)
  : m_length(length), m_mass(mass), m_rho(throatRadius) {}

TraceResult Geodesic::traceEquatorial(
    double x, double y, double z,
    double dirx, double diry, double dirz
) const {
    // Integration parameters
    const double h     = 0.01;    // step size
    const int    steps = 5000;   // RK4 steps

    double l     = std::sqrt(x*x + y*y + z*z);
    double theta = PI / 2.0;
    double phi   = std::atan2(y, x);   // 四象限 atan2

    bool outside = std::abs(l) > m_length;

    double X = 2.0 * (std::abs(l) - m_length) / (PI * m_mass);
    double r = m_rho + (outside
        ? m_mass * (X * std::atan(X) - 0.5 * std::log(1 + X*X))
        : 0.0);

    double n_l      = dirx;
    double n_theta  = -dirz;
    double n_phi    = diry;

    double p_l0     = n_l;
    double p_theta0 = r * n_theta;
    double p_phi    = r * n_phi;
    double B_sqr    = p_theta0*p_theta0+p_phi * p_phi;
    double b        = p_phi;

    double yState[5] = { l, theta, phi, p_l0, p_theta0 };
    double t = 0.0;

    for (int i = 0; i < steps; ++i) {
        double yNext[5];
        runge_kutta_step(
            t, yState, h, yNext,
            m_length, m_mass, m_rho,
            B_sqr, b
        );
        // advance all 5 components
        for (int j = 0; j < 5; ++j) {
            yState[j] = yNext[j];
        }
        t += h;
    }

    // {l, θ, φ, p_l, p_θ}
    return { yState[0], yState[1], yState[2], yState[3], yState[4] };
}

void Geodesic::geodesic_odes_2d(
    double /*t*/, const double y[5], double dydt[5],
    double length, double mass, double rho,
    double B_sqr, double b
) {
    double l       = y[0];
    bool   outside = std::abs(l) > length;
    double X       = 2.0 * (std::abs(l) - length) / (PI * mass);
    double r       = rho + (outside
        ? mass * (X * std::atan(X) - 0.5 * std::log(1 + X*X))
        : 0.0);
    double dr_dl   = outside
        ? (std::atan(X) * 2.0 * l / (PI * std::abs(l)))
        : 0.0;

    // A.7a: dl/dt = p_l
    double dl_dt       = y[3];
    // A.7b: dθ/dt = p_θ / r^2
    double dtheta_dt   = y[4] / (r * r);
    // A.7c: dφ/dt = b / (r^2 * sin^2 θ)
    double dphi_dt     = b / (r * r * std::sin(y[1]) * std::sin(y[1]));
    // A.7d: dp_l/dt = B^2 * dr_dl / r^3
    double dpl_dt      = B_sqr * dr_dl / (r * r * r);
    // A.7e: dp_θ/dt = b^2 cosθ / (r^2 sin^3 θ)
    double dptheta_dt  = b * b * std::cos(y[1])
                       / (r * r * std::pow(std::sin(y[1]), 3));

    dydt[0] = dl_dt;
    dydt[1] = dtheta_dt;
    dydt[2] = dphi_dt;
    dydt[3] = dpl_dt;
    dydt[4] = dptheta_dt;
    //std::cout<<"l = "<<l<<"dl_dt = "<<dl_dt<<"dpl_dt = "<<dpl_dt<<"\n";
}

void Geodesic::runge_kutta_step(
    double t, const double y[5], double h, double yOut[5],
    double length, double mass, double rho,
    double B_sqr, double b
) {
    double k1[5], k2[5], k3[5], k4[5], yt[5];

    // k1 = h * f(t, y)
    geodesic_odes_2d(t, y, k1, length, mass, rho, B_sqr, b);
    for (int i = 0; i < 5; ++i) k1[i] *= h;

    // k2 = h * f(t+h/2, y+k1/2)
    for (int i = 0; i < 5; ++i) yt[i] = y[i] + 0.5 * k1[i];
    geodesic_odes_2d(t + 0.5*h, yt, k2, length, mass, rho, B_sqr, b);
    for (int i = 0; i < 5; ++i) k2[i] *= h;

    // k3 = h * f(t+h/2, y+k2/2)
    for (int i = 0; i < 5; ++i) yt[i] = y[i] + 0.5 * k2[i];
    geodesic_odes_2d(t + 0.5*h, yt, k3, length, mass, rho, B_sqr, b);
    for (int i = 0; i < 5; ++i) k3[i] *= h;

    // k4 = h * f(t+h, y+k3)
    for (int i = 0; i < 5; ++i) yt[i] = y[i] + k3[i];
    geodesic_odes_2d(t + h, yt, k4, length, mass, rho, B_sqr, b);
    for (int i = 0; i < 5; ++i) k4[i] *= h;

    // Combine into yOut
    for (int i = 0; i < 5; ++i) {
        yOut[i] = y[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
    }
}


/*
TraceResult Geodesic::traceEquatorial(double phi0, double r0) const {
    // Integration parameters
    const double h = 0.01;        // step size in affine parameter
    const int    steps = 10000;    // number of RK4 steps

    // Initial momentum components for equatorial
    double l = r0;
    bool outside = std::abs(l) > m_length;

    double x = 2.0 * (std::abs(l) - m_length) / (M_PI * m_mass);
    double r = m_rho + (outside
        ? m_mass * (x * std::atan(x) - 0.5 * std::log(1 + x*x))
        : 0.0);

    
    double n_l   = -std::cos(phi0);
    double n_phi = -std::sin(phi0);
    double p_l0  = n_l;
    double p_phi = r * n_phi;
    double B_sqr = p_phi * p_phi;
    double b     = p_phi;

    // State vector y = {l, phi, p_l}
    double y[3] = { r0, phi0, p_l0 };
    double t = 0.0;

    for(int i = 0; i < steps; ++i) {
        double y_next[3];
        runge_kutta_step(t, y, h, y_next,
                         m_length, m_mass, m_rho,
                         B_sqr, b);
        // Advance
        y[0] = y_next[0];
        y[1] = y_next[1];
        y[2] = y_next[2];
        t    += h;
    }

    return { y[1], y[0],y[2] };
}

void Geodesic::geodesic_odes_2d(
    double t, const double y[3], double dydt[3],
    double length, double mass, double rho,
    double B_sqr, double b
) {
    double l = y[0];
    bool outside = std::abs(l) > length;

    double x = 2.0 * (std::abs(l) - length) / (M_PI * mass);
    double r = rho + (outside
        ? mass * (x * std::atan(x) - 0.5 * std::log(1 + x*x))
        : 0.0);

    double dr_dl = outside
        ? (std::atan(x) * 2.0 * l / (M_PI * std::abs(l)))
        : 0.0;

    // A.7a: dl/dt = p_l
    double dl_dt = y[2];
    // A.7c: dphi/dt = b / r^2
    double dphi_dt = b / (r * r);
    // A.7d: dp_l/dt = B^2 * dr_dl / r^3
    double dpl_dt = B_sqr * dr_dl / (r * r * r);

    dydt[0] = dl_dt;
    dydt[1] = dphi_dt;
    dydt[2] = dpl_dt;
}

void Geodesic::runge_kutta_step(
    double t, const double y[3], double h, double y_out[3],
    double length, double mass, double rho,
    double B_sqr, double b
) {
    double k1[3], k2[3], k3[3], k4[3];
    double yt[3];

    // k1 = h * f(t, y)
    geodesic_odes_2d(t, y, k1, length, mass, rho, B_sqr, b);
    for(int i = 0; i < 3; ++i) k1[i] *= h;

    // k2 = h * f(t+h/2, y+k1/2)
    for(int i = 0; i < 3; ++i) yt[i] = y[i] + 0.5 * k1[i];
    geodesic_odes_2d(t + 0.5*h, yt, k2, length, mass, rho, B_sqr, b);
    for(int i = 0; i < 3; ++i) k2[i] *= h;

    // k3 = h * f(t+h/2, y+k2/2)
    for(int i = 0; i < 3; ++i) yt[i] = y[i] + 0.5 * k2[i];
    geodesic_odes_2d(t + 0.5*h, yt, k3, length, mass, rho, B_sqr, b);
    for(int i = 0; i < 3; ++i) k3[i] *= h;

    // k4 = h * f(t+h, y+k3)
    for(int i = 0; i < 3; ++i) yt[i] = y[i] + k3[i];
    geodesic_odes_2d(t + h, yt, k4, length, mass, rho, B_sqr, b);
    for(int i = 0; i < 3; ++i) k4[i] *= h;

    // Combine
    for(int i = 0; i < 3; ++i)
        y_out[i] = y[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.0;
}


*/
