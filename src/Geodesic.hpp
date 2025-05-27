// Geodesic.hpp
#pragma once

struct TraceResult {
    double l_out;
    double theta_out;
    double phi_out;
    double pl_out;    // 新增這個欄位
    double ptheta_out; 
};

class Geodesic {
public:
    // length: throat half-length a
    // mass: wormhole mass M
    // throatRadius: throat radius rho
    Geodesic(double length, double mass, double throatRadius);

    // Equatorial mapping: input initial phi0 and radial distance r0
    // returns phi_out (exit azimuth) and l_out (exit proper radius)
    TraceResult traceEquatorial(double x, double y, double z, double dirx,double diry, double dirz) const;

private:
    double m_length;
    double m_mass;
    double m_rho;

    // ODE for 2D equatorial geodesic (A.7a, A.7c, A.7d)
    static void geodesic_odes_2d(
        double /*t*/, const double y[5], double dydt[5],
        double length, double mass, double rho,
        double B_sqr, double b
    );

    // Single RK4 step for the 3D state y = {l, phi, p_l}
    static void runge_kutta_step(
        double t, const double y[5], double h, double y_out[5],
        double length, double mass, double rho,
        double B_sqr, double b
    );
};


/*
#pragma once

struct TraceResult {
    double phi_out;
    double l_out;
    double pl_out;   // 新增這個欄位

};

class Geodesic {
public:
    // length: throat half-length a
    // mass: wormhole mass M
    // throatRadius: throat radius rho
    Geodesic(double length, double mass, double throatRadius);

    // Equatorial mapping: input initial phi0 and radial distance r0
    // returns phi_out (exit azimuth) and l_out (exit proper radius)
    TraceResult traceEquatorial(double phi0, double r0) const;

private:
    double m_length;
    double m_mass;
    double m_rho;

    // ODE for 2D equatorial geodesic (A.7a, A.7c, A.7d)
    static void geodesic_odes_2d(
        double t, const double y[3], double dydt[3],
        double length, double mass, double rho,
        double B_sqr, double b
    );

    // Single RK4 step for the 3D state y = {l, phi, p_l}
    static void runge_kutta_step(
        double t, const double y[3], double h, double y_out[3],
        double length, double mass, double rho,
        double B_sqr, double b
    );
};






*/