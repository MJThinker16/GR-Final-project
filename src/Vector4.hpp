// Vector4.hpp
#pragma once
#include <cmath>

struct Vector4 {
    double x, y, z, w;

    // Constructors
    Vector4();
    Vector4(double _x, double _y, double _z, double _w = 0.0);

    // Arithmetic operators
    Vector4 operator+(const Vector4& o) const;
    Vector4 operator-(const Vector4& o) const;
    Vector4 operator*(double s) const;
    Vector4 operator/(double s) const;

    // Dot and cross (w ignored for dot, w set to 0 for cross)
    double dot(const Vector4& o) const;
    Vector4 cross(const Vector4& o) const;

    // Length (Euclidean) and normalization
    double length() const;
    Vector4 normalize() const;
};

