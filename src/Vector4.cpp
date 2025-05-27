// Vector4.cpp
#include "Vector4.hpp"

Vector4::Vector4()
    : x(0.0), y(0.0), z(0.0), w(0.0) {}

Vector4::Vector4(double _x, double _y, double _z, double _w)
    : x(_x), y(_y), z(_z), w(_w) {}

Vector4 Vector4::operator+(const Vector4& o) const {
    return { x + o.x, y + o.y, z + o.z, w + o.w };
}

Vector4 Vector4::operator-(const Vector4& o) const {
    return { x - o.x, y - o.y, z - o.z, w - o.w };
}

Vector4 Vector4::operator*(double s) const {
    return { x * s, y * s, z * s, w * s };
}

Vector4 Vector4::operator/(double s) const {
    return { x / s, y / s, z / s, w / s };
}

double Vector4::dot(const Vector4& o) const {
    return x * o.x + y * o.y + z * o.z;
}

Vector4 Vector4::cross(const Vector4& o) const {
    return {
        y * o.z - z * o.y,
        z * o.x - x * o.z,
        x * o.y - y * o.x,
        0.0
    };
}

double Vector4::length() const {
    return std::sqrt(dot(*this));
}

Vector4 Vector4::normalize() const {
    double len = length();
    if (len == 0.0) return {0,0,0,w};
    return *this / len;
}
