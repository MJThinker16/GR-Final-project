// Matrix4x4.hpp
#pragma once

#include "Vector4.hpp"

// 4x4 同質座標矩陣 (Affine: rotation + translation)
struct Matrix4x4 {
    double m[4][4];

    // 建立單位矩陣
    static Matrix4x4 identity() {
        Matrix4x4 I;
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
                I.m[i][j] = (i == j ? 1.0 : 0.0);
        return I;
    }

    // 矩陣乘向量 (同質座標)
    Vector4 mul(const Vector4& v) const;

    // 矩陣相乘
    Matrix4x4 operator*(const Matrix4x4& o) const;

    // 反矩陣 (僅針對 rotation + translation)
    Matrix4x4 inverse() const;
};
