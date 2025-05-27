// Matrix4x4.cpp
#include "Matrix4x4.hpp"

Vector4 Matrix4x4::mul(const Vector4& v) const {
    return {
        m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z + m[0][3]*v.w,
        m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z + m[1][3]*v.w,
        m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z + m[2][3]*v.w,
        m[3][0]*v.x + m[3][1]*v.y + m[3][2]*v.z + m[3][3]*v.w
    };
}

Matrix4x4 Matrix4x4::operator*(const Matrix4x4& o) const {
    Matrix4x4 R;
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            R.m[i][j] = 0.0;
            for(int k = 0; k < 4; ++k)
                R.m[i][j] += m[i][k] * o.m[k][j];
        }
    }
    return R;
}

Matrix4x4 Matrix4x4::inverse() const {
    // 假設為仿射矩陣，分離 R (3x3) 和 t
    Matrix4x4 inv = identity();
    // 轉置前3x3作為逆
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            inv.m[i][j] = m[j][i];
    // 平移分量
    Vector4 t{m[0][3], m[1][3], m[2][3], 1.0};
    // inv translation = -R^T * t
    Vector4 t2 = {0,0,0,0};
    t2.x = -(inv.m[0][0]*t.x + inv.m[0][1]*t.y + inv.m[0][2]*t.z);
    t2.y = -(inv.m[1][0]*t.x + inv.m[1][1]*t.y + inv.m[1][2]*t.z);
    t2.z = -(inv.m[2][0]*t.x + inv.m[2][1]*t.y + inv.m[2][2]*t.z);
    inv.m[0][3] = t2.x;
    inv.m[1][3] = t2.y;
    inv.m[2][3] = t2.z;
    // 下行保持為 0,0,0,1
    inv.m[3][0] = inv.m[3][1] = inv.m[3][2] = 0.0;
    inv.m[3][3] = 1.0;
    return inv;
}

