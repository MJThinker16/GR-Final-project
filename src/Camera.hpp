// Camera.hpp
// Camera.cpp
/*Code contains following ability
1. Create Camera object with camera location, camera veiwing point, FoV angle
*/
#pragma once

#include "Vector4.hpp"
#include "Matrix4x4.hpp"


class Camera{
public:
    Camera(int width, int height, double fovDeg);
    ~Camera() = default;
    void lookAt(const Vector4& eye, const Vector4& target, const Vector4& up); //Set Camera pos, pointing, frame
    void setFOV(double fovDeg);
    void generateRay(        
        double px, double py,
        Vector4& origin,
        Vector4& dir,
        Matrix4x4& camToWorld,
        Matrix4x4& worldToCam
    )const;
private:
    int     m_width;
    int     m_height;
    double  m_fovX;
    double  m_fovY;    

    Vector4 m_position;
    Vector4 m_lookDir;
    Vector4 m_right;
    Vector4 m_up;
};



/*

#pragma once

#include "Vector4.hpp"
#include "Matrix4x4.hpp"

class Camera {
public:
    Camera(int width, int height, double fovDeg);
    ~Camera() = default;

    // 設定相機位置、目標與上方向量
    void lookAt(const Vector4& eye, const Vector4& target, const Vector4& up);
    // 設定水平視野角度（度）
    void setFOV(double fovDeg);

    // 產生一條光線，並同時取得 camera→world 與 world→camera 矩陣
    // px, py ∈ [-1,1] 對應螢幕歸一化座標
    void generateRay(
        double px, double py,
        Vector4& origin,
        Vector4& dir,
        Matrix4x4& camToWorld,
        Matrix4x4& worldToCam
    ) const;

private:
    double  m_fovX;
    double  m_fovY;

    Vector4 m_position;
    Vector4 m_lookDir;
    Vector4 m_up;
    Vector4 m_right;
};
*/