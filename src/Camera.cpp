// Camera.cpp
/*Code contains following ability
1. Create Camera object with camera location, camera veiwing point
2.
*/

#include "Camera.hpp"
#include <cmath>
#include<iostream>
#include<iomanip>
double M_PI1 = 3.14;

Camera::Camera(int width, int height, double fovDeg)
    : m_width(width), m_height(height) {
    setFOV(fovDeg);
}

void Camera::lookAt(const Vector4& eye, const Vector4& target, const Vector4& up) {
    m_position = eye;
    m_lookDir  = (target - eye).normalize();
    m_right    = m_lookDir.cross(up).normalize();
    m_up       = m_right.cross(m_lookDir).normalize();
    std::cout <<std::setprecision(2)<< std::fixed<<"lookDir ("<<m_lookDir.x<<","<<m_lookDir.y<<","<<m_lookDir.z<<") \n";
    std::cout <<std::setprecision(2)<< std::fixed<<"lookDiright ("<<m_right.x<<","<<m_right.y<<","<<m_right.z<<") \n";
    std::cout <<std::setprecision(2)<< std::fixed<<"lookDiup ("<<m_up.x<<","<<m_up.y<<","<<m_up.z<<") \n";

}

void Camera::setFOV(double fovDeg) {
    m_fovX = fovDeg * (M_PI1 / 180.0);
    m_fovY = m_fovX * (static_cast<double>(m_height) / m_width);
}

void Camera::generateRay(
    double px, double py,
    Vector4& origin,
    Vector4& dir,
    Matrix4x4& worldToCam,
    Matrix4x4& camToWorld
) const {
    // 光線起點
    origin = Vector4(m_position.x, m_position.y, m_position.z, 0.0);

    // 計算 camera 空間的偏移
    Vector4 r = m_right * (px * std::tan(0.5 * m_fovX));//小角度近似，應該要乘以tan(m_fovX/2)
    Vector4 u = m_up    * (py * std::tan(0.5 * m_fovY));
    Vector4 d = (m_lookDir + r + u).normalize();
    dir = Vector4(d.x, d.y, d.z, 0.0);

    // 構造 camera→world 矩陣
    camToWorld = Matrix4x4::identity();
	Vector4 new_x, new_y, new_z(0, 0, 1, 0);
	new_x.x = -origin.x;
	new_x.y = -origin.y;
	new_x.z = -origin.z;
	new_x.w = 0.0;
	new_x = new_x.normalize();
    
    double tmp = new_x.dot(dir);

	if (abs(abs(tmp) - 1.0) < 1e-8) // BUG: nan if ray direction and z axis are collinear
	{
		new_y = new_x.cross(new_z);
        new_z = new_y.cross(new_x);
	}
	else
	{
		new_z = dir.cross(new_x);
		new_y = new_x.cross(new_z);
	}

	new_y = new_y.normalize();
	new_y = new_y.normalize();


    camToWorld.m[0][0] = new_x.x;
    camToWorld.m[0][1] = new_x.y;
    camToWorld.m[0][2] = new_x.z;
    camToWorld.m[0][3] = 0;

    camToWorld.m[1][0] = new_y.x;
    camToWorld.m[1][1] = new_y.y;
    camToWorld.m[1][2] = new_y.z;
    camToWorld.m[1][3] = 0;

    camToWorld.m[2][0] = new_z.x;
    camToWorld.m[2][1] = new_z.y;
    camToWorld.m[2][2] = new_z.z;
    camToWorld.m[2][3] = 0;

    camToWorld.m[3][0] = 0.0;
    camToWorld.m[3][1] = 0.0;
    camToWorld.m[3][2] = 0.0;
    camToWorld.m[3][3] = 1.0;

    // 計算 world→camera 矩陣
    worldToCam = camToWorld.inverse(); //may error
}
/*
void build_camera_local_transformation(Math::Vector4 const &ray_origin, Math::Vector4 const &ray_dir, Math::Matrix4x4 &ans, Math::Matrix4x4 &inverse)
{
	using namespace Math;
	Vector4 new_x, new_y, new_z(0, 0, 1, 0);
	new_x.x = -ray_origin.x;
	new_x.y = -ray_origin.y;
	new_x.z = -ray_origin.z;
	new_x.w = 0.0;
	Vector4Normalize(new_x);

	Vector4 tmp;
	Vector4Dot(new_x, ray_dir, tmp);

	if (abs(abs(tmp.x) - 1.0) < 1e-8) // BUG: nan if ray direction and z axis are collinear
	{
		Vector4Cross(new_x, new_z, new_y);
		Vector4Cross(new_y, new_x, new_z);
	}
	else
	{
		Vector4Cross(ray_dir, new_x, new_z);
		Vector4Cross(new_x, new_z, new_y);
	}


	Vector4Normalize(new_y);
	Vector4Normalize(new_z);

	ans._11 = new_x.x; ans._12 = new_x.y; ans._13 = new_x.z; ans._14 = 0;
	ans._21 = new_y.x; ans._22 = new_y.y; ans._23 = new_y.z; ans._24 = 0;
	ans._31 = new_z.x; ans._32 = new_z.y; ans._33 = new_z.z; ans._34 = 0;
	ans._41 = 0; ans._42 = 0; ans._43 = 0; ans._44 = 1.0;

	Matrix4x4Inverse(ans, inverse);
}
*/