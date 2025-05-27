
// src/main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <omp.h>
#include "Vector4.hpp"
#include "Matrix4x4.hpp"
#include "Camera.hpp"
#include "Geodesic.hpp"
double M_PI3 = 3.141592653589793238462643383279502884197169399375105820974;

int main(){
    
/*
    Geodesic geo1(0.5, 0.5, 1);
    std::cout<<"res1\n";
    
    auto res1 = geo1.traceEquatorial(4.0, 0.0, 0.0, std::cos(M_PI3/6), std::sin(M_PI3/6),0);
    std::cout<<"res2\n";
    auto res2 = geo1.traceEquatorial(4.0, 0.0, 0.0, std::cos(M_PI3/6), -std::sin(M_PI3/6),0);
    std::cout<<"res3\n";
    auto res3 = geo1.traceEquatorial(4.0, 0.0, 0.0, std::cos(M_PI3/6), 0,-std::sin(M_PI3/6));
    std::cout<<"res4\n";
    auto res4 = geo1.traceEquatorial(4.0, 0.0, 0.0, std::cos(M_PI3/6), 0,std::sin(M_PI3/6));

    std::cout<<res1.l_out<<" "<<res2.l_out<<" "<<res3.l_out<<" "<<res4.l_out;
    return 0;
*/
    // —— 1. 硬編初始條件 —— 
    const int    width          = 1000;
    const int    height         = 1000;
    const double fovDeg         = 100.0;

    const double whLength       = 0.5;
    const double whMass         = 5;
    const double whThroatRadius = 1;
    
    const Vector4 eye           = { 6.5*whThroatRadius+whLength, 0.0, 0.0, 0.0 };
    const Vector4 lookAt        = { 0.0, 0.0, 0.0, 0.0 };
    const Vector4 upVec         = { 0.0, 0.0, 1.0, 0.0 };

    const int    samples        = 3;
    const int    threads        = omp_get_max_threads();
    const std::string outCsv    = "result.csv";
    
    // —— 2. 初始化模組 —— 
    Camera   cam(width, height, fovDeg);
    cam.lookAt(eye, lookAt, upVec);

    Geodesic geo(whLength, whMass, whThroatRadius);

    // 緩衝空間
    // —— 3. 緩衝空間：phi、l、dirX、dirY、dirZ —— 
    std::vector<double> phiBuf (width * height);
    std::vector<double> thetaBuf (width * height);
    std::vector<double> lBuf   (width * height);

    const double wUnit = 2.0 / width;
    const double hUnit = 2.0 / height;

    // —— 3. 平行渲染與數據累積 —— 
    omp_set_num_threads(threads);
    #pragma omp parallel for
    for(int row=0; row<height; ++row) {
        std::mt19937_64 rng(row ^ 0xDEADBEEF);
        std::uniform_real_distribution<double> dist(0.0,1.0);

        for(int col=0; col<width; ++col) {
            double x0 = (col/(double)width)*2.0 - 1.0 - wUnit/2.0;
            double y0 = (row/(double)height)*2.0 - 1.0 - hUnit/2.0;
            //std::cout <<std::setprecision(2)<< std::fixed<<"("<<x0<<","<<y0<<") ";
            double accPhi = 0.0, accL = 0.0,acctheta = 0.0;

            for(int s=0; s<samples; ++s) {
                double xo = dist(rng)*wUnit;
                double yo = dist(rng)*hUnit;
                // 產生光線與座標變換矩陣
                Vector4 origin, dir;
                Matrix4x4 cam2w, w2c;
                cam.generateRay(x0+xo, y0+yo,
                                origin, dir,
                                cam2w, w2c); //這裡的direction 是 outcoming ray, 和origin 都在global coordinate.
                //std::cout <<std::setprecision(2)<< std::fixed<<"("<<dir.x<<","<<dir.y<<","<<dir.z<<") ";
                //std::cout <<std::setprecision(2)<< std::fixed<<"("<<dirCam.x<<","<<dirCam.y<<","<<dirCam.z<<") ";
   
                // 赤道面映射
                auto res = geo.traceEquatorial(origin.x, origin.y, origin.z, dir.x, dir.y, dir.z);
                accPhi += res.phi_out;
                acctheta += res.theta_out;

                accL   += res.l_out;
                //std::cout <<std::setprecision(2)<< std::fixed<<"("<<res.l_out<<","<<res.theta_out<<","<<res.phi_out<<") ";

            }

            int idx = row * width + col;
            // 平均
            phiBuf[idx] = accPhi / samples;
            lBuf  [idx] = accL   / samples;
            thetaBuf  [idx] = acctheta   / samples;

        }
        //std::cout<<std::endl;
    }

    // —— 5. 批次寫入 CSV (含 dx,dy,dz) —— 
    std::ofstream ofs(outCsv);
    if(!ofs) {
        std::cerr << "無法建立 " << outCsv << "\n";
        return 1;
    }
    ofs << "row,col,l_out,theta_out,phi_out\n";
    for(int row = 0; row < height; ++row) {
        for(int col = 0; col < width; ++col) {
            int idx = row*width + col;
            ofs << row << ',' << col << ','
                << lBuf[idx] << ',' << thetaBuf[idx] << ','
                << phiBuf[idx]   << '\n';
        }
    }
    ofs.close();

    std::cout << "完成CSV 已儲存於 " << outCsv << "\n";
    return 0;
}






















/*
// src/main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <omp.h>
#include "Vector4.hpp"
#include "Matrix4x4.hpp"
#include "Camera.hpp"
#include "Geodesic.hpp"
double M_PI3 = 3.141592653589793238462643383279502884197169399375105820974;

int main(){


    // —— 1. 硬編初始條件 —— 
    const int    width          = 8;
    const int    height         = 8;
    const double fovDeg         = 100.0;

    const double whLength       = 0.5;
    const double whMass         = 0.5;
    const double whThroatRadius = 1;
    
    const Vector4 eye           = { 6.5*whThroatRadius+whLength, 0.0, 0.0, 0.0 };
    const Vector4 lookAt        = { 0.0, 0.0, 0.0, 0.0 };
    const Vector4 upVec         = { 0.0, 0.0, 1.0, 0.0 };

    const int    samples        = 1;
    const int    threads        = omp_get_max_threads();
    const std::string outCsv    = "result.csv";
    
    // —— 2. 初始化模組 —— 
    Camera   cam(width, height, fovDeg);
    cam.lookAt(eye, lookAt, upVec);

    Geodesic geo(whLength, whMass, whThroatRadius);

    // 緩衝空間
    // —— 3. 緩衝空間：phi、l、dirX、dirY、dirZ —— 
    std::vector<double> phiBuf (width * height);
    std::vector<double> lBuf   (width * height);
    std::vector<double> dirX   (width * height);
    std::vector<double> dirY   (width * height);
    std::vector<double> dirZ   (width * height);

    const double wUnit = 2.0 / width;
    const double hUnit = 2.0 / height;

    // —— 3. 平行渲染與數據累積 —— 
    omp_set_num_threads(threads);
    //#pragma omp parallel for
    for(int row=0; row<height; ++row) {
        std::mt19937_64 rng(row ^ 0xDEADBEEF);
        std::uniform_real_distribution<double> dist(0.0,1.0);

        for(int col=0; col<width; ++col) {
            double x0 = (col/(double)width)*2.0 - 1.0 + wUnit/2.0;
            double y0 = (row/(double)height)*2.0 - 1.0 + hUnit/2.0;
            //std::cout <<std::setprecision(2)<< std::fixed<<"("<<x0<<","<<y0<<") ";
            double accPhi = 0.0, accL = 0.0;
            double avgDx = 0.0, avgDy = 0.0, avgDz = 0.0;

            for(int s=0; s<samples; ++s) {
                double xo = 0;//dist(rng)*wUnit;
                double yo = 0;//dist(rng)*hUnit;
                // 產生光線與座標變換矩陣
                Vector4 origin, dir;
                Matrix4x4 cam2w, w2c;
                cam.generateRay(x0+xo, y0+yo,
                                origin, dir,
                                cam2w, w2c); //這裡的direction 是 outcoming ray, 和origin 都在global coordinate.
                //std::cout <<std::setprecision(2)<< std::fixed<<"("<<dir.x<<","<<dir.y<<","<<dir.z<<") ";
                Vector4 dirCam = (w2c.mul(dir)).normalize(); //這裡的direction 是 outcoming ray, 和origin 都在global coordinate.
                //std::cout <<std::setprecision(2)<< std::fixed<<"("<<dirCam.x<<","<<dirCam.y<<","<<dirCam.z<<") ";
                std::cout <<std::setprecision(2)<< std::fixed<<"("<<dir.x<<","<<dirCam.x<<") ";
   
                double phi0 = std::atan2(dirCam.y, dirCam.x);
                double r0   = std::sqrt(origin.x*origin.x
                                      + origin.y*origin.y
                                      + origin.z*origin.z);

                // 赤道面映射
                auto res = geo.traceEquatorial(phi0, r0);
                accPhi += res.phi_out;
                accL   += res.l_out;
                //std::cout <<std::setprecision(2)<< std::fixed<<"("<<res.phi_out<<","<<res.l_out<<") ";

                // 重建 local 出射方向、轉到 global
                Vector4 dirCamOut{
                    std::cos(res.phi_out),
                    std::sin(res.phi_out),
                     0.0, 0.0
                };
                Vector4 dirWorld = cam2w.mul(dirCamOut);
                //std::cout <<std::setprecision(2)<< std::fixed<<"("<<dirWorld.x<<","<<dirWorld.y<<","<<dirWorld.z<<") ";

                avgDx += dirWorld.x;
                avgDy += dirWorld.y;
                avgDz += dirWorld.z;
            }

            int idx = row * width + col;
            // 平均
            phiBuf[idx] = accPhi / samples;
            lBuf  [idx] = accL   / samples;
            dirX  [idx] = avgDx  / samples;
            dirY  [idx] = avgDy  / samples;
            dirZ  [idx] = avgDz  / samples;
        }
        std::cout<<std::endl;
    }

    // —— 5. 批次寫入 CSV (含 dx,dy,dz) —— 
    std::ofstream ofs(outCsv);
    if(!ofs) {
        std::cerr << "無法建立 " << outCsv << "\n";
        return 1;
    }
    ofs << "row,col,phi_out,l_out,dx,dy,dz\n";
    for(int row = 0; row < height; ++row) {
        for(int col = 0; col < width; ++col) {
            int idx = row*width + col;
            ofs << row << ',' << col << ','
                << phiBuf[idx] << ',' << lBuf[idx] << ','
                << dirX[idx]   << ',' << dirY[idx]   << ','
                << dirZ[idx]   << '\n';
        }
    }
    ofs.close();

    std::cout << "完成CSV 已儲存於 " << outCsv << "\n";
    return 0;
}
*/