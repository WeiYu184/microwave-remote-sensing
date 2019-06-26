// 类SatelliteOrbitFit：通过三次多项式拟合法实现卫星轨道参数内插
#pragma once
#include "Matrix.h"
#include "ImgData_R_W.h"

class SatelliteOrbitFit
{
public:
	SatelliteOrbitFit();
	~SatelliteOrbitFit();

	// 函数声明：
	static void GetOrbitParameters(double * Xs, double * Ys, double * Zs, double & t0, double & dt, double & a, double & b, double & T_azi0, double & T_r0);
	// 三次多项式拟合内插轨道参数：求解系数矩阵Ax/Ay/Az
	static void OrbitInterpolation(double *X, double *T, double *A);
	// 根据系数矩阵计算卫星任意时刻轨道矢量v：
	static void StateliteOrbitVector(double t, double *vectors, double *A);
};

