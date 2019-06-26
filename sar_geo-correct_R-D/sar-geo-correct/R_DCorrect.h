// 利用R-D方程进行几何校正
#pragma once
#include "math.h"
#include "memory.h"
#include "Matrix.h"
#include "ImgData_R_W.h"
#include "SatelliteOrbitFit.h"

class R_DCorrect
{
public:
	R_DCorrect();
	~R_DCorrect();

	static void CalculateT_azi(double & T_azi, double Xp, double Yp, double Zp, double * Ax, double * Ay, double * Az);
	static void CalculateDEMIJ(int & Dem_I, int & Dem_J, double T_azi, double * XA, double * YA, double * ZA, double T_start, double Xp, double Yp, double Zp);
};

