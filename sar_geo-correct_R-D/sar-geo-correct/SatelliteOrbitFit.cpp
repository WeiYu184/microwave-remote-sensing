// 类SatelliteOrbitFit：通过三次多项式拟合法实现卫星轨道参数内插

#include "SatelliteOrbitFit.h"


SatelliteOrbitFit::SatelliteOrbitFit()
{
}


SatelliteOrbitFit::~SatelliteOrbitFit()
{
}

// 函数定义：
// 轨道参数读取：
/*void SatelliteOrbitFit::GetOrbitParameters(double *X, double *T, double *Ax)

}*/

void SatelliteOrbitFit::GetOrbitParameters(double *Xs, double *Ys, double *Zs, double &t0, double &dt, double &a,double &b, double &T_azi0, double &T_r0)
{
	//////首先将需要的参数定义并赋值/////////////////////
//五个参考点坐标X,Y,Z（单位：m）
	FILE *fp;
	fopen_s(&fp,"parameter.txt", "r");
	if (fp == NULL)
	{
		printf("can't open file!\n");
	}
	for( int i = 0; i < 5; i++)
		fscanf_s(fp, "%lf %lf %lf", &Xs[i], &Ys[i], &Zs[i]);
	
	//第一成像时间，时间间隔（单位：s）
	t0 = 1.172766900000000E+04;
	dt = 4.167000000000000E+00;
	//椭球参数（单位：公里）
	a = 6378144;//长半轴
	b = 6356759;//短半轴

	//方位向时间初始值
	T_azi0 = 11728.060;
	//距离向时间初始值
	T_r0 = 5.5365520 / 1000;
	return;
}

// 三次多项式拟合内插轨道参数：求解系数矩阵Ax/Ay/Az
void SatelliteOrbitFit::OrbitInterpolation(double *X, double *T, double *Ax)
// 位置坐标X[5*1]，当前时刻T[5*4]，系数矩阵Ax[4*1]
{

	// T_T = trans(T)
	double T_T[20] = { 0 };
	Matrix::Mat_ATrans(T, T_T, 5, 4);

	// TT(4*4) = T .* T_T
	double TT[4 * 4] = { 0 };
	Matrix::Mat_ABMultiply(T_T, T, TT, 4, 5, 4);

	// TT_inv = inv(TT)
	double TT_inv[4 * 4] = { 0 };
	Matrix::Mat_AInv(TT, TT_inv, 4);

	// tmp = TT_inv .* T_T
	double tmp[4 * 5] = { 0 };
	Matrix::Mat_ABMultiply(TT_inv, T_T, tmp, 4, 4, 5);

	// Ax = tmp .* X
	Matrix::Mat_ABMultiply(tmp, X, Ax, 4, 5, 1);

}

// 根据系数矩阵计算卫星任意时刻轨道矢量v：
void SatelliteOrbitFit::StateliteOrbitVector(double t, double *v, double *Ax)
// 当前时刻t； 卫星轨道向量v； 系数矩阵Ax
{
	double tt = t - 1.172766900000000E+04;//以初始时间为0
	// 计算当前时刻坐标
	v[0] = Ax[0] + Ax[1] * tt + Ax[2] * tt * tt + Ax[3] * tt * tt * tt;
	// 计算当前时刻速度
	v[1] = Ax[1] + 2 * Ax[2] * tt + 3 * Ax[3] * tt * tt;
	// 计算当前时刻加速度
	v[2] = 2 * Ax[2] + 6 * Ax[3] * tt; 

}