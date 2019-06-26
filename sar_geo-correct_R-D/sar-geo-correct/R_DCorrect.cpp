// 利用R-D方程进行几何校正
#include "R_DCorrect.h"



R_DCorrect::R_DCorrect()
{
}


R_DCorrect::~R_DCorrect()
{
}

// 利用R-D模型迭代求解方位向时间、距离向时间、进而求出对应SAR影像上行列号i，j
void R_DCorrect::CalculateT_azi(double & T_azi, double Xp, double Yp, double Zp, double *Ax, double *Ay, double *Az)
{
	// T_azi是方位向时间初始值；Xp、Yp、Zp是DEM上的点对应大地坐标,A是算出来的系数矩阵
	// 迭代次数：
	double delta_T = 10;
	int counter = 0;//计数器
	while (counter <= 10 && abs(delta_T) > 1.0E-12)
	{
		//分别存放卫星X、Y、Z维度上的坐标、速度、加速度
		double vecX[3];
		double vecY[3];
		double vecZ[3];
		// 计算卫星轨道矢量
		SatelliteOrbitFit::StateliteOrbitVector(T_azi, vecX, Ax);
		SatelliteOrbitFit::StateliteOrbitVector(T_azi, vecY, Ay);
		SatelliteOrbitFit::StateliteOrbitVector(T_azi, vecZ, Az);
		//算出地面点与卫星相对位置矢量Rsp
		double Rsp[3] = { vecX[0] - Xp ,vecY[0] - Yp, vecZ[0] - Zp };
		//计算多普勒方程
		double Ft = Rsp[0] * vecX[1] + Rsp[1] * vecY[1] + Rsp[2] * vecZ[1];
		//对Ft关于t求偏微分
		double ft = (vecX[1] * vecX[1] + vecY[1] * vecY[1] + vecZ[1] * vecZ[1]) + (Rsp[0] * vecX[2] + Rsp[1] * vecY[2] + Rsp[2] * vecZ[2]);
		//计算方位向时间改正数delta_t
		delta_T = -(Ft / ft);

		//方位向时间经过改正后的值
		T_azi += delta_T;
		//计数加一
		counter++;
	}//while

	//判断循环结束的终止原因
	if (counter > 10)printf("迭代次数超限\n");
	if (delta_T <= 1.0E-12)
	{/*迭代终止*/
	}
	else printf("迭代异常终止\n");
}

//最邻近像元法求解当前点在DEM图像上的行列号Dem_I、Dem_J
//参数说明：分别：当前点在Dem上的行号、列号、当前点方位向时间
void R_DCorrect::CalculateDEMIJ(int &Dem_I, int &Dem_J, double T_azi, double * XA, double * YA, double * ZA, double T_start, double Xp, double Yp, double Zp)
{
	double C = 299792458;
	// 卫星坐标：
	double  S_XYZ[3] = { 0 };
	//多项式拟合法计算T_azi时刻的卫星坐标
	S_XYZ[0] = XA[0] + XA[1] * (T_azi - T_start) + XA[2] * (T_azi - T_start) * (T_azi - T_start) + XA[3] * (T_azi - T_start) * (T_azi - T_start) * (T_azi - T_start);
	S_XYZ[1] = YA[0] + YA[1] * (T_azi - T_start) + YA[2] * (T_azi - T_start) * (T_azi - T_start) + YA[3] * (T_azi - T_start) * (T_azi - T_start) * (T_azi - T_start);
	S_XYZ[2] = ZA[0] + ZA[1] * (T_azi - T_start) + ZA[2] * (T_azi - T_start) * (T_azi - T_start) + ZA[3] * (T_azi - T_start) * (T_azi - T_start) * (T_azi - T_start);
	//定义卫星与地面点的斜距
	double R = sqrt(pow(S_XYZ[0] - Xp, 2) + pow(S_XYZ[1] - Yp, 2) + pow(S_XYZ[2] - Zp, 2));
	//计算距离向时间
	double T_r = (2 * R) / C;

	//下面计算当前像元在Dem影像上的行列号
	//首先得到一些已有的数据
	//方位向时间初始值
	double T_azi0 = 11728.060;
	//距离向时间初始值
	double T_r0 = 5.5365520 / 1000;
	//方位向时间间隔
	double dT_azi = 1.0 / 1679.9020000;
	//距离向时间间隔
	double dT_r = 1.0 / (18.9624680*1E+06);

	//计算t_azi时刻对应的SAR影像行列号，采用最近邻元法
	Dem_I = (T_azi - T_azi0) / dT_azi + 0.5;//行号
	Dem_J = (T_r - T_r0) / dT_r + 0.5; //列号
}
