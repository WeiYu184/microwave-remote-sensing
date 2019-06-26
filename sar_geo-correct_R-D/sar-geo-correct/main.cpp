// Geo-correction Algorithm for SAR Image based on R-D Model
// By weiyu
// main.cpp : 定义控制台应用程序的入口点。
//

#include "ImgData_R_W.h"
#include "SatelliteOrbitFit.h"
#include "R_DCorrect.h"

constexpr auto PAI = 3.1415926535897932384626433832795;

// 函数声明：地心坐标到大地坐标的转换
void Coor_Trans(double a, double b, double B, double L, double H, double &Xp, double &Yp, double &Zp);

int main(void)
{
	/////////////////////////////////////////////////////   一、数据预处理   ////////////////////////////////
	double Xs[5], Ys[5], Zs[5];
	double t0, dt, a, b, T_azi0, T_r0;
	SatelliteOrbitFit::GetOrbitParameters(Xs, Ys, Zs, t0, dt, a, b, T_azi0, T_r0);
	// 参考点的当前时间（设初始时间为0）：
	double t[5] = { 0, dt * 1, dt * 2, dt * 3, dt * 4 };
	// 三次多项式拟合中t的矩阵：
	double T[20] = { 1, t[0], t[0] * t[0], t[0] * t[0] * t[0], 1, t[1], t[1] * t[1], t[1] * t[1] * t[1], 1, t[2], t[2] * t[2], t[2] * t[2] * t[2], 1, t[3], t[3] * t[3], t[3] * t[3] * t[3], 1, t[4], t[4] * t[4], t[4] * t[4] * t[4], };
	// 初始化系数矩阵
	double Ax[4] = { 0 }, Ay[4] = { 0 }, Az[4] = { 0 };

	/////////////////////////////////////////////////////   二、轨道拟合   ////////////////////////////////
	SatelliteOrbitFit::OrbitInterpolation(Xs, T, Ax);
	SatelliteOrbitFit::OrbitInterpolation(Ys, T, Ay);
	SatelliteOrbitFit::OrbitInterpolation(Zs, T, Az);

	//////////////////////////////////////////////////////  三、读取DEM   /////////////////////////////////
	//（dem.img是已在ERDAS中完成裁剪和加密后的ers2dem数据）
	string DEMImg = "dem.img";
	double DEMInform[6] = { 0 };
	int DEM_width, DEM_height, DEM_bands;
	unsigned short *Dem = ImgData_R_W::ImgRead(DEMImg, DEM_height, DEM_width, DEM_bands,DEMInform);
	//得到图像左上角点经纬度L0、B0：
	double L0 = DEMInform[0],B0 = DEMInform[3];

	//////////////////////////////////////////////////////  四、读取SAR   //////////////////////////////////
	string imgSAR = "dat_01.img";
	double SARInform[6] = { 0 };
	//定义图像宽度、长度、波段数
	int SAR_width, SAR_height, SAR_bands;
	unsigned short *Sar = ImgData_R_W::ImgRead(imgSAR, SAR_height, SAR_width, SAR_bands, SARInform);

	/////////////////////////////////////////////// 五、利用R-D方程进行几何校正 ////////////////////////////
	//************** Step1：参数定义及初始化工作
	//定义当前点的经纬度L、B,大地坐标Xp、Yp、Zp,行列号I,J;
	double B, L, Xp, Yp, Zp;
	int I, J;
	// 存储几何校正后的影像:
	unsigned short *correct_img = new unsigned short[DEM_width * DEM_height];
	// 创建数组方位向时间Tazi：
	double *Tazi = new double[DEM_width * DEM_height];
	// Tazi初值取SAR影像中间行的方位向时间:
	for (int i = 0; i < DEM_width; i++)
		for (int j = 0; j < DEM_height; j++)
			Tazi[i * DEM_height + j] = 13000 / 1679.9020000 + T_azi0; 
	cout << "processing...\n";
	//************** Step2：遍历DEM，迭代求解行列号及成像时间
	for (int i = 0; i < DEM_width; i++) {
		for (int j = 0; j < DEM_height; j++) {
			//（1） 求当前点的经纬度坐标：
			// 0.000139是一个像元的经/纬度宽度
			B = B0 - i * 0.000139;
			L = L0 + j * 0.000139;
			// （2）坐标转换：
			Coor_Trans(a, b, B, L, Dem[i * DEM_height + j], Xp, Yp, Zp);
			// （3）计算方位向时间T_azi：
			R_DCorrect::CalculateT_azi(Tazi[i*DEM_height + j], Xp, Yp, Zp, Ax, Ay, Az);
			// （4）最近邻像元法计算在DEM影像上的行列号：
			R_DCorrect::CalculateDEMIJ(I, J, Tazi[i*DEM_height + j], Ax, Ay, Az, t0, Xp, Yp, Zp);
			// （5）DEM生成：行列号在Sar影像范围内赋值，不在范围内的赋0
			if (I < SAR_width && J < SAR_height && I >= 0 && J >= 0)
				correct_img[i*DEM_height + j] = Sar[I * SAR_height + J];
			else  correct_img[i*DEM_height + j] = 0;
		}
		if( i%82 == 0)cout << "processing..."<< i/82<<"%\n";
	}
	cout<<"Success!";
	///////////////////////////////////////////////     六、校正后影像生成      ////////////////////////////
	int band = 1;
	unsigned short **cor_img = new unsigned short *[band];
	cor_img[0] = new unsigned short[DEM_width * DEM_height];
	//将correct_img写入img
	for (int i = 0; i < DEM_width; i++) 
		for (int j = 0; j < DEM_height; j++) 
			cor_img[0][i*DEM_height + j] = correct_img[i*DEM_height + j];
	// 存储校正后影像
	string filename = "geo-correct.tif";
	ImgData_R_W::ImgStore(cor_img, DEM_height, DEM_width,1, filename);
	return 0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////
}

// 函数定义：参心大地坐标转化为参心直角坐标（经纬度高程BLH->直角坐标XYZ）
void Coor_Trans(double a, double b, double B, double L, double H, double &Xp, double &Yp, double &Zp)
{
	//先将经纬度转化为弧度
	B *= PAI / 180;
	L *= PAI / 180;
	double e = sqrt(a * a - b * b) / a;
	//第一辅助系数
	double W = sqrt(1 - e * e * sin(B) * sin(B));
	//椭圆面卯酉圈曲率半径
	double N = a / W;
	Xp = (N + H) * cos(B) * cos(L);
	Yp = (N + H) * cos(B) * sin(L);
	Zp = (N * (1 - e * e) + H) * sin(B);
}