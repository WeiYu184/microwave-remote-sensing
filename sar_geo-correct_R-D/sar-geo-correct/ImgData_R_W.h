#pragma once
#include <iostream>
using namespace std;

class ImgData_R_W
{
public:
	ImgData_R_W();
	~ImgData_R_W();

	// img读写函数声明：
	// 读取img：
	static unsigned short * ImgRead(string imgPathSrc, int &width, int &height, int &bands, double *corner);
	// 存储img：
	static bool ImgStore(unsigned short **pImageData, int width, int height, int bands, string savePath);

};
