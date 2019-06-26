#include "ImgData_R_W.h"
#include <iostream>
#include "gdal_priv.h"
#include "gdal.h"

using namespace std;


ImgData_R_W::ImgData_R_W()
{
}


ImgData_R_W::~ImgData_R_W()
{
}

// img影像读写函数定义：

// img影像读取：


unsigned short * ImgData_R_W::ImgRead(string imgPathSrc, int &width, int &height, int &bands, double *corner)
{
	GDALAllRegister();//注册、读取图像
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	GDALDataset* pInDataset1 = (GDALDataset*)GDALOpen(imgPathSrc.c_str(), GA_ReadOnly);
	if (pInDataset1 == NULL)
	{
		printf("无法打开图像\n");
		GDALDestroyDriverManager();
	}
	//读取图像数据的参数并显示
	width = pInDataset1->GetRasterXSize();
	height = pInDataset1->GetRasterYSize();
	bands = pInDataset1->GetRasterCount();


	//开辟内存
	unsigned short *pImgData = new unsigned short[height * width];

	GDALRasterBand * pInRasterBand1 = pInDataset1->GetRasterBand(1);
	CPLErr error;
	error = pInRasterBand1->RasterIO(GF_Read, 0, 0, width, height, pImgData, width, height, GDT_UInt16, 0, 0);
	if (error == CE_Failure)
	{
		printf("图像读取出错\n");
		GDALDestroyDriverManager();
	}

	pInDataset1->GetGeoTransform(corner);

	//关闭gdal库相关的波段驱动和释放内存，读取图像结束
	GDALClose(pInDataset1);

	return pImgData;
}

// img影像存储：
bool ImgData_R_W::ImgStore(unsigned short **pImageData, int width, int height, int bands, string savePath)
// 图像信息pImageData，宽width，高height，波段数nChannels，保存路径savePath
{
	//注册驱动
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	GDALDriver *pMemDriver = NULL;
	pMemDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if (pMemDriver == NULL)
		return false;

	GDALDataset * pMemDataSet = pMemDriver->Create(savePath.c_str(), width, height, bands, GDT_UInt16, NULL);
	if (pMemDataSet == NULL)
	{
		printf("创建图像出错\n");
		return  false;
	}
	GDALRasterBand *pBand = NULL;

	for (int i = 1; i <= bands; i++)
	{
		pBand = pMemDataSet->GetRasterBand(i);
		CPLErr err = pBand->RasterIO
		(GF_Write, 0, 0, width, height, pImageData[i - 1], width, height, GDT_UInt16, 0, 0);
		if (err == CE_Failure)
		{
			printf("保存图像出错\n");
			return  false;
		}
	}

	//关闭驱动
	GDALClose(pMemDataSet);
	GetGDALDriverManager()->DeregisterDriver(pMemDriver);

	return true;
}