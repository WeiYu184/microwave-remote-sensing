#pragma once
#include<stdio.h>
#include<iostream>

using namespace std;

class Matrix
{
public:
	Matrix();
	~Matrix();


	// 矩阵运算函数声明：

	// 矩阵转置：AT = A^T
	static void Mat_ATrans(double *A, double *AT, int row, int col);
	// 矩阵相乘：C = A .* B
	static void Mat_ABMultiply(double *A, double *B, double *C, int row_1, int col_1, int col_2);
	// 生成单位阵I：
	static void Mat_AUnit(double *I, int n);
	// 交换行位置：
	static void Mat_ASwap(double *A, int i, int line, int n);
	// 形成上三角阵：C = uptriangle（A）；
	static void Mat_AUptriangle(double *A, double *C, int n);
	// A转化成单位矩阵：
	static void Mat_A2Unit(double*A, double*C, int n);
	// 矩阵求逆：C = inv（A）
	static void Mat_AInv(double *A, double *C, int n);
};

