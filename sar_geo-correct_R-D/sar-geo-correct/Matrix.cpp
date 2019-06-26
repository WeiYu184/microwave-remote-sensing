#include "Matrix.h"

Matrix::Matrix()
{
}


Matrix::~Matrix()
{
}


// 矩阵运算函数定义：

// 矩阵转置：AT = A^T
void Matrix::Mat_ATrans(double *A, double *AT, int row, int col)
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			AT[j*row + i] = A[i*col + j];
		}
	}

}

// 矩阵相乘：C = A .* B
void Matrix::Mat_ABMultiply(double *A, double *B, double *C, int row_1, int col_1, int col_2)
{
	double ret;
	for (int i = 0; i < row_1; i++)
	{
		for (int j = 0; j < col_2; j++)
		{
			ret = 0;
			for (int k = 0; k < col_1; k++)
			{
				ret += A[i * col_1 + k] * B[k * col_2 + j];  
			}
			C[i * col_2 + j] = ret;
		}
	}
}

// 生成单位阵I：
void Matrix::Mat_AUnit(double *I, int n)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j)
			{
				I[i*n + j] = 1;
			}
			else
				I[i*n + j] = 0;
		}
	}
}

// 交换行位置：
void Matrix::Mat_ASwap(double *A, int i, int line, int n)
{
	int j;
	double temp;
	for (j = 0; j < n; j++)
	{
		temp = A[i*n + j];
		A[i*n + j] = A[line*n + j];
		A[line*n + j] = temp;
	}
}

// 形成上三角阵：C = uptriangle（A）
void Matrix::Mat_AUptriangle(double *A, double *C, int n)
{
	int i, j, k, m, line;
	double max, temp, mmul;
	for (i = 0; i < n; i++)
	{
		max = fabs(A[i*n + i]);
		temp = A[i*n + i];
		line = i;
		for (j = i + 1; j < n; j++)
		{
			if (fabs(A[j*n + i]) > max)
			{
				max = fabs(A[j*n + i]);
				temp = A[j*n + i];
				line = j;
			}
		}
		if (max <= 1e-10)
		{
			printf("no inverse array\n");
			return;
		}
		if (line != i)
		{
			Mat_ASwap(A, i, line, n);
			Mat_ASwap(C, i, line, n);
		}
		for (k = 0; k < n; k++)
		{
			A[i*n + k] /= temp;
			C[i*n + k] /= temp;
		}
		for (k = i + 1; k < n; k++)
		{
			mmul = A[k*n + i];
			for (m = 0; m < n; m++)
			{
				A[k*n + m] -= A[i*n + m] * mmul;
				C[k*n + m] -= C[i*n + m] * mmul;
			}
		}
	}

}

// A转化成单位矩阵：
void Matrix::Mat_A2Unit(double*A, double*C, int n)
{
	int i, j, k;
	double mmul;
	for (i = n - 1; i > 0; i--)
	{
		//从下往上每一行进行计算
		for (j = i - 1; j >= 0; j--)
		{
			mmul = A[j*n + i];
			A[j*n + i] -= A[i*n + i] * mmul;
			for (k = 0; k < n; k++)
			{
				//每一行都减去行末元素*i行对应列的值
				C[j*n + k] -= C[i*n + k] * mmul;
			}
		}
	}

}

// 矩阵求逆：C = inv（A） 
void Matrix::Mat_AInv(double *A, double *C, int n)
{
	//将C初始化为单位阵
	Mat_AUnit(C, n);
	//将A化为了上三角阵
	Mat_AUptriangle(A, C, n);
	//将A化为了单位阵，这时C即为所求
	Mat_A2Unit(A, C, n);
}
