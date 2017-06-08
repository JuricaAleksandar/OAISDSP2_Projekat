#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <math.h>

#define PI 3.14159265

void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* Create buffers for YUV image */
	uchar* Y_buff = new uchar[xSize*ySize];
	char* U_buff = new char[xSize*ySize / 4];
	char* V_buff = new char[xSize*ySize / 4];

	uchar* NEW_Y_buff = new uchar[newXSize*newYSize];
	char* NEW_U_buff = new char[newXSize*newYSize / 4];
	char* NEW_V_buff = new char[newXSize*newYSize / 4];

	double hScalingFactor = (double)newXSize / xSize;
	double vScalingFactor = (double)newYSize / ySize;

	/* Convert input image to YUV420 image */
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	for(int i = 0; i < newYSize; i++)
		for (int j = 0; j < newXSize; j++)
		{
			int ii = (i - 1) / vScalingFactor;
			int jj = (j - 1) / hScalingFactor;

			if(ii < ySize - 1)
				ii++;

			if (jj < xSize - 1)
				jj++;

			NEW_Y_buff[i*newXSize + j] = Y_buff[ii*xSize + jj];
		}

	for (int i = 0; i < newYSize / 2; i++)
		for (int j = 0; j < newXSize / 2; j++)
		{
			int ii = (i - 1) / vScalingFactor;
			int jj = (j - 1) / hScalingFactor;

			if (ii<ySize/2 - 1)
				ii = (i - 1) / vScalingFactor + 1;

			if (jj<xSize - 1)
				jj = (j - 1) / hScalingFactor + 1;

			NEW_V_buff[i*newXSize / 2 + j] = V_buff[ii*xSize / 2 + jj];
			NEW_U_buff[i*newXSize / 2 + j] = U_buff[ii*xSize / 2 + jj];
		}

	/* Convert YUV image back to RGB */
	YUV420toRGB(NEW_Y_buff, NEW_U_buff, NEW_V_buff, newXSize, newYSize, output);

	/* Delete used memory buffers */
	delete[] Y_buff;
	delete[] U_buff;
	delete[] V_buff;

	delete[] NEW_Y_buff;
	delete[] NEW_U_buff;
	delete[] NEW_V_buff;
}

void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* Create buffers for YUV image */
	uchar* Y_buff = new uchar[xSize*ySize];
	char* U_buff = new char[xSize*ySize / 4];
	char* V_buff = new char[xSize*ySize / 4];

	uchar* NEW_Y_buff = new uchar[newXSize*newYSize];
	char* NEW_U_buff = new char[newXSize*newYSize / 4];
	char* NEW_V_buff = new char[newXSize*newYSize / 4];

	double hScalingFactor = (double)newXSize / xSize;
	double vScalingFactor = (double)newYSize / ySize;

	/* Convert input image to YUV420 image */
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	for (int i = 0; i < newYSize; i++)
		for (int j = 0; j < newXSize; j++){
			
			double a = j / hScalingFactor - floor(j / hScalingFactor);
			double b = i / vScalingFactor - floor(i / vScalingFactor);

			int ii = floor((i - 1) / vScalingFactor);
			int jj = floor((j - 1) / hScalingFactor);

			if (ii < ySize - 1)
				ii++;

			if (jj < xSize - 1)
				jj++;

			int iii = ii;
			int jjj = jj;

			if (iii < ySize - 1)
				iii++;

			if (jjj < xSize - 1)
				jjj++;

			uchar gorelevo = Y_buff[ii*xSize + jj];
			uchar goredesno = Y_buff[ii*xSize + jjj];
			uchar dolelevo = Y_buff[iii*xSize + jj];
			uchar doledesno = Y_buff[iii*xSize + jjj];

			NEW_Y_buff[i*newXSize + j] = (1 - a)*(1 - b)*Y_buff[ii*xSize + jj] +
				(1 - a)*b*Y_buff[iii*xSize + jj] +
				a*(1 - b)*Y_buff[ii*xSize + jjj] +
				a*b*Y_buff[iii*xSize + jjj];

			uchar novi = NEW_Y_buff[i*newXSize + j];
		}

	for (int i = 0; i < newYSize / 2; i++)
		for (int j = 0; j < newXSize / 2; j++) {
			uchar a = j / hScalingFactor - floor(j / hScalingFactor);
			uchar b = i / vScalingFactor - floor(i / vScalingFactor);

			int ii = floor((i - 1) / vScalingFactor);
			int jj = floor((j - 1) / hScalingFactor);

			if (ii < ySize/2 - 1)
				ii++;

			if (jj < xSize/2 - 1)
				jj++;

			int iii = ii;
			int jjj = jj;

			if (iii < ySize / 2 - 1)
				iii++;

			if (jjj < xSize / 2 - 1)
				jjj++;

			NEW_U_buff[i*newXSize/2 + j] = (1 - a)*(1 - b)*U_buff[ii*xSize/2 + jj] +
				(1 - a)*b*Y_buff[iii*xSize/2 + jj] +
				a*(1 - b)*Y_buff[ii*xSize/2 + jjj] +
				a*b*Y_buff[iii*xSize/2 + jjj];
			NEW_V_buff[i*newXSize / 2 + j] = (1 - a)*(1 - b)*V_buff[ii*xSize / 2 + jj] +
				(1 - a)*b*Y_buff[iii*xSize / 2 + jj] +
				a*(1 - b)*Y_buff[ii*xSize / 2 + jjj] +
				a*b*Y_buff[iii*xSize / 2 + jjj];
		}

	/* Convert YUV image back to RGB */
	YUV420toRGB(NEW_Y_buff, NEW_U_buff, NEW_V_buff, newXSize, newYSize, output);

	/* Delete used memory buffers */
	delete[] Y_buff;
	delete[] U_buff;
	delete[] V_buff;

	delete[] NEW_Y_buff;
	delete[] NEW_U_buff;
	delete[] NEW_V_buff;
}

double weight(double d)
{
	if (abs(d) < 1)
		return 3.0*abs(d*d*d) / 2.0 - 5.0*abs(d*d) / 2.0 + 1;
	else if (abs(d) >= 1 && abs(d) < 2)
		return -1.0*abs(d*d*d) / 2.0 + 5.0*abs(d*d) / 2.0 - 4.0*abs(d) + 2;
	else
		return 0;
}

uchar cubicInterpolateUchar(const uchar * input,double d)
{
	int value = ((input[0] * weight(d + 1) + input[1] * weight(d) + input[2] * weight(1 - d) + input[3] * weight(2 - d)));
	if (value > 255)
		return 255;
	else if (value < 0)
		return 0;
	return (uchar)value;
}

char cubicInterpolateChar(const char * input, double d)
{
	int value = ((input[0] * weight(d + 1) + input[1] * weight(d) + input[2] * weight(1 - d) + input[3] * weight(2 - d)));
	if (value > 127)
		return 127;
	else if (value < -128)
		return -128;
	return (char)value;
}

void bicubicInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* Create buffers for YUV image */
	uchar* Y_buff = new uchar[xSize*ySize];
	char* U_buff = new char[xSize*ySize / 4];
	char* V_buff = new char[xSize*ySize / 4];

	uchar* NEW_Y_buff = new uchar[newXSize*newYSize];
	char* NEW_U_buff = new char[newXSize*newYSize / 4];
	char* NEW_V_buff = new char[newXSize*newYSize / 4];

	double hScalingFactor = (double)newXSize / xSize;
	double vScalingFactor = (double)newYSize / ySize;

	/* Convert input image to YUV420 image */
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	for(int i = 0; i < newYSize; i++)
		for (int j = 0; j < newXSize; j++)
		{
			double dVertical = i / vScalingFactor - floor(i / vScalingFactor);
			double dHorizontal = j / hScalingFactor - floor(j / hScalingFactor);

			int ii = i / vScalingFactor;
			int jj = j / hScalingFactor;

			uchar cubicIntRes[4];
			for (int m = 0; m < 4; m++)
			{
				int dif;
				if(dif = ii - 1 + m - ySize > 0)
					cubicIntRes[m] = cubicInterpolateUchar(Y_buff + (ii - 1 + m - dif)*xSize + jj - 1, dHorizontal);
				else
					cubicIntRes[m] = cubicInterpolateUchar(Y_buff + (ii - 1 + m)*xSize + jj - 1, dHorizontal);
			}

			NEW_Y_buff[i*newXSize + j] = cubicInterpolateUchar(cubicIntRes, dVertical);
		}

	for (int i = 0; i < newYSize/2; i++)
		for (int j = 0; j < newXSize/2; j++)
		{
			double dVertical = i / vScalingFactor - floor(i / vScalingFactor);
			double dHorizontal = j / hScalingFactor - floor(j / hScalingFactor);

			int ii = i / vScalingFactor;
			int jj = j / hScalingFactor;

			char cubicIntRes[4];

			for (int m = 0; m < 4; m++)
			{
				int dif;
				if (dif = ii - 1 + m - ySize/2 > 0)
					cubicIntRes[m] = cubicInterpolateChar(U_buff + (ii - 1 + m - dif)*xSize/2 + jj - 1, dHorizontal);
				else
					cubicIntRes[m] = cubicInterpolateChar(U_buff + (ii - 1 + m)*xSize/2 + jj - 1, dHorizontal);
			}

			NEW_U_buff[i*newXSize / 2 + j] = cubicInterpolateChar(cubicIntRes, dVertical);

			for (int m = 0; m < 4; m++)
			{
				int dif;
				if (dif = ii - 1 + m - ySize / 2 > 0)
					cubicIntRes[m] = cubicInterpolateChar(V_buff + (ii - 1 + m - dif)*xSize / 2 + jj - 1, dHorizontal);
				else
					cubicIntRes[m] = cubicInterpolateChar(V_buff + (ii - 1 + m)*xSize / 2 + jj - 1, dHorizontal);
			}

			NEW_V_buff[i*newXSize / 2 + j] = cubicInterpolateChar(cubicIntRes, dVertical);
		}
	/* Convert YUV image back to RGB */
	YUV420toRGB(NEW_Y_buff, NEW_U_buff, NEW_V_buff, newXSize, newYSize, output);

	/* Delete used memory buffers */
	delete[] Y_buff;
	delete[] U_buff;
	delete[] V_buff;

	delete[] NEW_Y_buff;
	delete[] NEW_U_buff;
	delete[] NEW_V_buff;
}

void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* Create buffers for YUV image */
	uchar* Y_buff = new uchar[xSize*ySize];
	char* U_buff = new char[xSize*ySize / 4];
	char* V_buff = new char[xSize*ySize / 4];

	uchar* NEW_Y_buff = new uchar[xSize*ySize];
	char* NEW_U_buff = new char[xSize*ySize / 4];
	char* NEW_V_buff = new char[xSize*ySize / 4];

	/* Convert input image to YUV420 image */
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	for (int i = 0; i < ySize; i++)
		for (int j = 0; j < xSize; j++)
		{
			int ii = i*cos(PI*angle/180) + j*sin(PI*angle/180) - m*sin(PI*angle/180) - n*cos(PI*angle/180) + n;
			int jj = j*cos(PI*angle/180) - i*sin(PI*angle/180) - m*cos(PI*angle/180) + n*sin(PI*angle/180) + m;

			if (ii < 0 || ii > ySize || jj < 0 || jj > xSize)
				NEW_Y_buff[i*xSize + j] = 0;
			else
				NEW_Y_buff[i*xSize + j] = Y_buff[ii*xSize + jj];
		}

	for (int i = 0; i < ySize/2; i++)
		for (int j = 0; j < xSize/2; j++)
		{
			int ii = i*cos(PI*angle / 180) + j*sin(PI*angle / 180) - m*sin(PI*angle / 180)/2 - n*cos(PI*angle / 180)/2 + n/2;
			int jj = j*cos(PI*angle / 180) - i*sin(PI*angle / 180) - m*cos(PI*angle / 180)/2 + n*sin(PI*angle / 180)/2 + m/2;

			if (ii < 0 || ii > ySize / 2 || jj < 0 || jj > xSize / 2)
			{
				NEW_U_buff[i*xSize / 2 + j] = 0;
				NEW_V_buff[i*xSize / 2 + j] = 0;
			}
			else
			{
				NEW_U_buff[i*xSize / 2 + j] = U_buff[ii*xSize / 2 + jj];
				NEW_V_buff[i*xSize / 2 + j] = V_buff[ii*xSize / 2 + jj];
			}
		}

	/* Convert YUV image back to RGB */
	YUV420toRGB(NEW_Y_buff, NEW_U_buff, NEW_V_buff, xSize, ySize, output);

	/* Delete used memory buffers */
	delete[] Y_buff;
	delete[] U_buff;
	delete[] V_buff;

	delete[] NEW_Y_buff;
	delete[] NEW_U_buff;
	delete[] NEW_V_buff;
}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* Create buffers for YUV image */
	uchar* Y_buff = new uchar[xSize*ySize];
	char* U_buff = new char[xSize*ySize / 4];
	char* V_buff = new char[xSize*ySize / 4];

	uchar* NEW_Y_buff = new uchar[xSize*ySize];
	char* NEW_U_buff = new char[xSize*ySize / 4];
	char* NEW_V_buff = new char[xSize*ySize / 4];

	/* Convert input image to YUV420 image */
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	for (int i = 0; i < ySize; i++)
		for (int j = 0; j < xSize; j++)
		{
			int ii = i*cos(PI*angle / 180) + j*sin(PI*angle / 180) - m*sin(PI*angle / 180) - n*cos(PI*angle / 180) + n;
			int jj = j*cos(PI*angle / 180) - i*sin(PI*angle / 180) - m*cos(PI*angle / 180) + n*sin(PI*angle / 180) + m;

			if (ii < 0 || ii > ySize || jj < 0 || jj > xSize)
				NEW_Y_buff[i*xSize + j] = 0;
			else
				NEW_Y_buff[i*xSize + j] = Y_buff[ii*xSize + jj];
		}

	for (int i = 0; i < ySize / 2; i++)
		for (int j = 0; j < xSize / 2; j++)
		{
			int ii = i*cos(PI*angle / 180) + j*sin(PI*angle / 180) - m*sin(PI*angle / 180) / 2 - n*cos(PI*angle / 180) / 2 + n / 2;
			int jj = j*cos(PI*angle / 180) - i*sin(PI*angle / 180) - m*cos(PI*angle / 180) / 2 + n*sin(PI*angle / 180) / 2 + m / 2;

			if (ii < 0 || ii > ySize / 2 || jj < 0 || jj > xSize / 2)
			{
				NEW_U_buff[i*xSize / 2 + j] = 0;
				NEW_V_buff[i*xSize / 2 + j] = 0;
			}
			else
			{
				NEW_U_buff[i*xSize / 2 + j] = U_buff[ii*xSize / 2 + jj];
				NEW_V_buff[i*xSize / 2 + j] = V_buff[ii*xSize / 2 + jj];
			}
		}

	/* Convert YUV image back to RGB */
	YUV420toRGB(NEW_Y_buff, NEW_U_buff, NEW_V_buff, xSize, ySize, output);

	/* Delete used memory buffers */
	delete[] Y_buff;
	delete[] U_buff;
	delete[] V_buff;

	delete[] NEW_Y_buff;
	delete[] NEW_U_buff;
	delete[] NEW_V_buff;
}