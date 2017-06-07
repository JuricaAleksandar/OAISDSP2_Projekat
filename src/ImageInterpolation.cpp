#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <math.h>


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

			if(ii<ySize - 1)
				ii = (i - 1) / vScalingFactor + 1;

			if (jj<xSize - 1)
				jj = (j - 1) / hScalingFactor + 1;

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
			uchar a = j / hScalingFactor - floor(j / hScalingFactor);
			uchar b = i / vScalingFactor - floor(i / vScalingFactor);

			int ii = (i - 1) / vScalingFactor;
			int jj = (j - 1) / hScalingFactor;

			int iii = ii;
			int jjj = jj;

			if (ii < ySize - 1)
				iii = ii + 1;

			if (jj < xSize - 1)
				jjj = jj + 1;

			NEW_Y_buff[i*newXSize + j] = (1 - a)*(1 - b)*Y_buff[ii*xSize + jj] +
				(1 - a)*b*Y_buff[ii*xSize + jjj] +
				a*(1 - b)*Y_buff[iii*xSize + jj] +
				a*b*Y_buff[iii*xSize + jjj];
		}

	for (int i = 0; i < newYSize / 2; i++)
		for (int j = 0; j < newXSize / 2; j++) {
			uchar a = j / hScalingFactor - floor(j / hScalingFactor);
			uchar b = i / vScalingFactor - floor(i / vScalingFactor);

			int ii = (i - 1) / vScalingFactor;
			int jj = (j - 1) / hScalingFactor;

			int iii = ii;
			int jjj = jj;

			if (ii < ySize/2 - 1)
				iii = ii + 1;

			if (jj < xSize/2 - 1)
				jjj = jj + 1;

			NEW_U_buff[i*newXSize/2 + j] = (1 - a)*(1 - b)*U_buff[ii*xSize/2 + jj] +
				(1 - a)*b*Y_buff[ii*xSize/2 + jjj] +
				a*(1 - b)*Y_buff[iii*xSize/2 + jj] +
				a*b*Y_buff[iii*xSize/2 + jjj];
			NEW_V_buff[i*newXSize / 2 + j] = (1 - a)*(1 - b)*V_buff[ii*xSize / 2 + jj] +
				(1 - a)*b*Y_buff[ii*xSize / 2 + jjj] +
				a*(1 - b)*Y_buff[iii*xSize / 2 + jj] +
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
			cubicIntRes[0] = cubicInterpolateUchar(Y_buff+(ii-1)*xSize + jj - 1, dHorizontal);
			cubicIntRes[1] = cubicInterpolateUchar(Y_buff+ii*xSize + jj - 1, dHorizontal);
			cubicIntRes[2] = cubicInterpolateUchar(Y_buff+(ii+1)*xSize + jj - 1, dHorizontal);
			cubicIntRes[3] = cubicInterpolateUchar(Y_buff+(ii+2)*xSize + jj - 1, dHorizontal);

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

			cubicIntRes[0] = cubicInterpolateChar(U_buff + (ii - 1)*xSize/2 + jj - 1, dHorizontal);
			cubicIntRes[1] = cubicInterpolateChar(U_buff + ii*xSize/2 + jj - 1, dHorizontal);
			cubicIntRes[2] = cubicInterpolateChar(U_buff + (ii + 1)*xSize/2 + jj - 1, dHorizontal);
			cubicIntRes[3] = cubicInterpolateChar(U_buff + (ii + 2)*xSize/2 + jj - 1, dHorizontal);

			NEW_U_buff[i*newXSize/2 + j] = cubicInterpolateChar(cubicIntRes, dVertical);

			cubicIntRes[0] = cubicInterpolateChar(V_buff + (ii - 1)*xSize / 2 + jj - 1, dHorizontal);
			cubicIntRes[1] = cubicInterpolateChar(V_buff + ii*xSize / 2 + jj - 1, dHorizontal);
			cubicIntRes[2] = cubicInterpolateChar(V_buff + (ii + 1)*xSize / 2 + jj - 1, dHorizontal);
			cubicIntRes[3] = cubicInterpolateChar(V_buff + (ii + 2)*xSize / 2 + jj - 1, dHorizontal);

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
	/* TO DO */
}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}