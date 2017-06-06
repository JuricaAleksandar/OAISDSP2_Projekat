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