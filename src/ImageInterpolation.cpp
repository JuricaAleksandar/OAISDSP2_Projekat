#define _USE_MATH_DEFINES
#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <thread>
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

			int ii = i / vScalingFactor;
			int jj = j / hScalingFactor;

			int iii = ii;
			int jjj = jj;

			if (iii < ySize - 1)
				iii++;

			if (jjj < xSize - 1)
				jjj++;

			NEW_Y_buff[i*newXSize + j] = (1 - a)*(1 - b)*Y_buff[ii*xSize + jj] +
				(1 - a)*b*Y_buff[iii*xSize + jj] +
				a*(1 - b)*Y_buff[ii*xSize + jjj] +
				a*b*Y_buff[iii*xSize + jjj];
		}

	for (int i = 0; i < newYSize / 2; i++)
		for (int j = 0; j < newXSize / 2; j++) {
			uchar a = j / hScalingFactor - floor(j / hScalingFactor);
			uchar b = i / vScalingFactor - floor(i / vScalingFactor);

			int ii = i / vScalingFactor;
			int jj = j / hScalingFactor;

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
	int value = ((input[0] * weight(d) + input[1] * weight(d-1) + input[2] * weight(2 - d) + input[3] * weight(3 - d)));
	if (value > 255)
		return 255;
	else if (value < 0)
		return 0;
	return (uchar)value;
}

char cubicInterpolateChar(const char * input, double d)
{
	int value = ((input[0] * weight(d) + input[1] * weight(d-1) + input[2] * weight(2-d) + input[3] * weight(3 - d)));
	if (value > 127)
		return 127;
	else if (value < -128)
		return -128;
	return (char)value;
}

void processingUV(int xSize,int ySize,int newXSize,int newYSize,double vScalingFactor,double hScalingFactor,char * input,char * output)
{
	for (int i = 0; i < newYSize; i++)
		for (int j = 0; j < newXSize; j++)
		{
			double dVertical = i / vScalingFactor - floor(i / vScalingFactor);
			double dHorizontal = j / hScalingFactor - floor(j / hScalingFactor);

			int ii = i / vScalingFactor;
			int jj = j / hScalingFactor;

			char funInput[4];
			char cubicIntRes[4];

			for (int m = 0; m < 4; m++)
			{
				int tmp = m;
				if (ii + tmp == 0)
					tmp++;
				else if (ii + tmp > ySize - 1)
					tmp = ySize - ii;

				if (jj == 0)
				{
					funInput[0] = input[(ii - 1 + tmp)*xSize];
					funInput[1] = input[(ii - 1 + tmp)*xSize];
					funInput[2] = input[(ii - 1 + tmp)*xSize + 1];
					funInput[3] = input[(ii - 1 + tmp)*xSize + 2];
				}
				else if (jj == xSize - 2)
				{
					funInput[0] = input[(ii - 1 + tmp)*xSize + jj - 1];
					funInput[1] = input[(ii - 1 + tmp)*xSize + jj];
					funInput[2] = input[(ii - 1 + tmp)*xSize + jj + 1];
					funInput[3] = input[(ii - 1 + tmp)*xSize + jj + 1];
				}
				else if (jj == xSize - 1)
				{
					funInput[0] = input[(ii - 1 + tmp)*xSize + jj - 1];
					funInput[1] = input[(ii - 1 + tmp)*xSize + jj];
					funInput[2] = input[(ii - 1 + tmp)*xSize + jj];
					funInput[3] = input[(ii - 1 + tmp)*xSize + jj];
				}
				else
				{
					funInput[0] = input[(ii - 1 + tmp)*xSize + jj - 1];
					funInput[1] = input[(ii - 1 + tmp)*xSize + jj];
					funInput[2] = input[(ii - 1 + tmp)*xSize + jj + 1];
					funInput[3] = input[(ii - 1 + tmp)*xSize + jj + 2];
				}
				cubicIntRes[m] = cubicInterpolateChar(funInput, dHorizontal + 1);
			}
			output[i*newXSize + j] = cubicInterpolateChar(cubicIntRes, dVertical + 1);
		}
}

void processingY(int xSize, int ySize, int newXSize, int newYSize, double vScalingFactor, double hScalingFactor, uchar * input, uchar * output)
{
	for (int i = newYSize/2; i < newYSize; i++)
		for (int j = 0; j < newXSize; j++)
		{
			double dVertical = i / vScalingFactor - floor(i / vScalingFactor);
			double dHorizontal = j / hScalingFactor - floor(j / hScalingFactor);

			int ii = i / vScalingFactor;
			int jj = j / hScalingFactor;

			uchar funInput[4];
			uchar cubicIntRes[4];

			for (int m = 0; m < 4; m++)
			{
				int tmp = m;
				if (ii + tmp == 0)
					tmp++;
				else if (ii + tmp > ySize - 1)
					tmp = ySize - ii;

				if (jj == 0)
				{
					funInput[0] = input[(ii - 1 + tmp)*xSize];
					funInput[1] = input[(ii - 1 + tmp)*xSize];
					funInput[2] = input[(ii - 1 + tmp)*xSize + 1];
					funInput[3] = input[(ii - 1 + tmp)*xSize + 2];
				}
				else if (jj == xSize - 2)
				{
					funInput[0] = input[(ii - 1 + tmp)*xSize + jj - 1];
					funInput[1] = input[(ii - 1 + tmp)*xSize + jj];
					funInput[2] = input[(ii - 1 + tmp)*xSize + jj + 1];
					funInput[3] = input[(ii - 1 + tmp)*xSize + jj + 1];
				}
				else if (jj == xSize - 1)
				{
					funInput[0] = input[(ii - 1 + tmp)*xSize + jj - 1];
					funInput[1] = input[(ii - 1 + tmp)*xSize + jj];
					funInput[2] = input[(ii - 1 + tmp)*xSize + jj];
					funInput[3] = input[(ii - 1 + tmp)*xSize + jj];
				}
				else
				{
					funInput[0] = input[(ii - 1 + tmp)*xSize + jj - 1];
					funInput[1] = input[(ii - 1 + tmp)*xSize + jj];
					funInput[2] = input[(ii - 1 + tmp)*xSize + jj + 1];
					funInput[3] = input[(ii - 1 + tmp)*xSize + jj + 2];
				}
				cubicIntRes[m] = cubicInterpolateUchar(funInput, dHorizontal + 1);
			}
			output[i*newXSize + j] = cubicInterpolateUchar(cubicIntRes, dVertical + 1);
		}
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

	/* U,V and first half of Y buffer processing */
	std::thread uProcessingThread(processingUV, xSize/2,ySize/2,newXSize/2,newYSize/2,vScalingFactor,hScalingFactor,U_buff,NEW_U_buff);
	std::thread vProcessingThread(processingUV, xSize/2, ySize/2, newXSize/2, newYSize/2, vScalingFactor, hScalingFactor, V_buff, NEW_V_buff);
	std::thread halfYProcessingThread(processingY, xSize, ySize, newXSize, newYSize, vScalingFactor, hScalingFactor, Y_buff, NEW_Y_buff);
	
	/* Second half of Y buffer processing */
	for (int i = 0; i < newYSize /2; i++)
		for (int j = 0; j < newXSize; j++)
		{
			double dVertical = i / vScalingFactor - floor(i / vScalingFactor);
			double dHorizontal = j / hScalingFactor - floor(j / hScalingFactor);

			int ii = i / vScalingFactor;
			int jj = j / hScalingFactor;

			uchar funInput[4];
			uchar cubicIntRes[4];
			for (int m = 0; m < 4; m++)
			{
				int tmp = m;
				if (ii + tmp == 0)
					tmp++;
				else if (ii + tmp > ySize - 1)
					tmp = ySize - ii;

				if (jj == 0)
				{
					funInput[0] = Y_buff[(ii - 1 + tmp)*xSize];
					funInput[1] = Y_buff[(ii - 1 + tmp)*xSize];
					funInput[2] = Y_buff[(ii - 1 + tmp)*xSize + 1];
					funInput[3] = Y_buff[(ii - 1 + tmp)*xSize + 2];
				}
				else if (jj == xSize - 2)
				{
					funInput[0] = Y_buff[(ii - 1 + tmp)*xSize + jj - 1];
					funInput[1] = Y_buff[(ii - 1 + tmp)*xSize + jj];
					funInput[2] = Y_buff[(ii - 1 + tmp)*xSize + jj + 1];
					funInput[3] = Y_buff[(ii - 1 + tmp)*xSize + jj + 1];
				}
				else if (jj == xSize - 1)
				{
					funInput[0] = Y_buff[(ii - 1 + tmp)*xSize + jj - 1];
					funInput[1] = Y_buff[(ii - 1 + tmp)*xSize + jj];
					funInput[2] = Y_buff[(ii - 1 + tmp)*xSize + jj];
					funInput[3] = Y_buff[(ii - 1 + tmp)*xSize + jj];
				}
				else
				{
					funInput[0] = Y_buff[(ii - 1 + tmp)*xSize + jj - 1];
					funInput[1] = Y_buff[(ii - 1 + tmp)*xSize + jj];
					funInput[2] = Y_buff[(ii - 1 + tmp)*xSize + jj + 1];
					funInput[3] = Y_buff[(ii - 1 + tmp)*xSize + jj + 2];
				}
				cubicIntRes[m] = cubicInterpolateUchar(funInput, dHorizontal + 1);
			}
			NEW_Y_buff[i*newXSize + j] = cubicInterpolateUchar(cubicIntRes, dVertical + 1);
		}

	uProcessingThread.join();
	vProcessingThread.join();
	halfYProcessingThread.join();

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

	double angleRad = M_PI*angle / 180;

	for (int i = 0; i < ySize; i++)
		for (int j = 0; j < xSize; j++)
		{
			int ii = round(i*cos(angleRad) + j*sin(angleRad) - m*sin(angleRad) - n*cos(angleRad) + n);
			int jj = round(j*cos(angleRad) - i*sin(angleRad) - m*cos(angleRad) + n*sin(angleRad) + m);

			if (ii == ySize)
				ii--;

			if (jj == xSize)
				jj--;

			if (ii < 0 || ii >= ySize || jj < 0 || jj >= xSize)
				NEW_Y_buff[i*xSize + j] = 0;
			else
				NEW_Y_buff[i*xSize + j] = Y_buff[ii*xSize + jj];
		}

	for (int i = 0; i < ySize/2; i++)
		for (int j = 0; j < xSize/2; j++)
		{
			int ii = round(i*cos(angleRad) + j*sin(angleRad) - m*sin(angleRad)/2 - n*cos(angleRad)/2 + n/2);
			int jj = round(j*cos(angleRad) - i*sin(angleRad) - m*cos(angleRad)/2 + n*sin(angleRad)/2 + m/2);

			if (ii == ySize/2)
				ii--;

			if (jj == xSize/2)
				jj--;

			if (ii < 0 || ii >= ySize / 2 || jj < 0 || jj >= xSize / 2)
			{
				NEW_U_buff[i*xSize / 2 + j] = 0;
				NEW_V_buff[i*xSize / 2 + j] = 0;
			}
			else
			{
				NEW_U_buff[i*xSize / 2 + j] = U_buff[ii*xSize/2 + jj];
				NEW_V_buff[i*xSize / 2 + j] = V_buff[ii*xSize/2 + jj];
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

	double angleRad = M_PI*angle / 180;

	for (int i = 0; i < ySize; i++)
		for (int j = 0; j < xSize; j++)
		{
			double a = j*cos(angleRad) - i*sin(angleRad) - m*cos(angleRad) + n*sin(angleRad) + m -
					floor(j*cos(angleRad) - i*sin(angleRad) - m*cos(angleRad) + n*sin(angleRad) + m);

			double b = i*cos(angleRad) + j*sin(angleRad) - m*sin(angleRad) - n*cos(angleRad) + n -
				floor(i*cos(angleRad) + j*sin(angleRad) - m*sin(angleRad) - n*cos(angleRad) + n);

			int ii = i*cos(angleRad) + j*sin(angleRad) - m*sin(angleRad) - n*cos(angleRad) + n;
			int jj = j*cos(angleRad) - i*sin(angleRad) - m*cos(angleRad) + n*sin(angleRad) + m;

			if (ii == ySize)
				ii--;

			if (jj == xSize)
				jj--;

			int iii = ii;
			int jjj = jj;

			if (iii < ySize - 1)
				iii++;

			if (jjj < xSize - 1)
				jjj++;

			if (ii < 0 || ii >= ySize || jj < 0 || jj >= xSize)
				NEW_Y_buff[i*xSize + j] = 0;
			else
				NEW_Y_buff[i*xSize + j] = (1 - a)*(1 - b)*Y_buff[ii*xSize + jj] +
				(1 - a)*b*Y_buff[iii*xSize + jj] +
				a*(1 - b)*Y_buff[ii*xSize + jjj] +
				a*b*Y_buff[iii*xSize + jjj];
		}

	for (int i = 0; i < ySize / 2; i++)
		for (int j = 0; j < xSize / 2; j++)
		{
			double a = j*cos(angleRad) - i*sin(angleRad) - m*cos(angleRad)/2.0 + n*sin(angleRad)/2.0 + m/2.0 -
				floor(j*cos(angleRad) - i*sin(angleRad) - m*cos(angleRad) / 2.0 + n*sin(angleRad) / 2.0 + m / 2.0);

			double b = i*cos(angleRad) + j*sin(angleRad) - m*sin(angleRad)/2.0 - n*cos(angleRad)/2.0 + n/2.0 -
				floor(i*cos(angleRad) + j*sin(angleRad) - m*sin(angleRad) / 2.0 - n*cos(angleRad) / 2.0 + n / 2.0);

			int ii = i*cos(angleRad) + j*sin(angleRad) - m*sin(angleRad) / 2 - n*cos(angleRad) / 2 + n / 2;
			int jj = j*cos(angleRad) - i*sin(angleRad) - m*cos(angleRad) / 2 + n*sin(angleRad) / 2 + m / 2;

			if (ii == ySize/2)
				ii--;

			if (jj == xSize/2)
				jj--;

			int iii = ii;
			int jjj = jj;

			if (iii < ySize/2 - 1)
				iii++;

			if (jjj < xSize/2 - 1)
				jjj++;

			if (ii < 0 || ii >= ySize / 2 || jj < 0 || jj >= xSize / 2)
			{
				NEW_U_buff[i*xSize / 2 + j] = 0;
				NEW_V_buff[i*xSize / 2 + j] = 0;
			}
			else
			{
				NEW_U_buff[i*xSize / 2 + j] = (1 - a)*(1 - b)*U_buff[ii*xSize/2 + jj] +
					(1 - a)*b*U_buff[iii*xSize/2 + jj] +
					a*(1 - b)*U_buff[ii*xSize/2 + jjj] +
					a*b*U_buff[iii*xSize/2 + jjj];
				NEW_V_buff[i*xSize / 2 + j] = (1 - a)*(1 - b)*V_buff[ii*xSize/2 + jj] +
					(1 - a)*b*V_buff[iii*xSize/2 + jj] +
					a*(1 - b)*V_buff[ii*xSize/2 + jjj] +
					a*b*V_buff[iii*xSize/2 + jjj];
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