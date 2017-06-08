
#include "ImageProcessing.h"
#include "ImageInterpolation.h"

#include <QDebug>

int upto(int x, int y) {
	return (x + y - 1) & ~(y - 1);
}

void imageProcessingFun(const QString& progName, QImage* const outImgs, const QImage* const inImgs, const QVector<double>& params) 
{
	int X_SIZE = inImgs->width();
	int Y_SIZE = inImgs->height();

	/* NOTE: Calculate output image resolution and construct output image object */

	if(progName == "Sample and hold") 
	{	
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */

		/* Create empty output image */

		int NEW_X_SIZE = upto(X_SIZE * params[1], 4);
		int NEW_Y_SIZE = upto(X_SIZE * params[0], 4);

		*outImgs = *(new QImage(NEW_X_SIZE, NEW_Y_SIZE, inImgs->format()));

		sampleAndHold(inImgs->bits(),X_SIZE,Y_SIZE,outImgs->bits(),NEW_X_SIZE,NEW_Y_SIZE);
	}
	else if (progName == "Bilinear") 
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */

		/* Create empty output image */

		int NEW_X_SIZE = upto(X_SIZE * params[1], 4);
		int NEW_Y_SIZE = upto(X_SIZE * params[0], 4);

		*outImgs = *(new QImage(NEW_X_SIZE, NEW_Y_SIZE, inImgs->format()));

		bilinearInterpolate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), NEW_X_SIZE, NEW_Y_SIZE);
	}
	else if (progName == "Bicubic")
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */

		int NEW_X_SIZE = upto(X_SIZE * params[1], 4);
		int NEW_Y_SIZE = upto(X_SIZE * params[0], 4);

		*outImgs = *(new QImage(NEW_X_SIZE, NEW_Y_SIZE, inImgs->format()));

		bicubicInterpolate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), NEW_X_SIZE, NEW_Y_SIZE);

	}
	else if(progName == "Rotation") 
	{	
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Rotation angle in degrees is params[0]*/
		/* Center of rotation coordinates are (XSIZE/2, YSIZE/2) */

		int m = X_SIZE/2;
		int n = Y_SIZE/2;

		*outImgs = *(new QImage(X_SIZE, Y_SIZE, inImgs->format()));

		imageRotate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(),m,n,params[0]);
	
	}
	else if (progName == "Rotation Bilinear") 
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Rotation angle in degrees is params[0]*/
		/* Center of rotation coordinates are (XSIZE/2, YSIZE/2) */

		int m = X_SIZE / 2;
		int n = Y_SIZE / 2;

		*outImgs = *(new QImage(X_SIZE, Y_SIZE, inImgs->format()));

		imageRotateBilinear(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), m, n, params[0]);
	}

}

