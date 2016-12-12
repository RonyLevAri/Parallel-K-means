#define _CRT_SECURE_NO_DEPRECATE
#include "PointsCalc.h"
#include "FileHandle.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


Point* calcPoints(long numCircles, double theta, double *r, double *a, double *b, int rank)
{
	Point *points = NULL; /*The array of points calculated from circles*/
	int i;

	points = (Point *)malloc(numCircles * sizeof(Point));

	for (i = 0; i < numCircles; i++) {
		
		points[i].x = a[i] + r[i] * cos(theta);
		points[i].y = b[i] + r[i] * sin(theta);
	}
	return points;
}