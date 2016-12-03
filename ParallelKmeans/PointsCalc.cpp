#define _CRT_SECURE_NO_DEPRECATE
#include "PointsCalc.h"
#include "FileHandle.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


Point* calcPoints(Circle *circles, long numCircles, double theta)
{
	Point *points = NULL; /*The array of points calculated from circles*/
	int i;

	points = (Point *)malloc(numCircles * sizeof(Point));

	for (i = 0; i < numCircles; i++) {
		
		points[i].x = circles[i].a + circles[i].radius * cos(theta);
		points[i].y = circles[i].b + circles[i].radius * sin(theta);
		printf("point # %d: x = %lf y = %lf\n", i, points[i].x, points[i].y);
	}
	return points;
}