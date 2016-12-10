#define _CRT_SECURE_NO_DEPRECATE
#include "PointsCalc.h"
#include "FileHandle.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


Point* calcPoints(Circle *circles, long numCircles, double theta, double *r, double *a, double *b, int rank)
{
	Point *points = NULL; /*The array of points calculated from circles*/
	int i;

	//printf("#%d numCircles = %ld, theta = %lf\n", rank, numCircles, theta); fflush(stdout);

	points = (Point *)malloc(numCircles * sizeof(Point));

	for (i = 0; i < numCircles; i++) {
		
		//points[i].x = circles[i].a + circles[i].radius * cos(theta);
		//points[i].y = circles[i].b + circles[i].radius * sin(theta);
		points[i].x = a[i] + r[i] * cos(theta);
		points[i].y = b[i] + r[i] * sin(theta);
		//printf("#%d a = %lf, b = %lf, r = %lf\n", rank, a[i], b[i], r[i]); fflush(stdout);
		//printf("#%d point # %d: x = %lf y = %lf\n", rank, i, points[i].x, points[i].y); fflush(stdout);
	}
	return points;
}