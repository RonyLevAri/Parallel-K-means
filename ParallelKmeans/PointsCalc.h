#pragma once
#define calcPoints PointsCalcFunc

typedef struct Point
{
	double x;
	double y;
} Point;

extern Point* calcPoints(long numCircles, double theta, double *r, double *a, double *b, int rank);