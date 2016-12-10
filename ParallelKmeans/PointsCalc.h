#pragma once
#define calcPoints PointsCalcFunc
#include "FileHandle.h"

typedef struct Point
{
	double x;
	double y;
} Point;

extern Point* calcPoints(Circle *circles, long numCircles, double theta, double *r, double *a, double *b, int rank);