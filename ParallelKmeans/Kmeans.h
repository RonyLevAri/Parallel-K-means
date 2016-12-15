#pragma once
#define runKmeans KmeansMainFunc
#include "PointsCalc.h"

typedef struct Cluster
{
	long index;
	Point center;
	long numClustPoints;
	Point *clustPoints;
}Cluster;

typedef struct KmeansAns
{
	double minDistance;
	double timeStep;
	double *CentersX;
	double *CentersY;
}KmeansAns;

extern KmeansAns* runKmeans(Point *points, long numPoints, long numClusters, long maxIter, double step);