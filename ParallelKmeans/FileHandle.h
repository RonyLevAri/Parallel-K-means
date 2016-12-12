#pragma once
#define readFile FileHandleReadFile
#define WriteToFile FileHandleWriteToFile
#include "Kmeans.h"

typedef struct Input
{
	long numCircles;
	long clusters;
	double deltaT;
	double interval;
	long maxItr;
	double *r;
	double *a;
	double *b;
}Input;

extern Input* readFile(char rout[]);

extern void WriteToFile(char rout[], KmeansAns *ans, long numClusters);



