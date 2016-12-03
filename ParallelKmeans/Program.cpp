#define _CRT_SECURE_NO_DEPRECATE
#include "FileHandle.h"
#include "PointsCalc.h"
#include "Kmeans.h"
#include <stdlib.h>
#include <stdio.h>
#include "omp.h"

#define PI 3.14159265358979323846

int main(int argc, char *argv[])
{
	long numSteps; /*Number of "world images" to run k-means on*/
	double step = 0;
	double theta;
	Point *points = NULL;
	KmeansAns **ans = NULL;
	long kmeansMindistIndex = 0;
	long j = 0;

	Input *input = readFile("C:\\Users\\ronylevari\\Documents\\Visual Studio 2015\\Projects\\\ParallelKmeans\\points.txt");

	numSteps = (long)((*input).interval / (*input).timeStep);

	ans = (KmeansAns **)malloc(numSteps * sizeof(KmeansAns*));

	do {
		printf("step: %lf\n", step);
		theta = (2 * PI * step / (*input).interval);
		points = calcPoints((*input).circles, (*input).numCircles, theta);
		ans[j] = runKmeans(points, (*input).numCircles, (*input).clusters, (*input).maxItr, step);
		step += (*input).timeStep;
		if (j == 0 || (*(ans[j])).minDistance < (*(ans[kmeansMindistIndex])).minDistance ) {
			kmeansMindistIndex = j;
		}
		j++;
	} while (numSteps > j);
	
	printf("*************************************\n");
	printf("General min dist: %lf\n", (*ans[kmeansMindistIndex]).minDistance);
	printf("Time : %lf\n", (*(ans[kmeansMindistIndex])).timeStep);
	printf("Centers\n");
	for (int i = 0; i < 2; i++) {
		printf("(x = %lf , y = %lf)\n", (*(ans[kmeansMindistIndex])).centers[i].x, (*(ans[kmeansMindistIndex])).centers[i].y);
	}
	printf("*************************************\n");
	free(points);
	free(ans);
	free(input);
}