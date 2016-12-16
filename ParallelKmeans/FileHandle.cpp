#define _CRT_SECURE_NO_DEPRECATE
#include "FileHandle.h"
#include "Kmeans.h"
#include <stdio.h>
#include <stdlib.h>

#define HEADER_INFO 5
#define CIRCLE_INFO 4

static void readHeaderLine(FILE *fp, struct Input *input);
static void readCircles(FILE *fp, struct Input *input);

Input* readFile(char rout[])
{
	FILE *fp = NULL; 
	struct Input *input = NULL;
		
	if ((fp = fopen(rout, "r")) == NULL) {
		
		printf("Cannot open file.\n"); fflush(stdout);
		exit(EXIT_FAILURE);

	}
	else {
		input = (Input *) malloc(sizeof(Input));
		readHeaderLine(fp, input);
		readCircles(fp, input);
		fclose(fp);
	}

	return input;
}

static void readHeaderLine(FILE *fp, struct Input *input)
{
	int n; /*number of items scaned by fscanf*/
	
	if ((n = fscanf(fp, "%ld %ld %lf %lf %ld", &(*input).numCircles, &(*input).clusters, &(*input).deltaT, &(*input).interval, &(*input).maxItr)) != HEADER_INFO) {

		printf("File header does not contian the requested data.\n"); fflush(stdout);
		exit(EXIT_FAILURE);

	}
}

static void readCircles(FILE *fp, struct Input *input)
{
	int n; /*number of items scaned by fscanf*/
	int i;
	long tmp;

	(*input).r = (double *)malloc((*input).numCircles * sizeof(double));
	(*input).a = (double *)malloc((*input).numCircles * sizeof(double));
	(*input).b = (double *)malloc((*input).numCircles * sizeof(double));

	for (i = 0; i < (*input).numCircles; i++) {

		if ((n = fscanf(fp, "%ld %lf %lf %lf", &tmp, &((*input).a[i]), &((*input).b[i]), &((*input).r[i])) != CIRCLE_INFO)) {
			printf("Could not read circle - file text error.\n"); fflush(stdout);
			exit(EXIT_FAILURE);
		}
	}
}

void WriteToFile(char rout[], KmeansAns *ans, long numClusters)
{

	FILE *fp = NULL;

	if ((fp = fopen(rout, "w")) == NULL) {

		printf("Cannot open file.\n"); fflush(stdout);
		exit(EXIT_FAILURE);

	}

	fprintf(fp, "The Global minimum distance was found to be %lf at time step %lf.\nThe Centers found:\n", (*ans).minDistance, (*ans).timeStep);

	for (int i = 0; i < numClusters; i++) {
		fprintf(fp, "#%d: (x = %lf, y = %lf)\n", i, (*ans).CentersX[i], (*ans).CentersY[i]);
	}
	
	fclose(fp);
}




