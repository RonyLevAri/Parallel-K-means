#define _CRT_SECURE_NO_DEPRECATE
#include "FileHandle.h"
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
		
		printf("Cannot open file.\n");
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
	
	if ((n = fscanf(fp, "%ld %ld %lf %lf %ld", &(*input).numCircles, &(*input).clusters, &(*input).timeStep, &(*input).interval, &(*input).maxItr)) != HEADER_INFO) {

		printf("File header does not contian the requested data.\n");
		exit(EXIT_FAILURE);

	}
}

static void readCircles(FILE *fp, struct Input *input)
{
	int n; /*number of items scaned by fscanf*/
	Circle *circles = NULL; /*The array of circles read from file*/
	int i;

	circles = (Circle *)malloc((*input).numCircles * sizeof(Circle));

	for (i = 0; i < (*input).numCircles; i++) {

		if ((n = fscanf(fp, "%ld %lf %lf %lf", &(circles[i].index), &(circles[i].a), &(circles[i].b), &(circles[i].radius))) != CIRCLE_INFO) {
			printf("Could not read circle - file text error.\n");
			exit(EXIT_FAILURE);
		}
		(*input).circles = circles;
	}
}





