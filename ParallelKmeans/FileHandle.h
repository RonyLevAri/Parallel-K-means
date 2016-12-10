#pragma once
#define readFile FileHandlereadFile

typedef struct Circle
{
	long index;
	double radius;
	double a;
	double b;
} Circle;

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
	Circle *circles;
}Input;

extern Input* readFile(char rout[]);



