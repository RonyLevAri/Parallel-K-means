#define _CRT_SECURE_NO_DEPRECATE
#include "FileHandle.h"
#include "PointsCalc.h"
#include "Kmeans.h"
#include <stdlib.h>
#include <stdio.h>
#include "omp.h"
#include <mpi.h>

#define PI 3.14159265358979323846

int main(int argc, char *argv[])
{
	int ierr, world_size, world_rank, name_len, tag, master;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;
	long numJobsInTotal; /*Number of "world images" to run k-means on*/
	long numJobsPerProc;
	int jobsModulos;
	double startStep, endStep;
	double theta; /*"world image theta"*/
	Point *points = NULL; /*pionts, based on theta calculation*/
	KmeansAns **ans = NULL; /*"world image" k-means answer*/
	long kmeansMindistIndex = 0; /*the minimal k-means distance of all "world images"*/
	long j;
	int i;
	Input *input = NULL; /*data read from file*/
	double *r = NULL;
	double *a = NULL;
	double *b = NULL;
	double startTime, finishTime;
		
	// Initialize the MPI environment
	master = tag = 0;
	ierr =  MPI_Init(&argc, &argv);
	if (ierr != MPI_SUCCESS) {
		printf("ERROR: colud not initialize MPI_Init\n"); fflush(stdout);
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Get the number of processes
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	if (ierr != MPI_SUCCESS) {
		printf("ERROR: colud not initialize MPI_Comm_size\n"); fflush(stdout);
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Get the rank of the process
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	if (ierr != MPI_SUCCESS) {
		printf("ERROR: colud not initialize MPI_Comm_rank\n"); fflush(stdout);
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Get the name of the processor
	ierr =  MPI_Get_processor_name(processor_name, &name_len);
	if (ierr != MPI_SUCCESS) {
		printf("ERROR: colud not initialize MPI_Get_processor_name\n"); fflush(stdout);
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	printf("world_size = %d, rank = %d\n", world_size, world_rank); fflush(stdout);

	if (world_size < 2) {
		printf("ERROR: number of processes must be of size 2 at least\n"); fflush(stdout);
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (world_rank == master) {
		// record work start time
		startTime = MPI_Wtime();
		input = readFile("C:\\Users\\ronylevari\\Documents\\Visual Studio 2015\\Projects\\\ParallelKmeans\\points.txt"); fflush(stdout);
		//printf("read from file: numcircles = %ld, delta = %lf, clusters = %ld, maxIter = %ld\n", (*input).numCircles, (*input).deltaT, (*input).clusters, (*input).maxItr); fflush(stdout);
		//printf("read from file after population: numcircles = %ld, delta = %lf, clusters = %ld, maxIter = %ld\n", numCircles, deltaT, clusters, maxItr); fflush(stdout);
	} // only master reads data from file
	else {
		//printf("I am process #%d, initializing input\n", world_rank); fflush(stdout);
		input = (Input *)malloc(sizeof(Input*));
	} // other processes should initialize the struct

	// bradcast common data to all processes
	MPI_Bcast(&((*input).deltaT), 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
	MPI_Bcast(&((*input).numCircles), 1, MPI_LONG, master, MPI_COMM_WORLD);
	MPI_Bcast(&((*input).clusters), 1, MPI_LONG, master, MPI_COMM_WORLD);
	MPI_Bcast(&((*input).maxItr), 1, MPI_LONG, master, MPI_COMM_WORLD);
	MPI_Bcast(&((*input).interval), 1, MPI_DOUBLE, master, MPI_COMM_WORLD);

	if (world_rank != master) {
		(*input).r = (double *)malloc((*input).numCircles * sizeof(double));
		(*input).a = (double *)malloc((*input).numCircles * sizeof(double));
		(*input).b = (double *)malloc((*input).numCircles * sizeof(double));
	}

	MPI_Bcast(&((*input).r[0]), (*input).numCircles, MPI_DOUBLE, master, MPI_COMM_WORLD);
	MPI_Bcast(&((*input).a[0]), (*input).numCircles, MPI_DOUBLE, master, MPI_COMM_WORLD);
	MPI_Bcast(&((*input).b[0]), (*input).numCircles, MPI_DOUBLE, master, MPI_COMM_WORLD);

	if (world_rank == master) {
		
		// allocateJobs()
		numJobsInTotal = (long)((*input).interval / (*input).deltaT) + 1;
		jobsModulos = numJobsInTotal % world_size;
		printf("from process #%d, the number of jobs in total is: %ld, the number of processes is: %d and the modulos is: %ld\n", world_rank, numJobsInTotal, world_size, jobsModulos); fflush(stdout);

		startStep = 0;
		for (i = 1; i < world_size && startStep < (*input).interval; i++) {

			numJobsPerProc = numJobsInTotal / (long)world_size;

			if (i <= jobsModulos) {
				numJobsPerProc++;
			}

			endStep = startStep + (numJobsPerProc * (*input).deltaT - (*input).deltaT);
			printf("process #%d is going to do %ld jobs from %lf to %lf\n", i, numJobsPerProc, startStep, endStep); fflush(stdout);
			MPI_Send(&endStep, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MPI_Send(&startStep, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MPI_Send(&numJobsPerProc, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
			startStep = endStep + (*input).deltaT;
		}
		numJobsPerProc = numJobsInTotal / (long)world_size;
		endStep = (*input).interval;
		printf("process #%d will do %ld jobs from %lf to %lf\n", world_rank, numJobsPerProc, startStep, endStep); fflush(stdout);
	}
	else { // add condition for procesess < jobs - so that no nodes will be waiting
		MPI_Recv(&endStep, 1, MPI_DOUBLE, master, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&startStep, 1, MPI_DOUBLE, master, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&numJobsPerProc, 1, MPI_LONG, master, 0, MPI_COMM_WORLD, &status);
		printf("process #%d  will do %ld jobs from %lf to %lf\n", world_rank, numJobsPerProc, startStep, endStep); fflush(stdout);
	}

	if (numJobsPerProc > 0) {
		ans = (KmeansAns **)malloc(numJobsPerProc * sizeof(KmeansAns*));
		j = 0;
		do {
			//printf("#%d in step: %lf\n", world_rank, startStep); fflush(stdout);
			theta = (2 * PI * startStep / (*input).interval);
			points = calcPoints((*input).circles, (*input).numCircles, theta, (*input).r, (*input).a, (*input).b, world_rank);
			ans[j] = runKmeans(points, (*input).numCircles, (*input).clusters, (*input).maxItr, startStep, world_rank);
			startStep += (*input).deltaT;
			if (j == 0 || (*(ans[j])).minDistance < (*(ans[kmeansMindistIndex])).minDistance) {
				kmeansMindistIndex = j;
			}
			j++;
		} while (numJobsPerProc > j);

		printf("*************************************\n"); fflush(stdout);
		printf("#%d General min dist is %lf, the time is %lf\n", world_rank, (*ans[kmeansMindistIndex]).minDistance, (*(ans[kmeansMindistIndex])).timeStep); fflush(stdout);
		for (int i = 0; i < 3; i++) {
			printf("#%d Centers (x = %lf , y = %lf)\n", world_rank, (*(ans[kmeansMindistIndex])).centers[i].x, (*(ans[kmeansMindistIndex])).centers[i].y); fflush(stdout);
		}
		printf("*************************************\n"); fflush(stdout);

		//MPI_Datatype KmeansAnsType;
		//BuildKmeansAnsType(&((*(ans[kmeansMindistIndex])).minDistance), &((*(ans[kmeansMindistIndex])).timeStep), (*(ans[kmeansMindistIndex])).CentersX, (*(ans[kmeansMindistIndex])).CentersY, (*input).clusters, &KmeansAnsType);
		double globalMinDist;
		int globalMinDistInMaster = 0; // 0 = false, 1 = true

		MPI_Allreduce(&((*(ans[kmeansMindistIndex])).minDistance), &globalMinDist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		
		if (globalMinDist == ((*(ans[kmeansMindistIndex])).minDistance)) {

			// if minimum found by none master process , send best snapshoot to master
			if (world_rank != master) {
				MPI_Send(&((*(ans[kmeansMindistIndex])).timeStep), 1, MPI_DOUBLE, master, 0, MPI_COMM_WORLD);
				MPI_Send(&((*(ans[kmeansMindistIndex])).CentersX), (*input).clusters, MPI_DOUBLE, master, 0, MPI_COMM_WORLD);
				MPI_Send(&((*(ans[kmeansMindistIndex])).CentersY), (*input).clusters, MPI_DOUBLE, master, 0, MPI_COMM_WORLD);
			}
			else {
				globalMinDistInMaster = 1;
			}
		}

		if (world_rank == master) {
			
			if (globalMinDistInMaster == 0) {
				MPI_Recv(&((*(ans[kmeansMindistIndex])).timeStep), 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&((*(ans[kmeansMindistIndex])).CentersX), (*input).clusters, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&((*(ans[kmeansMindistIndex])).CentersY), (*input).clusters, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			}

			printf("*************************************\n"); fflush(stdout);
			printf("Final min dist is %lf, the time is %lf\n", (*ans[kmeansMindistIndex]).minDistance, (*(ans[kmeansMindistIndex])).timeStep); fflush(stdout);
			for (int i = 0; i < 3; i++) {
				printf("Final Centers (x = %lf , y = %lf)\n", (*(ans[kmeansMindistIndex])).centers[i].x, (*(ans[kmeansMindistIndex])).centers[i].y); fflush(stdout);
			}
			printf("*************************************\n"); fflush(stdout);
			finishTime = MPI_Wtime();
			printf("total time : %lf\n  ", finishTime - startTime);
		}

	}

	free(points);
	free(ans);
	free(input);
	free(a);
	free(b);
	free(r);

	// Finalize the MPI environment.
	ierr = MPI_Finalize();
}

void BuildKmeansAnsType(double *minDist, double *timeStep, double *centerX, double *centerY, long numClusters, MPI_Datatype *mpiKmensAnsPtr) {
	
	int blocklens[4] = {1, 1, numClusters, numClusters};
	MPI_Aint displacement[4]; // begining address of elements
	MPI_Datatype oldTypes[4] = {MPI_DOUBLE};
	MPI_Datatype newTypes[4];

	MPI_Aint startAddress;
	MPI_Aint address;

	displacement[0] = 0;

	MPI_Address(minDist, &startAddress);

	MPI_Address(timeStep, &address);
	displacement[1] = address - startAddress;

	MPI_Address(centerX, &address);
	displacement[2] = address - startAddress;

	MPI_Address(centerY, &address);
	displacement[3] = address - startAddress;

	MPI_Type_struct(4, blocklens, displacement, oldTypes, mpiKmensAnsPtr);
	MPI_Type_commit(mpiKmensAnsPtr);
}