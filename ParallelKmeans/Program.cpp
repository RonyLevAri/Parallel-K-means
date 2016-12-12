#define _CRT_SECURE_NO_DEPRECATE
#include "FileHandle.h"
#include "PointsCalc.h"
#include "Kmeans.h"
#include <stdlib.h>
#include <stdio.h>
#include "omp.h"
#include <mpi.h>

#define PI 3.14159265358979323846

static void allocateJobRange(Input *input, int world_size, double *startStep, double *endStep, long *numJobsToProc);
static void buildMpiKmeansAnsType(KmeansAns *ans, long numClusters, MPI_Datatype *mpiKmensAnsPtr);
static void buildMpiInputType(Input *input, MPI_Datatype *mpiKmensInputPtr);

int main(int argc, char *argv[])
{
	// MPI variables
	int ierr, world_size, world_rank, tag, master;
	MPI_Status status;
	double startTime, finishTime;
	MPI_Datatype mpiInput, mpiKmeansAns; /*mpi datatype for a more effecient communication*/

	// local variables
	long numJobsToProc; /*number of "world images" to run k-means on*/
	int jobCounter;
	double startStep, endStep; /*start deltaT, end deltaT*/
	double theta; /*"world image" theta*/
	Point *points = NULL; /*pionts, based on theta calculation*/
	KmeansAns **ans = NULL; /*"world image" k-means answers array*/
	long kmeansMindistIndex = 0; /*the minimal k-means distance of all calculated "world images"*/
	Input *input = NULL; /*data read from file*/
	double globalMinDist; /*the minimum k-means distance calculated by all nodes*/
	int globalMinDistInMaster = 0; // 0 = false, 1 = true
	
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

	if (world_size < 2) {
		printf("ERROR: number of processes must be of size 2 at least\n"); fflush(stdout);
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (world_rank == master) {
		startTime = MPI_Wtime();
		input = readFile("C:\\Users\\ronylevari\\Documents\\Visual Studio 2015\\Projects\\\ParallelKmeans\\points.txt"); fflush(stdout);
	} // only master reads data from file
	else {
		input = (Input *)malloc(sizeof(Input*));
	} // other processes should initialize the input struct

	// bradcast common data to all processes
	MPI_Bcast(&((*input).numCircles), 1, MPI_LONG, master, MPI_COMM_WORLD);
	if (world_rank != master) {
		(*input).r = (double *)malloc((*input).numCircles * sizeof(double));
		(*input).a = (double *)malloc((*input).numCircles * sizeof(double));
		(*input).b = (double *)malloc((*input).numCircles * sizeof(double));
	}

	buildMpiInputType(input, &mpiInput);
	MPI_Bcast(input, 1, mpiInput, master, MPI_COMM_WORLD);

	//MPI_Bcast(&((*input).deltaT), 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
	//MPI_Bcast(&((*input).clusters), 1, MPI_LONG, master, MPI_COMM_WORLD);
	//MPI_Bcast(&((*input).maxItr), 1, MPI_LONG, master, MPI_COMM_WORLD);
	//MPI_Bcast(&((*input).interval), 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
	//MPI_Bcast(&((*input).r[0]), (*input).numCircles, MPI_DOUBLE, master, MPI_COMM_WORLD);
	//MPI_Bcast(&((*input).a[0]), (*input).numCircles, MPI_DOUBLE, master, MPI_COMM_WORLD);
	//MPI_Bcast(&((*input).b[0]), (*input).numCircles, MPI_DOUBLE, master, MPI_COMM_WORLD);

	// allocate deltaT's jobs
	if (world_rank == master) {
		allocateJobRange(input, world_size, &startStep, &endStep, &numJobsToProc);
	}
	else { 
		MPI_Recv(&endStep, 1, MPI_DOUBLE, master, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&startStep, 1, MPI_DOUBLE, master, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&numJobsToProc, 1, MPI_LONG, master, tag, MPI_COMM_WORLD, &status);
		printf("process #%d  will do %ld jobs from %lf to %lf\n", world_rank, numJobsToProc, startStep, endStep); fflush(stdout);
	}

	// do k-means calculations
	if (numJobsToProc > 0) {
		ans = (KmeansAns **)malloc(numJobsToProc * sizeof(KmeansAns*));
		jobCounter = 0;
		do {
			//printf("#%d in step: %lf\n", world_rank, startStep); fflush(stdout);
			theta = (2 * PI * startStep / (*input).interval);
			points = calcPoints((*input).circles, (*input).numCircles, theta, (*input).r, (*input).a, (*input).b, world_rank);
			ans[jobCounter] = runKmeans(points, (*input).numCircles, (*input).clusters, (*input).maxItr, startStep, world_rank);
			startStep += (*input).deltaT;
			if (jobCounter == 0 || (*(ans[jobCounter])).minDistance < (*(ans[kmeansMindistIndex])).minDistance) {
				kmeansMindistIndex = jobCounter;
			}
			jobCounter++;
		} while (numJobsToProc > jobCounter);

		printf("*************************************\n"); fflush(stdout);
		printf("#%d General min dist is %lf, the time is %lf\n", world_rank, (*ans[kmeansMindistIndex]).minDistance, (*(ans[kmeansMindistIndex])).timeStep); fflush(stdout);
		for (int i = 0; i < 3; i++) {
			printf("#%d Centers (x = %lf , y = %lf)\n", world_rank, (*(ans[kmeansMindistIndex])).centers[i].x, (*(ans[kmeansMindistIndex])).centers[i].y); fflush(stdout);
		}
		printf("*************************************\n"); fflush(stdout);
		
		// gather globalMinDist 
		MPI_Allreduce(&((*(ans[kmeansMindistIndex])).minDistance), &globalMinDist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		buildMpiKmeansAnsType(ans[kmeansMindistIndex], (*input).clusters, &mpiKmeansAns);
		
		if (globalMinDist == ((*(ans[kmeansMindistIndex])).minDistance)) {

			// if minimum found by none master process, send best "world image" to master
			if (world_rank != master) {
				MPI_Request request1, request2, request3;
				MPI_Isend(ans[kmeansMindistIndex], 1, mpiKmeansAns, master, tag, MPI_COMM_WORLD, &request1);
			}
			else {
				globalMinDistInMaster = 1;
			}
		}

		
		if (world_rank == master) {

			if (globalMinDistInMaster == 0) {
				printf("master recieving\n");
				MPI_Recv(ans[kmeansMindistIndex], 1, mpiKmeansAns, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
			}

			printf("*************************************\n"); fflush(stdout);
			printf("Final min dist is %lf, the time is %lf\n", globalMinDist, (*(ans[kmeansMindistIndex])).timeStep); fflush(stdout);
			for (int i = 0; i < (*input).clusters; i++) {
				printf("Final Centers (x = %lf , y = %lf)\n", (*(ans[kmeansMindistIndex])).CentersX[i], (*(ans[kmeansMindistIndex])).CentersY[i]); fflush(stdout);
			}
			printf("*************************************\n"); fflush(stdout);
			finishTime = MPI_Wtime();
			printf("total time : %lf\n  ", finishTime - startTime);
		}
	}

	free(points);
	free(ans);
	free(input);

	printf("process #%d bye bye\n", world_rank);
	// Finalize the MPI environment.
	ierr = MPI_Finalize();
}

static void allocateJobRange(Input *input, int world_size, double *startStep, double *endStep, long *numJobsToProc)
{
	long numJobsInTotal; /*Total number of "world images" to run k-means on*/
	int jobsModulos; /*Modulos of number of "world images" to run k-means on*/

	numJobsInTotal = (long)((*input).interval / (*input).deltaT) + 1;
	jobsModulos = numJobsInTotal % world_size;
	printf("from master, the number of jobs in total is: %ld, the number of processes is: %d and the modulos is: %ld\n", numJobsInTotal, world_size, jobsModulos); fflush(stdout);

	(*startStep) = 0;
	for (int i = 1; i < world_size && (*startStep) < (*input).interval; i++) {

		(*numJobsToProc) = numJobsInTotal / (long)world_size;

		if (i <= jobsModulos) {
			(*numJobsToProc)++;
		}

		(*endStep) = (*startStep) + ((*numJobsToProc) * (*input).deltaT - (*input).deltaT);
		printf("process #%d is going to do %ld jobs from %lf to %lf\n", i, (*numJobsToProc), (*startStep), (*endStep)); fflush(stdout);
		MPI_Send(&(*endStep), 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		MPI_Send(&(*startStep), 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		MPI_Send(&(*numJobsToProc), 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
		(*startStep) = (*endStep) + (*input).deltaT;
	}
	(*numJobsToProc) = numJobsInTotal / (long)world_size;
	(*endStep) = (*input).interval;
	printf("master will do %ld jobs from %lf to %lf\n", (*numJobsToProc), (*startStep), (*endStep)); fflush(stdout);
}

static void buildMpiInputType(Input *input, MPI_Datatype *mpiKmensInputPtr)
{
	int blocklens[8] = {1, 1, 1, 1, 1, (*input).numCircles, (*input).numCircles, (*input).numCircles};
	MPI_Aint displacement[8]; // begining address of elements
	MPI_Datatype oldTypes[8] = {MPI_LONG, MPI_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

	MPI_Aint startAddress;
	MPI_Aint address;

	displacement[0] = 0;

	MPI_Address(&((*input).numCircles), &startAddress);

	MPI_Address(&((*input).clusters), &address);
	displacement[1] = address - startAddress;

	MPI_Address(&((*input).deltaT), &address);
	displacement[2] = address - startAddress;

	MPI_Address(&((*input).interval), &address);
	displacement[3] = address - startAddress;

	MPI_Address(&((*input).maxItr), &address);
	displacement[4] = address - startAddress;

	MPI_Address(&((*input).r[0]), &address);
	displacement[5] = address - startAddress;

	MPI_Address(&((*input).a[0]), &address);
	displacement[6] = address - startAddress;

	MPI_Address(&((*input).b[0]), &address);
	displacement[7] = address - startAddress;

	MPI_Type_struct(8, blocklens, displacement, oldTypes, mpiKmensInputPtr);
	MPI_Type_commit(mpiKmensInputPtr);
}

static void buildMpiKmeansAnsType(KmeansAns *ans, long numClusters, MPI_Datatype *mpiKmensAnsPtr)
{
	
	int blocklens[4] = {1, 1, numClusters, numClusters};
	MPI_Aint displacement[4]; // begining address of elements
	MPI_Datatype oldTypes[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE , MPI_DOUBLE};

	MPI_Aint startAddress;
	MPI_Aint address;

	displacement[0] = 0;

	MPI_Address(&((*ans).minDistance), &startAddress);

	MPI_Address(&((*ans).timeStep), &address);
	displacement[1] = address - startAddress;

	MPI_Address(&((*ans).CentersX[0]), &address);
	displacement[2] = address - startAddress;

	MPI_Address(&((*ans).CentersY[0]), &address);
	displacement[3] = address - startAddress;

	MPI_Type_struct(4, blocklens, displacement, oldTypes, mpiKmensAnsPtr);
	MPI_Type_commit(mpiKmensAnsPtr);
}