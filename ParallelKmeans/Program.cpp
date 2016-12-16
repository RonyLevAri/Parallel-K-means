#define _CRT_SECURE_NO_DEPRECATE
#include "FileHandle.h"
#include "PointsCalc.h"
#include "Kmeans.h"
#include <stdlib.h>
#include <stdio.h>
#ifdef _OPENMP
	#include "omp.h"
#else
	#define omp_get_thread_num() 0
#endif
#include <mpi.h>

#define PI 3.14159265358979323846

#define INPUT_ROUTE "C:\\Users\\afeka.ACADEMIC\\Desktop\\ParallelKmeans\\kmean100.txt"
#define OUTPUT_ROUTE "C:\\Users\\afeka.ACADEMIC\\Desktop\\ParallelKmeans\\answer.txt"

static void allocateJobRange(Input *input, int world_size, double *startStep, double *endStep, long *numJobsToProc);
static void buildMpiKmeansAnsType(KmeansAns *ans, long numClusters, MPI_Datatype *mpiKmensAnsPtr);
static void buildMpiInputType(Input *input, MPI_Datatype *mpiKmensInputPtr);
static long findMinDistIndex(KmeansAns **ans, long arrSize);

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
	long kmeansMinDistIndex = 0; /*the minimal k-means distance of all calculated "world images"*/
	Input *input = NULL; /*data read from file*/
	double globalMinDist; /*the minimum k-means distance calculated by all nodes*/
	int isGlobalMinDistInMaster = 0; // 0 = false, 1 = true
	
	// Initialize the MPI environment
	master = tag = 0;
	ierr =  MPI_Init(&argc, &argv);
	if (ierr != MPI_SUCCESS) {
		printf("ERROR: colud not initialize MPI_Init\n"); fflush(stdout);
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	ierr = MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	if (ierr != MPI_SUCCESS) {
		printf("ERROR: colud not initialize MPI_Comm_size\n"); fflush(stdout);
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	if (ierr != MPI_SUCCESS) {
		printf("ERROR: colud not initialize MPI_Comm_rank\n"); fflush(stdout);
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (world_rank == master) {
		startTime = MPI_Wtime();
		input = readFile(INPUT_ROUTE); 
	} // only master reads data from file
	else {
		input = (Input *)malloc(sizeof(Input*));
	} // other processes initialize the input struct

	if (world_size > 1) {
		// bradcast common data to all processes only if there are more than 1 node
		MPI_Bcast(&((*input).numCircles), 1, MPI_LONG, master, MPI_COMM_WORLD);
		// all processes that are not the master should initialize their input circle data arrays before recieving them 
		if (world_rank != master) {
			(*input).r = (double *)malloc((*input).numCircles * sizeof(double));
			(*input).a = (double *)malloc((*input).numCircles * sizeof(double));
			(*input).b = (double *)malloc((*input).numCircles * sizeof(double));
		}
		// processes build derived datatype for more efficient process communication 
		buildMpiInputType(input, &mpiInput);
		MPI_Bcast(input, 1, mpiInput, master, MPI_COMM_WORLD);
	}
	
	// allocate deltaT's jobs (logic will work also with a single node)
	if (world_rank == master) {
		allocateJobRange(input, world_size, &startStep, &endStep, &numJobsToProc);
	}
	else { 
		MPI_Recv(&startStep, 1, MPI_DOUBLE, master, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&numJobsToProc, 1, MPI_LONG, master, tag, MPI_COMM_WORLD, &status);
	}

	// all nodes that were given jobs, participate in k-means calculations and communication
	if (numJobsToProc > 0) {
		
		ans = (KmeansAns **)malloc(numJobsToProc * sizeof(KmeansAns*));

		// procersses k-means loop can, after modification, be distributed among thrads using OpenMP
		
		theta = 0;

		#pragma omp parallel default(none) shared(ans, numJobsToProc, startStep, world_rank, input) private(jobCounter, theta, points) 
		{
			#pragma omp for 
			for (jobCounter = 0; jobCounter < numJobsToProc; jobCounter++) {
				theta = (2 * PI * (startStep + (jobCounter * (*input).deltaT)) / (*input).interval);
				points = calcPoints((*input).numCircles, theta, (*input).r, (*input).a, (*input).b);
				ans[jobCounter] = runKmeans(points, (*input).numCircles, (*input).clusters, (*input).maxItr, (startStep + (jobCounter * (*input).deltaT)));
			}
		} /*end of parallel construct - barrier*/

		//findMinDistIndex() was initilly processed during the k-means loop but had to be taken outside 
		// in order to parallel the loop (prevent loop iteration dependency) 
		kmeansMinDistIndex = findMinDistIndex(ans, numJobsToProc);

		MPI_Barrier(MPI_COMM_WORLD);

		printf("*************************************\n"); fflush(stdout);
		printf("#%d General min dist is %lf, the time is %lf\n", world_rank, (*ans[kmeansMinDistIndex]).minDistance, (*(ans[kmeansMinDistIndex])).timeStep); fflush(stdout);
		for (int i = 0; i < (*input).clusters; i++) {
			printf("Final Centers (x = %lf , y = %lf)\n", (*(ans[kmeansMinDistIndex])).CentersX[i], (*(ans[kmeansMinDistIndex])).CentersY[i]); fflush(stdout);
		}
		printf("*************************************\n"); fflush(stdout);
		
		// gather globalMinDist information (logic will work also with a single process)
		if (world_size > 1) {
			MPI_Allreduce(&((*(ans[kmeansMinDistIndex])).minDistance), &globalMinDist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
			buildMpiKmeansAnsType(ans[kmeansMinDistIndex], (*input).clusters, &mpiKmeansAns);
		} // activate Allreduce and build derived datatype only if more than one process is invovled
		else {
			globalMinDist = ((*(ans[kmeansMinDistIndex])).minDistance);
		} // if only one process invovled, set the global min distance to the distance found by the single process
		
		if (globalMinDist == ((*(ans[kmeansMinDistIndex])).minDistance)) {

			// if minimum found by none master process, send best "world image" to master
			if (world_rank != master) {
				MPI_Request request1;
				// use MPI_Isend in case more than one process have the nim distance, so they won't block 
				MPI_Isend(ans[kmeansMinDistIndex], 1, mpiKmeansAns, master, tag, MPI_COMM_WORLD, &request1);
			} // END OF SLAVE THREADS WORK
			else {
				isGlobalMinDistInMaster = 1;
			}
		}
	}

	// master should colect the min distance data for final result
	if (world_rank == master) {

		if (isGlobalMinDistInMaster == 0) {
			MPI_Recv(ans[kmeansMinDistIndex], 1, mpiKmeansAns, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
		}

		WriteToFile(OUTPUT_ROUTE, ans[kmeansMinDistIndex], (*input).clusters);
		finishTime = MPI_Wtime();
		printf("total time : %lf\n", finishTime - startTime);
	}

	free(points);
	free(ans);
	free(input);

	printf("process #%d bye bye\n", world_rank);
	// Finalize the MPI environment.
	ierr = MPI_Finalize();
}

static long findMinDistIndex(KmeansAns **ans, long arrSize)
{
	long minDistIndex = 0;
	for (int k = 1; k < arrSize; k++) {
		if ((*(ans[k])).minDistance < (*(ans[minDistIndex])).minDistance) {
			minDistIndex = k;
		}
	}
	return minDistIndex;
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