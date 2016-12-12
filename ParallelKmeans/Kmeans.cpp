#define _CRT_SECURE_NO_DEPRECATE
#include "PointsCalc.h"
#include "Kmeans.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


static void initiateClusters(Point *points, Cluster *clusters, long numPoints, long numClusters);
static void resetClusters(Cluster *clusters, long numClusters);

KmeansAns* runKmeans(Point *points, long numPoints, long numClusters, long maxIter, double step, int rank)
{
	KmeansAns *ans = (KmeansAns *)malloc(sizeof(KmeansAns));
	long clusterIndex;
	double newCenterX = 0, newCenterY = 0;
	long i, j;
	bool changed = true;
	double generalMinDist = -1;
	long iter = 0;
	Cluster *clusters = (Cluster *)malloc(numClusters * sizeof(Cluster));

	initiateClusters(points, clusters, numPoints, numClusters);

	while (changed && iter < maxIter) {

		//printf("#%d, Iter #%ld\n", rank, iter);

		if (iter != 0) {
			resetClusters(clusters, numClusters);
		}


		generalMinDist = -1;

		// calcDistances();
		for (i = 0; i < numPoints; i++) {
			//printf("#%d point #%ld: x = %lf y = %lf\n", rank, i, points[i].x, points[i].y); fflush(stdout);
			double d = 0, minD = -1;
			for (j = 0; j < numClusters; j++) {
				//printf("#%d cluster #%ld with center: x = %lf y = %lf\n", rank, j, clusters[j].center.x, clusters[j].center.y); fflush(stdout);
				d = fabs(sqrt(pow(points[i].x - clusters[j].center.x, 2) + pow(points[i].y - clusters[j].center.y, 2)));
				//printf("#%d calculated distance from center %ld (%lf, %lf) = %lf calculated numD = %lf\n", rank, j, clusters[j].center.x, clusters[j].center.y, d, minD); fflush(stdout);
				if (minD == -1 || d < minD) {
					minD = d;
					clusterIndex = j;
				}
			}
			//printf("#%d Point #%ld is in center #%ld\n", rank, i, clusterIndex); fflush(stdout);
			//assign point to cluster
			clusters[clusterIndex].clustPoints[clusters[clusterIndex].numClustPoints] = points[i];
			clusters[clusterIndex].numClustPoints++;
			//printf("#%d Cluster #%ld has %ld points\n", rank, clusterIndex, clusters[clusterIndex].numClustPoints); fflush(stdout);
			if ((generalMinDist == -1) || (minD < generalMinDist)) {
				generalMinDist = minD;
			}
			//printf("#%d General min distance is %lf\n", rank, generalMinDist); fflush(stdout);
		}

		// calcCenters();
		changed = false;
		for (i = 0; i < numClusters; i++) {
			newCenterX = 0;
			newCenterY = 0;
			//printf("#%d cluster #%ld with center: x = %lf y = %lf\n", rank, i, clusters[i].center.x, clusters[i].center.y); fflush(stdout);
			for (j = 0; j < clusters[i].numClustPoints; j++) {
				newCenterX += clusters[i].clustPoints[j].x;
				newCenterY += clusters[i].clustPoints[j].y;
			}
			newCenterX = newCenterX / clusters[i].numClustPoints;
			newCenterY = newCenterY / clusters[i].numClustPoints;
			//printf("#%d cluster #%ld new center: x = %lf y = %lf\n", rank, i, newCenterX, newCenterY); fflush(stdout);
			if ((clusters[i].center.x != newCenterX) && (clusters[i].center.y != newCenterY)) {
				changed = true;
				clusters[i].center.x = newCenterX;
				clusters[i].center.y = newCenterY;
			}
		}
		iter++;
	}
	// Kmeans Answer
	(*ans).minDistance = generalMinDist;
	(*ans).timeStep = step;
	(*ans).CentersX = (double*)malloc(numClusters * sizeof(double));
	(*ans).CentersY = (double*)malloc(numClusters * sizeof(double));
	for (i = 0; i < numClusters; i++) {
		(*ans).CentersX[i] = clusters[i].center.x;
		(*ans).CentersY[i] = clusters[i].center.y;
	}

	return ans;
}

static void initiateClusters(Point *points, Cluster *clusters, long numPoints, long numClusters)
{
	long i = 0;

	for (i = 0; i < numClusters; i++) {
		clusters[i].center = points[i];
		clusters[i].index = i;
		clusters[i].numClustPoints = 0;
		clusters[i].clustPoints = (Point *)malloc(numPoints * sizeof(Point));
	}
}

static void resetClusters(Cluster *clusters, long numClusters)
{
	long i = 0;

	for (i = 0; i < numClusters; i++) {
		clusters[i].numClustPoints = 0;
	}
}